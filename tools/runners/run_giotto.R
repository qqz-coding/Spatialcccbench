library(Giotto)
library(data.table)
library(optparse)

option_list <- list(
  make_option("--analysis_dataset", type = "character"),
  make_option("--count_path", type = "character"),
  make_option("--meta_path", type = "character"),
  make_option("--LR_ref_path", type = "character"),
  make_option("--output_root", type = "character", default = "result"),
  make_option("--python_path", type = "character", default = ""),
  make_option("--n_perms", type = "integer", default = 1000),
  make_option("--n_cores", type = "integer", default = 2)
)
opt <- parse_args(OptionParser(option_list = option_list))
required <- c("analysis_dataset", "count_path", "meta_path", "LR_ref_path")
missing <- required[vapply(required, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
if (length(missing) > 0) stop("Missing required options: ", paste(missing, collapse = ", "))

output_dir <- file.path(opt$output_root, opt$analysis_dataset, "giotto")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
instruction_args <- list(save_dir = output_dir, save_plot = FALSE, show_plot = FALSE)
if (nzchar(opt$python_path)) instruction_args$python_path <- opt$python_path
myinstructions <- do.call(createGiottoInstructions, instruction_args)

# CSV contract: rows are spots, columns are genes; Giotto receives genes by spots.
count_df <- read.csv(opt$count_path, row.names = 1, check.names = FALSE)
meta_df <- read.csv(opt$meta_path, row.names = 1, check.names = FALSE)
needed_meta <- c("spatial1", "spatial2", "cell_type")
if (!all(needed_meta %in% colnames(meta_df))) stop("Metadata must contain spatial1, spatial2, and cell_type.")
common_spots <- intersect(rownames(count_df), rownames(meta_df))
if (length(common_spots) == 0) stop("Count and metadata CSV files have no shared spot IDs.")
count_df <- count_df[common_spots, , drop = FALSE]
meta_df <- meta_df[common_spots, , drop = FALSE]
st_count_df <- t(as.matrix(count_df))
st_pos_df <- meta_df[, c("spatial1", "spatial2"), drop = FALSE]
colnames(st_pos_df) <- c("row", "col")
st_meta_df <- data.frame(celltypes = as.character(meta_df$cell_type), row.names = common_spots)

st.obj <- createGiottoObject(raw_exprs = st_count_df,
                             spatial_locs = st_pos_df,
                             cell_metadata = st_meta_df,
                             instructions = myinstructions)

# preprocessing
st.obj <- filterGiotto(gobject = st.obj,
                       expression_threshold = 1,
                       gene_det_in_min_cells = 1,
                       min_det_genes_per_cell = 100,
                       expression_values = c('raw'),
                       verbose = T)
st.obj <- normalizeGiotto(gobject = st.obj)

st.obj <- calculateHVG(gobject = st.obj)
st.obj <- runPCA(gobject = st.obj)
st.obj <- runUMAP(st.obj, dimensions_to_use = 1:15)


st.obj = createSpatialNetwork(gobject = st.obj)

# identify genes with a spatial coherent expression profile
km_spatialgenes = binSpect(st.obj, bin_method = 'kmeans')

# codes below is just for proximity
# cell_proximities = cellProximityEnrichment(gobject = st.obj,
#                                            cluster_column = 'celltypes',
#                                            spatial_network_name = 'Delaunay_network',
#                                            adjust_method = 'fdr',
#                                            number_of_simulations = 1000)


# ligand & receptor database can be changed
LR_data = data.table::fread(opt$LR_ref_path)
if (!all(c("ligand.symbol", "receptor.symbol") %in% colnames(LR_data))) {
  stop("LR reference must contain ligand.symbol and receptor.symbol.")
}
LR_data_1 = LR_data

LR_data[, ligand_det := ifelse(ligand.symbol %in% st.obj@gene_ID, T, F)]
LR_data[, receptor_det := ifelse(receptor.symbol %in% st.obj@gene_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$ligand.symbol
select_receptors = LR_data_det$receptor.symbol


## get statistical significance of gene pair expression changes based on expression ##
# expr_only_scores = exprCellCellcom(gobject = st.obj,
#                                    cluster_column = 'celltypes',
#                                    random_iter = 1000,
#                                    gene_set_1 = select_ligands,
#                                    gene_set_2 = select_receptors,
#                                    verbose = F)

giotto_lapply = function(X, cores = NA, fun, ...) {

  # get type of os
  os = .Platform$OS.type

  # set number of cores automatically, but with limit of 10
  if(is.na(cores) | !is.numeric(cores)) {

    cores = parallel::detectCores()
    if(cores <= 2) {
      cores = cores
    } else {
      cores = cores - 2
      cores = ifelse(cores > 10, 10, cores)
    }
  }

  if(os == 'unix') {
    save_list = parallel::mclapply(X = X, mc.cores = cores,
                                   FUN = fun, ...)
  } else if(os == 'windows') {
    save_list = parallel::mclapply(X = X, mc.cores = 1,
                                   FUN = fun, ...)

    # !! unexplainable errors are returned for some nodes !! #
    # currently disabled #
    #cl <- parallel::makeCluster(cores)
    #save_list = parallel::parLapply(cl = cl, X = X,
    #                                fun = fun, ...)
  }

  return(save_list)
}

spatCellCellcom = function(gobject,
                           spatial_network_name = 'Delaunay_network',
                           cluster_column = 'cell_types',
                           random_iter = 1000,
                           gene_set_1,
                           gene_set_2,
                           log2FC_addendum = 0.1,
                           min_observations = 2,
                           detailed = FALSE,
                           adjust_method = c("fdr", "bonferroni","BH", "holm", "hochberg", "hommel",
                                             "BY", "none"),
                           adjust_target = c('genes', 'cells'),
                           do_parallel = TRUE,
                           cores = NA,
                           set_seed = TRUE,
                           seed_number = 1234,
                           verbose = c('a little', 'a lot', 'none')) {

  verbose = match.arg(verbose, choices = c('a little', 'a lot', 'none'))

  ## check if spatial network exists ##
  spat_networks = showNetworks(gobject = gobject, verbose = F)
  if(!spatial_network_name %in% spat_networks) {
    stop(spatial_network_name, ' is not an existing spatial network \n',
         'use showNetworks() to see the available networks \n',
         'or create a new spatial network with createSpatialNetwork() \n')
  }


  cell_metadata = pDataDT(gobject)

  ## get all combinations between cell types
  all_uniq_values = unique(cell_metadata[[cluster_column]])
  same_DT = data.table(V1 = all_uniq_values, V2 = all_uniq_values)
  combn_DT = as.data.table(t(combn(all_uniq_values, m = 2)))
  # combn_DT = rbind(same_DT, combn_DT)

  ## parallel option ##
  if(do_parallel == TRUE) {


    savelist = giotto_lapply(X = 1:nrow(combn_DT), cores = cores, fun = function(row) {

      cell_type_1 = combn_DT[row][['V1']]
      cell_type_2 = combn_DT[row][['V2']]

      specific_scores = specificCellCellcommunicationScores(gobject = gobject,
                                                            cluster_column = cluster_column,
                                                            random_iter = random_iter,
                                                            cell_type_1 = cell_type_1,
                                                            cell_type_2 = cell_type_2,
                                                            gene_set_1 = gene_set_1,
                                                            gene_set_2 = gene_set_2,
                                                            spatial_network_name = spatial_network_name,
                                                            log2FC_addendum = log2FC_addendum,
                                                            min_observations = min_observations,
                                                            detailed = detailed,
                                                            adjust_method = adjust_method,
                                                            adjust_target = adjust_target,
                                                            set_seed = set_seed,
                                                            seed_number = seed_number)

    })


  } else {

    ## for loop over all combinations ##
    savelist = list()
    countdown = nrow(combn_DT)

    for(row in 1:nrow(combn_DT)) {

      cell_type_1 = combn_DT[row][['V1']]
      cell_type_2 = combn_DT[row][['V2']]

      if(verbose == 'a little' | verbose == 'a lot') cat('\n\n PROCESS nr ', countdown,': ', cell_type_1, ' and ', cell_type_2, '\n\n')

      if(verbose %in% c('a little', 'none')) {
        specific_verbose = F
      } else {
        specific_verbose = T
      }

      specific_scores = specificCellCellcommunicationScores(gobject = gobject,
                                                            cluster_column = cluster_column,
                                                            random_iter = random_iter,
                                                            cell_type_1 = cell_type_1,
                                                            cell_type_2 = cell_type_2,
                                                            gene_set_1 = gene_set_1,
                                                            gene_set_2 = gene_set_2,
                                                            spatial_network_name = spatial_network_name,
                                                            log2FC_addendum = log2FC_addendum,
                                                            min_observations = min_observations,
                                                            detailed = detailed,
                                                            adjust_method = adjust_method,
                                                            adjust_target = adjust_target,
                                                            set_seed = set_seed,
                                                            seed_number = seed_number,
                                                            verbose = specific_verbose)
      savelist[[row]] = specific_scores
      countdown = countdown - 1
    }

  }

  finalDT = do.call('rbind', savelist)

  # data.table variables
  LR_comb = LR_expr = NULL

  data.table::setorder(finalDT, LR_comb, -LR_expr)

  return(finalDT)
}

start_time = Sys.time()
## get statistical significance of gene pair expression changes upon cell-cell interaction
spatial_all_scores = spatCellCellcom(gobject = st.obj,
                                     spatial_network_name = 'Delaunay_network',
                                     cluster_column = 'celltypes',
                                     random_iter = opt$n_perms,
                                     gene_set_1 = select_ligands,
                                     gene_set_2 = select_receptors,
                                     adjust_method = 'fdr',
                                     do_parallel = T,
                                     cores = opt$n_cores,
                                     verbose = 'a little')
end_time = Sys.time()
total_time =  difftime(end_time, start_time, units = "secs")

selected_spat = spatial_all_scores[pvalue<0.05 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
#selected_spat = read.csv("D:/benchmark/test/spatial_scores_lymphnode.tsv",sep="\t")

output_path <- file.path(output_dir, "result.csv")
write.csv(selected_spat, file = output_path, row.names = FALSE)
cat(sprintf("runtime_seconds: %.3f\n", as.numeric(total_time)))
cat("output:", normalizePath(output_path), "\n")
