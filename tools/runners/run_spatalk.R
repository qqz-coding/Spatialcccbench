suppressPackageStartupMessages({
  library(optparse)
  library(tidyr)
  library(SpaTalk)
})

option_list <- list(
  make_option("--analysis_dataset", type = "character"),
  make_option("--count_path", type = "character"),
  make_option("--meta_path", type = "character"),
  make_option("--LR_ref_path", type = "character"),
  make_option("--output_root", type = "character", default = "result"),
  make_option("--species", type = "character", default = "Human"),
  make_option("--spot_max_cell", type = "integer", default = 30)
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("analysis_dataset", "count_path", "meta_path", "LR_ref_path")
missing <- required[vapply(required, function(x) is.null(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  stop("Missing required options: ", paste(missing, collapse = ", "))
}

output_dir <- file.path(opt$output_root, opt$analysis_dataset, "spatalk")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# CSV contract: rows are spots, columns are genes; the first column is the spot ID.
count_df <- read.csv(opt$count_path, row.names = 1, check.names = FALSE)
meta_df <- read.csv(opt$meta_path, row.names = 1, check.names = FALSE)
needed_meta <- c("spatial1", "spatial2", "cell_type")
if (!all(needed_meta %in% colnames(meta_df))) {
  stop("Metadata must contain: ", paste(needed_meta, collapse = ", "))
}
common_spots <- intersect(rownames(count_df), rownames(meta_df))
if (length(common_spots) == 0) stop("Count and metadata CSV files have no shared spot IDs.")
count_df <- count_df[common_spots, , drop = FALSE]
meta_df <- meta_df[common_spots, , drop = FALSE]

spot_meta <- data.frame(
  X = common_spots,
  cell = common_spots,
  x = as.integer(meta_df$spatial1),
  y = as.integer(meta_df$spatial2),
  celltype = as.character(meta_df$cell_type),
  spot = common_spots,
  row.names = common_spots,
  check.names = FALSE
)

x_span <- max(spot_meta$x) - min(spot_meta$x)
y_span <- max(spot_meta$y) - min(spot_meta$y)
spot_data <- generate_spot(
  st_data = as.matrix(count_df),
  st_meta = spot_meta,
  x_min = min(spot_meta$x),
  x_res = max(1, ceiling(x_span / 100)),
  x_max = max(spot_meta$x),
  y_min = min(spot_meta$y),
  y_res = max(1, ceiling(y_span / 100)),
  y_max = max(spot_meta$y)
)

obj <- createSpaTalk(
  st_data = as.matrix(spot_data$st_data),
  st_meta = spot_data$st_meta[, c(1, 2, 3)],
  species = opt$species,
  if_st_is_sc = FALSE,
  spot_max_cell = opt$spot_max_cell
)
obj <- dec_celltype(
  object = obj,
  sc_data = as.matrix(count_df),
  sc_celltype = spot_meta$celltype
)

lr_data <- read.csv(opt$LR_ref_path, check.names = FALSE)
required_lr <- c("ligand.symbol", "receptor.symbol")
if (!all(required_lr %in% colnames(lr_data))) {
  stop("LR reference must contain ligand.symbol and receptor.symbol.")
}
lr_data <- lr_data[, required_lr]
colnames(lr_data) <- c("ligand", "receptor")
lr_data <- separate_rows(lr_data, receptor, sep = ",\\s*")
lr_data$species <- opt$species

started <- Sys.time()
obj <- find_lr_path(object = obj, lrpairs = lr_data, pathways = pathways)
obj <- dec_cci_all(object = obj)
runtime <- as.numeric(difftime(Sys.time(), started, units = "secs"))

output_path <- file.path(output_dir, "result.csv")
write.csv(obj@lrpair, file = output_path)
cat(sprintf("runtime_seconds: %.3f\n", runtime))
cat("output:", normalizePath(output_path), "\n")
