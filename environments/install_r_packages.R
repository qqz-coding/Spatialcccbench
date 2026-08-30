args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[[1]] %in% c("spatalk", "giotto")) {
  stop("Usage: Rscript install_r_packages.R <spatalk|giotto>")
}

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

if (args[[1]] == "spatalk") {
  if (!requireNamespace("NNLM", quietly = TRUE)) {
    remotes::install_github("linxihui/NNLM", upgrade = "never")
  }
  if (!requireNamespace("SpaTalk", quietly = TRUE)) {
    remotes::install_github("ZJUFanLab/SpaTalk", upgrade = "never")
  }
  cat("SpaTalk version:", as.character(utils::packageVersion("SpaTalk")), "\n")
}

if (args[[1]] == "giotto") {
  wanted <- "3.3.2"
  installed <- if (requireNamespace("Giotto", quietly = TRUE)) {
    as.character(utils::packageVersion("Giotto"))
  } else {
    ""
  }
  if (installed != wanted) {
    remotes::install_github("giotto-suite/Giotto@v3.3.2", upgrade = "never")
  }
  cat("Giotto version:", as.character(utils::packageVersion("Giotto")), "\n")
}
