.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("GRAB", quietly = TRUE)) {
    message("Installing GRAB package from GitHub...")
    devtools::install_github("GeneticAnalysisinBiobanks/GRAB",ref="main")
  }
}
