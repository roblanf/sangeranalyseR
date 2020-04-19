.onAttach <- function(libname, pkgname) {
  log_layout(layout_glue_colors)
  log_threshold(TRACE)
  packageStartupMessage("Welcome to sangeranalyseR")
}
