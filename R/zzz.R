#' @useDynLib DISQ, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.onAttach <- function(libname, pkgname) {

  # Avoid attaching extra packages during devtools load/document steps
  # (keeps devtools::document() stable and fast)
  if (nzchar(Sys.getenv("DEVTOOLS_LOAD"))) {
    return(invisible(NULL))
  }

  attach_if_installed <- function(pkg) {
    if (paste0("package:", pkg) %in% search()) return(invisible(TRUE))

    if (!requireNamespace(pkg, quietly = TRUE)) {
      packageStartupMessage(
        sprintf("DISQ: package '%s' is not installed. Install it to enable related features.", pkg)
      )
      return(invisible(FALSE))
    }

    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    )
    invisible(TRUE)
  }

  attach_if_installed("BPST")
  attach_if_installed("Triangulation")

  invisible(NULL)
}
