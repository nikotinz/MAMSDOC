#' .onAttach
#'
#' Print MAMS package info
#' @param ... Further arguments passed to or from other methods.
.onAttach <- function(...) {
  packageStartupMessage(paste("********** MAMS Version",
                      packageDescription("MAMS")$Version), " ********** \n")
  packageStartupMessage(paste0(
    "Type MAMSNews() to see new features/changes/bug fixes.\n",
    "Type help(MAMS) for an overview of the package.\n"))
}
