#' Run the MCMChybridGP demonstration
#'
#' This is the user-facing demo function. It wraps the internal
#' \code{.runMCMCdemo()} helper and provides a clean, exported entry point.
#'
#' @param ... Passed to the internal demo function.
#'
#' @return Invisibly returns the demo output.
#' @export
Demo <- function(...) {
    .runMCMCdemo(...)
}

