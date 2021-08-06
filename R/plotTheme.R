#' Theme for validation plots
#'
#' Based on ggpubr::theme_pubr but with added major y axis grid lines
#'
#' @param ... passed to ggpubr::theme_pubr

plotTheme <- function(majorYlines = TRUE, ...) {
  .theme <- theme_pubr(...)
  if (majorYlines) {
    .theme <- .theme +
      theme(panel.grid.major.y = element_line(colour = "grey", linetype = "dotted", size = 0.5))
  }
  .theme
}
