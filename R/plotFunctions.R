MADplots <- function(ggData, xvar = "speciesCode", yvar = "MAD", colourvar = "variable",
                     xlabs = NULL, collabs = NULL) {
  gg <- ggplot(ggData, aes_string(x = xvar, y = yvar, colour = colourvar)) +
    stat_summary(fun.data = "mean_sdl", position = position_dodge(width = 0.5))

  if (!is.null(xlabs)) {
   gg <- gg + scale_x_discrete(labels = xlabs)
  }

  if (!is.null(collabs)) {
    gg <- gg + scale_color_brewer(palette = "Dark2", labels = collabs)
  } else {
    gg <- gg + scale_color_brewer(palette = "Dark2")
  }
  gg <- gg + plotTheme(base_size = 12, margin = FALSE) +
    labs(colour = "", x = "") +
    theme(strip.text.x = element_blank(), strip.background = element_blank()) +
    facet_wrap(colourvar, ncol = 1, scales = "free_y")

  return(gg)
}
