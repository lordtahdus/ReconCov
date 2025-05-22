#' @export
plot_heatmap <- function(mat) {
  ggplot(reshape2::melt(mat), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1, 1), name = "Corr") +
    theme_minimal() +
    coord_fixed()
}
