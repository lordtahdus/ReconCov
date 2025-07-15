#' Plot a heatmap of a matrix
#'
#' @param mat A matrix to be plotted as a heatmap.
#' @param is_unit_spherical Logical, if TRUE, the color scale is set to [-1, 1].
#' @return A ggplot object representing the heatmap.
#'
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10)
#' plot_heatmap(mat)
#' plot_heatmap(mat, is_unit_spherical = TRUE)
#'
#' @import ggplot2
#' @export
plot_heatmap <- function(mat, is_unit_spherical = FALSE) {
  if (is_unit_spherical) {
    limit <- c(-1, 1)
    midpoint <- 0
  }
  else {
    max_val <- max(abs(mat), na.rm = TRUE)
    limit <- c(-max_val, max_val)
    midpoint <- 0
  }
  if (!is.null(colnames(mat)) || !is.null(rownames(mat))) {
    colnames(mat) <- NULL
    rownames(mat) <- NULL
  }
  ggplot(reshape2::melt(mat), aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = midpoint, limit = limit, name = "Corr") +
    scale_y_reverse() +
    theme_minimal() +
    coord_fixed()
}
