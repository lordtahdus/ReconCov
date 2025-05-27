#' Transform a list of MSE matrices into a tsibble/tibble format
#'
#' @param MSE A list of matrices where each matrix corresponds to a model's MSE values.
#' @param to_ts Logical indicating whether to convert the result to a tsibble.
#' @return A tsibble or tibble containing the MSE values, model names, series names, and forecast horizons.
#'
#' @examples
#' MSE <- list(
#'   model1 = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 3),
#'   model2 = matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), nrow = 3)
#' )
#' result <- transform_sim_MSE(MSE, to_ts = TRUE)
#'
#' @importFrom purrr imap_dfr
#' @importFrom dplyr mutate select row_number
#' @importFrom tidyr pivot_longer
#' @importFrom tsibble as_tsibble
#' @importFrom tibble as_tibble
#' @export
transform_sim_MSE <- function(MSE, to_ts = TRUE) {
  df <- imap_dfr(MSE, function(mat, model_name) {
    as_tibble(mat) %>%
      mutate(h = row_number()) %>%
      pivot_longer(-h, names_to = "series", values_to = "MSE") %>%
      mutate(.model = model_name) %>%
      select(.model, series, h, MSE)
  })
  if (to_ts) {
    return(as_tsibble(df, key = c(.model, series), index = h))
  }
  return(df)
}
