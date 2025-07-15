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


#' Transform an length-M list of models' error2 matrices into df format
#'
#' @param error_list A length-M list of list of error matrices where each matrix
#'                   corresponds to a model's error squared values.
#' @return A data frame containing the error squared values, model names,
#'         series names, and forecast horizons, with id represent simulation rep.
#'
#' @examples
#' error_list <- list(
#'   base = matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 3),
#'   mint_shr = matrix(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7), nrow = 3)
#' )
#' error_df <- transform_error_list(error_list)
#'
#' @importFrom purrr imap_dfr map_dfr
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#'
#' @export
transform_error_list <- function(error_list) {
  imap_dfr(
    error_list, ~ {
      sim_id   <- .y                  # 1, 2, ..., M
      models   <- .x                  # list(base = matrix, mint_shr = matrix, …)
      map_dfr(names(models), function(model) {
        mat <- models[[model]]
        as.data.frame(mat) %>%        # convert matrix to data.frame, columns = series
          mutate(h = seq_len(nrow(mat))) %>%
          pivot_longer(-h,
                       names_to  = "series",
                       values_to = "e2") %>%
          mutate(.model = model,
                 id     = sim_id)
      })
    }
  )
}
