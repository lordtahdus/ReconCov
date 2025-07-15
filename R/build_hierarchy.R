#' Construct a summing matrix for hierarchical time series
#'
#' Accepts a single `structure` describing all levels: the bottom block sizes and
#' how to merge them successively into upper levels.
#'
#' @param structure  List whose first element is an integer vector of bottom block sizes,
#'                   and each subsequent element is a list of integer vectors indicating
#'                   which blocks to merge at that level.
#' @param ascending Logical; if TRUE (default), the S matrix is ordered from bottom to top.
#' @param sparse     Logical; if TRUE (default), returns a sparse dgCMatrix (from Matrix).
#'
#' @return  The S matrix (rows = all aggregation levels, top down, then bottom identity).
#' @examples
#' S <- construct_S(
#'   structure = list(
#'     # bottom: 4 groups of size 3
#'     c(3,3,3,3),
#'     # next: merge blocks 1+2, merge 3+4
#'     list(c(1,2), c(3,4)),
#'     # top: merge the two groups
#'     list(c(1,2))
#'   )
#' )
#'
#' @export
construct_S <- function(
    structure,
    ascending = TRUE,
    sparse = TRUE
) {
  if(!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install the 'Matrix' package for sparse support.")
  }

  # Helper: split a pipe-delimited string
  split_pipe <- function(x) strsplit(x, "\\|")[[1]]

  # Generate bottom labels A1, A2, ... for each block of given size
  bottom_sizes <- structure[[1]]
  bottom_blocks <- mapply(
    function(size, idx) paste0("A", idx, seq_len(size)),
    bottom_sizes, seq_along(bottom_sizes), SIMPLIFY = FALSE
  )
  bottom_labels <- unlist(bottom_blocks)

  # Build one sparse block-per-level, top-down
  levels <- list()
  prev_blocks <- bottom_blocks

  for(level_groups in structure[-1]) {
    parent_ids <- seq_along(level_groups)
    parent_labels <- paste0(LETTERS[length(levels)+2], parent_ids)
    # map each parent to the children labels
    mapping <- setNames(
      lapply(level_groups, function(g) unlist(prev_blocks[g])),
      parent_labels
    )
    # make one sparse submatrix
    rows <- rep(seq_along(mapping), vapply(mapping, length, 0L))
    cols <- match(unlist(mapping), unlist(prev_blocks))
    blk <- sparseMatrix(
      i = rows, j = cols, x = 1L,
      dims = c(length(mapping), length(unlist(prev_blocks))),
      dimnames = list(names(mapping), unlist(prev_blocks))
    )
    levels[[length(levels)+1]] <- blk
    prev_blocks <- mapping
  }

  # Identity for bottom level
  Ibot <- Diagonal(length(bottom_labels))
  dimnames(Ibot) <- list(bottom_labels, bottom_labels)

  # Combine: all upper blocks first, then the bottom identity
  if (ascending) {
    S <- do.call(rbind, c(list(Ibot), levels))
  } else {
    S <- do.call(rbind, c(levels[length(levels) : 1], list(Ibot)))
  }

  if(!sparse) S <- as.matrix(S)
  S
}


############### ARCHIVED #################


# build_2level_hierarchy <- function(bottom_data, groups) {
#
#   T <- nrow(bottom_data)
#   p <- ncol(bottom_data)
#
#   # Number of groups
#   G <- length(groups)
#
#   # Level 1 (Group-level) aggregation:
#   group_indices <- split(1:p, rep(1:G, groups))
#   S_group <- matrix(0, nrow = G, ncol = p)
#   for (i in seq_len(G)) {
#     S_group[i, group_indices[[i]]] <- 1
#   }
#
#   # Level 2 (Top-level) aggregation:
#   S_top <- matrix(1, nrow = 1, ncol = p)
#   # Bottom-level is the identity matrix.
#   S_bottom <- diag(1, nrow = p, ncol = p)
#
#   S <- rbind(S_bottom, S_group, S_top)
#   hier_data <- bottom_data %*% t(S)
#
#   # Create series names.
#   series_names <- c(paste0("B", 1:p), paste0("G", seq_len(G)), "Top")
#   colnames(hier_data) <- series_names
#   rownames(S) <- series_names
#
#   return(list(
#     hier_data = hier_data,
#     S = S,
#     series_names = series_names
#   ))
# }



#' #' Build a hierarchy from bottom level series
#' #'
#' #' @param bottom_data A numeric matrix of dimension T x p containing the bottom-level series
#' #' @param groups Integer vector specifying how many bottom series belong to each group
#' #'   (e.g. c(2,3,4,3)).
#' #'
#' #' @export
#' build_hierarchy <- function(bottom_data, groups) {
#'   if (!is.matrix(bottom_data)) {
#'     stop("bottom_data must be a matrix of dimension T x p.")
#'   }
#'   T <- nrow(bottom_data)
#'   p <- ncol(bottom_data)
#'
#'   if (sum(groups) != p) {
#'     stop("sum(groups) must match ncol(bottom_data).")
#'   }
#'
#'   # 1) Bottom level: we keep these as identity
#'   #    S so far is the identity (p x p).
#'   S_current <- diag(1, nrow = p, ncol = p)
#'   # hier_data so far is just bottom_data
#'   hier_data <- bottom_data
#'   total_series_count <- p  # so far
#'
#'   # We'll store some metadata about how many new series we add at each step
#'   level_info <- list()
#'   level_info[["Level_bottom"]] <- p
#'
#'   # 2) Aggregate each group => one series per group.
#'   #    We'll create k = length(groups) new aggregator rows in S.
#'   group_starts <- c(1, cumsum(groups) + 1)
#'   group_starts <- group_starts[-length(group_starts)]  # remove last
#'   k <- length(groups)
#'
#'   group_agg_indices <- list()
#'   for(i in seq_len(k)) {
#'     idx_start <- group_starts[i]
#'     idx_end   <- idx_start + groups[i] - 1
#'     group_agg_indices[[i]] <- idx_start:idx_end
#'   }
#'
#'   # We'll build new aggregator rows for each group.
#'   group_rows <- matrix(0, nrow = k, ncol = p)
#'
#'   for(i in seq_len(k)) {
#'     group_rows[i, group_agg_indices[[i]]] <- 1
#'   }
#'
#'   # Bind these aggregator rows to the existing S_current
#'   S_new1 <- rbind(S_current, group_rows)
#'   # Summation => we get T x k aggregator columns
#'   group_data <- bottom_data %*% t(group_rows)
#'
#'   # Update hier_data
#'   hier_data <- cbind(hier_data, group_data)
#'   total_series_count <- total_series_count + k
#'   level_info[["Level_groups"]] <- k
#'
#'   # 3) Suppose we want to aggregate those "k group-level series" in pairs
#'   #    For example, if k=4, we do pairs: (1,2), (3,4).
#'   #    Adjust code to handle odd number of groups if needed.
#'   #    We'll produce new aggregator rows in S for these pairs.
#'
#'   if (k %% 2 != 0) {
#'     stop("For illustration, we assume k is even, so we can pair them up.")
#'   }
#'   nPairs <- k / 2
#'
#'   pair_rows <- matrix(0, nrow = nPairs, ncol = p)
#'   # which portion of S_new1 do these aggregator rows correspond to? They refer to the newly made group rows.
#'
#'   # The aggregator rows for groups are appended at the bottom of S_new1 => row indices (p+1):(p+k)
#'   # We want to sum consecutive pairs of those aggregator rows, e.g. row p+1 + row p+2 => new aggregator, ...
#'   # We'll do: pairIndices for the new group-level series:
#'   group_row_indices <- (p+1) : (p+k)
#'   # We form pairs: c(1,2), c(3,4), ...
#'   # for each pair, we sum those aggregator rows
#'   for(pair_i in seq_len(nPairs)) {
#'     # This is 2 consecutive group-level aggregator rows
#'     r1 <- group_row_indices[2*pair_i - 1]
#'     r2 <- group_row_indices[2*pair_i]
#'     # Summation in S => row r1 + row r2
#'     pair_rows[pair_i, ] <- S_new1[r1, ] + S_new1[r2, ]
#'   }
#'
#'   S_new2 <- rbind(S_new1, pair_rows)
#'   # The actual aggregator data for these pairs is just group_data[ ,2*pair_i-1] + group_data[,2*pair_i].
#'   # We'll do that in code:
#'   pair_data <- matrix(0, nrow = T, ncol = nPairs)
#'   for(pair_i in seq_len(nPairs)) {
#'     pair_data[, pair_i] <- group_data[, 2*pair_i - 1] + group_data[, 2*pair_i]
#'   }
#'
#'   # Update hier_data
#'   hier_data <- cbind(hier_data, pair_data)
#'   total_series_count <- total_series_count + nPairs
#'   level_info[["Level_pairs"]] <- nPairs
#'
#'   # 4) Finally, aggregate the nPairs super-groups => 1 top series
#'   if (nPairs %% 2 != 0) {
#'     stop("For illustration, we assume the next step is to combine the 2 pairs => top series.")
#'   }
#'   # We'll do 1 aggregator row that sums all nPairs or do we do pairs again?
#'   # Let's follow the user example: if nPairs=2 => we sum them => top series
#'   top_row <- matrix(0, nrow = 1, ncol = p)
#'   # The row indices for the pairs are (p+k+1) : (p+k+nPairs) in S_new2
#'   # Sum them all
#'   pair_row_indices <- (p+k+1) : (p+k+nPairs)
#'   top_row[1, ] <- colSums(S_new2[pair_row_indices, , drop=FALSE])
#'
#'   S_final <- rbind(S_new2, top_row)
#'
#'   top_data <- rowSums(pair_data)
#'   # add to hier_data
#'   hier_data <- cbind(hier_data, top_data)
#'   total_series_count <- total_series_count + 1
#'   level_info[["Level_top"]] <- 1
#'
#'   # Summation matrix is S_final => dimension is total_series_count x p
#'   # hier_data is T x total_series_count
#'
#'   # Return
#'   out <- list(
#'     hier_data = hier_data,
#'     S = S_final,
#'     info = level_info
#'   )
#'   return(out)
#' }
