# =============================================================================
# Sample Splitting Utilities
# Functions for creating train/test splits and cross-validation folds
# =============================================================================

#' Create simple train/test split
#'
#' @param n Sample size
#' @param train.proportion Proportion of data for training (default 0.5)
#' @param seed Random seed (optional)
#'
#' @return List with components:
#'   \item{train}{Indices for training set}
#'   \item{test}{Indices for test set}
#'
#' @keywords internal
make_split <- function(n, train.proportion = 0.5, seed = NULL, clusters = NULL) {

  validate_split_params(train.proportion)

  if (!is.null(seed)) set.seed(seed)

  if (!is.null(clusters)) {
    unique_cl  <- unique(clusters)
    n_cl       <- length(unique_cl)
    train_size <- floor(train.proportion * n_cl)
    train_cl   <- sample(unique_cl, size = train_size, replace = FALSE)
    train_idx  <- which(clusters %in% train_cl)
    test_idx   <- which(!clusters %in% train_cl)
  } else {
    train_size <- floor(train.proportion * n)
    train_idx  <- sample(seq_len(n), size = train_size, replace = FALSE)
    test_idx   <- setdiff(seq_len(n), train_idx)
  }

  list(train = train_idx, test = test_idx)
}


#' Create K-fold cross-validation assignments
#'
#' @param n Sample size
#' @param nfolds Number of folds (default 5)
#' @param seed Random seed (optional)
#'
#' @return Integer vector of fold assignments (1, 2, ..., nfolds)
#'
#' @keywords internal
make_folds <- function(n, nfolds = 5, seed = NULL, clusters = NULL) {

  validate_cv_params(nfolds, n)

  if (!is.null(seed)) set.seed(seed)

  if (!is.null(clusters)) {
    unique_cl <- unique(clusters)
    n_cl      <- length(unique_cl)
    if (nfolds > n_cl)
      stop(sprintf("nfolds (%d) cannot exceed number of clusters (%d)", nfolds, n_cl))
    cl_folds          <- sample(rep(seq_len(nfolds), length.out = n_cl))
    names(cl_folds)   <- as.character(unique_cl)
    return(cl_folds[as.character(clusters)])
  } else {
    return(sample(rep(seq_len(nfolds), length.out = n)))
  }
}

