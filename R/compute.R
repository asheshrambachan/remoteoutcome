# =============================================================================
# Core RSV Computation Functions
# Mathematical functions for computing RSV treatment effect estimates
# =============================================================================

#' ATE-specific relevance: (lambda' J-inv lambda)-inv, optionally pi-weighted over strata
#'
#' For a single stratum (no X): \eqn{(\lambda' J^{-1} \lambda)^{-1}}.
#' For K=2 with outcomes in \eqn{\{0,1\}} and yK=0 this reduces to J (the existing scalar).
#' For multiple strata: pi-weighted average of per-stratum relevances,
#'   \eqn{\sum_x \pi(x) (\lambda' J(x)^{-1} \lambda)^{-1}}
#' where pi_weights must be supplied (same weights used for theta_ATE).
#' Returns NA when J is non-invertible or all strata are NA.
#'
#' @keywords internal
compute_relevance_ate <- function(denom_raw, y_levels, yK, pi_weights = NULL) {
  yj_levels <- setdiff(y_levels, yK)
  lambda    <- as.numeric(yj_levels) - as.numeric(yK)
  scalar_rel <- function(J) {
    if (!is.matrix(J) || !all(is.finite(J))) return(NA_real_)
    tryCatch(
      as.numeric(1 / drop(t(lambda) %*% solve(J) %*% lambda)),
      error = function(e) NA_real_
    )
  }
  if (length(denom_raw) == 1) {
    scalar_rel(denom_raw[[1]])
  } else if (!is.null(pi_weights) && length(pi_weights) == length(denom_raw)) {
    rels  <- vapply(denom_raw, scalar_rel, FUN.VALUE = numeric(1))
    valid <- !is.na(rels)
    if (!any(valid)) return(NA_real_)
    sum(pi_weights[valid] * rels[valid]) / sum(pi_weights[valid])
  } else {
    NA_real_
  }
}

#' Ensure S_e and S_o are numeric
#'
#' Converts logical sample indicators to numeric (0/1) if needed.
#'
#' @param data Data frame with S_e and S_o columns
#'
#' @return Data frame with numeric S_e and S_o
#'
#' @keywords internal
.ensure_numeric_indicators <- function(data) {
  if (is.logical(data$S_e)) {
    data$S_e <- as.numeric(data$S_e)
  }
  if (is.logical(data$S_o)) {
    data$S_o <- as.numeric(data$S_o)
  }
  data
}

.get_yK_value <- function(y_levels, theta_init){
  yj_cols <- names(theta_init) |>
    grep(pattern = "^theta_Y", value = TRUE) |>
    sub(pattern = "^theta_Y", replacement = "") |>
    type.convert(as.is = TRUE)
  yK <- setdiff(y_levels, yj_cols)
  yK
}

#' Compute stratified counts for RSV estimation
#'
#' Computes counts of observations within experimental and observational samples,
#' stratified by discrete covariates X (if provided). Returns normalized counts
#' used in delta computations.
#'
#' @param df Data frame with Y, D, S_e, S_o, and optionally X columns
#' @param X_cols Character vector of column names for discrete covariates (optional)
#'
#' @return Data frame with normalized counts:
#' \describe{
#'   \item{count_D1_Se}{Proportion with D=1 in experimental sample, by X stratum}
#'   \item{count_D0_Se}{Proportion with D=0 in experimental sample, by X stratum}
#'   \item{count_Y*_So}{Proportion with Y=* in observational sample, by X stratum}
#' }
#'
#' @keywords internal
compute_counts <- function(df, y_levels, X_cols = character(0)) {
  # Base denominator: COUNT(X=x) or total n if no X
  base <- df |>
    group_by(across(all_of(X_cols))) |>
    summarise(denom = dplyr::n(), .groups = "drop")
  
  # DSe counts within S_e==1
  DSe <- df |>
    dplyr::filter(.data$S_e == 1) |>
    group_by(across(all_of(X_cols))) |>
    summarise(
      D1 = sum(.data$D == 1, na.rm = TRUE),
      D0 = sum(.data$D == 0, na.rm = TRUE),
      .groups = "drop"
    )

  # YSo counts within S_o==1 (wide format)
  YSo <- df |>
    dplyr::filter(.data$S_o == 1) |>
    mutate(
      Y = factor(.data$Y, levels = y_levels)
      ) |>
    count(across(all_of(c(X_cols, "Y"))), name = "cnt", .drop = FALSE) |>
    tidyr::pivot_wider(
      names_from   = "Y",
      values_from  = "cnt",
      values_fill  = 0,
      names_prefix = "Y"
    )
  # Set zero counts to NA, but only on the count columns — never on X_cols which may equal 0
  cnt_cols <- setdiff(names(YSo), X_cols)
  YSo[cnt_cols][YSo[cnt_cols] == 0] <- NA

  # Join; if no X_cols, cross join of one-row tables
  if (length(X_cols) == 0) {
    out <- base |>
      cross_join(DSe) |>
      cross_join(YSo)
  } else {
    out <- base |>
      left_join(DSe, by = X_cols) |>
      left_join(YSo, by = X_cols)
  }

  out <- out |>
    mutate(
      # Normalize. If S_e==1 subset had no rows for some group, D1/D0 will be NA -> treat as 0
      across(starts_with("D", ignore.case = FALSE), ~ .x / .data$denom, .names = "count_{.col}_Se"),

      # Normalize all Y* columns (which are counts) into COUNT(Y=y,S_o|X=x)
      across(starts_with("Y"), ~ .x / .data$denom, .names = "count_{.col}_So")
    ) |>
    select(-all_of("denom"), -starts_with("D", ignore.case = FALSE), -starts_with("Y", ignore.case = FALSE))

  return(out)
}


#' Compute delta functions for RSV estimation
#'
#' Computes delta_e (experimental delta) and delta_o (observational deltas)
#' used in the RSV estimating equations. Can use either predicted or observed
#' values.
#'
#' @param df Data frame with observations and predictions
#' @param yK Reference category for Y
#' @param X_cols Character vector of column names for discrete covariates (optional)
#' @param use_pred Logical; if TRUE, use predictions (pred_*), if FALSE use observations
#'
#' @return Data frame with delta columns:
#' \describe{
#'   \item{delta_e}{Experimental delta (single column)}
#'   \item{delta_o_Y*}{Observational deltas (one per non-reference Y category)}
#'   \item{X_cols}{X columns if provided}
#' }
#'
#' @details
#' For use_pred = TRUE:
#' \itemize{
#'   \item delta_e = (pred_D / count_D1 - (1-pred_D) / count_D0) * pred_S_e
#'   \item delta_o_j = (pred_Yj / count_Yj - pred_YK / count_YK) * pred_S_o
#' }
#'
#' For use_pred = FALSE: Uses observed Y, D instead of predictions.
#'
#' @keywords internal
compute_deltas <- function(df, y_levels, yK, X_cols = character(0), use_pred = TRUE) {

  count_df <- compute_counts(df, y_levels = y_levels, X_cols = X_cols)
  
  if (length(X_cols) == 0) {
    merged <- bind_cols(df, count_df)
  } else {
    merged <- left_join(df, count_df, by = X_cols)
  }
  
  # Δ^e(x)
  if (use_pred) {
    merged <- merged |>
      mutate(delta_e = (.data$pred_D / .data$count_D1_Se - (1 - .data$pred_D) / .data$count_D0_Se) * .data$pred_S_e)
  } else {
    merged <- merged |>
      mutate(delta_e = (.data$D == 1 & .data$S_e == 1) / .data$count_D1_Se - (.data$D == 0 & .data$S_e == 1) / .data$count_D0_Se)
  }

  # Split levels: j = all except K
  yj_levels <- setdiff(y_levels, yK)

  # Δ^o_j(x)
  if (use_pred) {
    pred_yj_cols <- paste0("pred_Y", yj_levels)
    pred_yK_col <- paste0("pred_Y", yK)
    count_yj_cols <- paste0("count_Y", yj_levels, "_So")
    count_yK_col  <- paste0("count_Y", yK, "_So")

    pred_yj <- as.matrix(merged[, pred_yj_cols, drop = FALSE])
    count_yj <- as.matrix(merged[, count_yj_cols, drop = FALSE])
    pred_yK <- merged[[pred_yK_col]]
    count_yK <- merged[[count_yK_col]]
    delta_o <- (pred_yj / count_yj - pred_yK / count_yK) * merged$pred_S_o
    delta_o <- delta_o  # Ensure matrix even with 1 column
    colnames(delta_o) <- gsub("^pred_", "delta_o_", pred_yj_cols)
  } else {
    merged$Y <- factor(merged$Y, levels = y_levels)
    Y_df <- model.matrix(~ Y - 1, model.frame(~ Y, merged, na.action = na.pass))
    Y_df[is.na(Y_df)] <- 0
    merged <- cbind(merged, as.data.frame(Y_df))

    true_yj_cols <- paste0("Y", yj_levels)
    true_yK_col <- paste0("Y", yK)
    count_yj_cols <- paste0("count_Y", yj_levels, "_So")
    count_yK_col  <- paste0("count_Y", yK, "_So")

    true_yj <- as.matrix(merged[, true_yj_cols, drop = FALSE])
    true_yK <- merged[[true_yK_col]]
    count_yj <- as.matrix(merged[, count_yj_cols, drop = FALSE])
    count_yK   <- merged[[count_yK_col]]
    delta_o <- (true_yj == 1 & merged$S_o == 1) / count_yj -  (true_yK == 1 & merged$S_o == 1) / count_yK
    delta_o <- delta_o # Ensure matrix even with 1 column
    colnames(delta_o) <- paste0("delta_o_", true_yj_cols)
  }

  cbind(merged, delta_o) |> 
    select(all_of(X_cols), starts_with("delta_"))
}


#' Solve least squares for theta estimation
#'
#' Internal helper function to solve least squares regression without intercept.
#' Used in computing initial theta estimates.
#'
#' @param y Response vector
#' @param X Design matrix (without intercept)
#'
#' @return One-row data frame with theta estimates as named columns
#'
#' @keywords internal
solve_ls <- function(y, X) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  
  # Drop rows where all X entries are NA
  ok_rows <- rowSums(!is.na(X)) > 0
  y <- y[ok_rows]
  X <- X[ok_rows, , drop = FALSE]
  
  # If nothing left
  if (length(y) == 0 || nrow(X) == 0) {
    theta <- rep(NA_real_, ncol(X))
    names(theta) <- colnames(X)
    return(as.data.frame(as.list(theta), check.names = FALSE))
  }
  
  # Drop columns that are all NA
  ok_cols <- colSums(!is.na(X)) > 0
  
  # Fit on remaining columns
  if (any(ok_cols)) {
    fit <- lm.fit(x = X[, ok_cols, drop = FALSE], y = y)
    theta <- rep(NA_real_, ncol(X))
    theta[ok_cols] <- fit$coefficients
  } else {
    theta <- rep(NA_real_, ncol(X))
  }
  
  names(theta) <- colnames(X)
  as.data.frame(as.list(theta), check.names = FALSE)
}


#' Compute initial theta estimate
#'
#' Computes the initial treatment effect estimate theta using observed sample
#' indicators (S_e, S_o) as predictions for the sample selection probabilities.
#' This provides a starting value for the iterative RSV procedure.
#'
#' @param df Data frame with observations and predictions (pred_Y*, pred_D)
#' @param yK Reference category for Y
#' @param X_cols Character vector of column names for discrete covariates (optional)
#'
#' @return Data frame with initial theta estimates:
#' \describe{
#'   \item{theta_Y*}{Initial theta for each non-reference Y category}
#'   \item{X_cols}{X columns if provided}
#' }
#'
#' @keywords internal
compute_theta_init <- function(df, y_levels, yK = NULL, X_cols = character(0)) {
  if (is.null(yK)) {
    yK <- y_levels[1]
  }
  delta_df <- compute_deltas(df, y_levels = y_levels, yK = yK, X_cols = X_cols, use_pred = TRUE)
  
  y_col <- "delta_e"
  B_cols <- delta_df |>
    select(starts_with("delta_o_Y")) |>
    colnames()
  theta_init <- delta_df |>
    group_by(across(all_of(X_cols))) |>
    group_modify(~ solve_ls(.x[[y_col]], .x[, B_cols, drop = FALSE])) |>
    ungroup()
  
  # Rename theta columns to "theta_*" for easy merging later
  theta_init |>
    rename_with(~ gsub("^delta_o_", "theta_", .x), .cols = all_of(B_cols))
}


#' Compute sigma-squared (variance) for efficient weighting
#'
#' Computes the variance sigma-squared used in constructing efficient weights
#' for the RSV estimator. Applies lower bound stabilization.
#'
#' @param df Data frame with observations and predictions
#' @param theta_init Data frame with theta estimates (by X stratum if applicable)
#' @param sigma2.lower Lower bound for sigma-squared stabilization (optional)
#' @param sigma2.quantile Quantile for sigma-squared lower bound if sigma2.lower not provided
#'
#' @return Vector of sigma-squared values (length n)
#'
#' @details
#' Computes sigma^2 = term1 + term2 + term3 where:
#' \itemize{
#'   \item term1: Experimental sample variance component
#'   \item term2: Observational sample variance component (sum of theta^2)
#'   \item term3: Observational sample variance component (cross term)
#' }
#'
#' Stabilization: sigma^2 is bounded below by its sigma2.quantile-th quantile
#' to prevent numerical instability.
#'
#' @keywords internal
compute_sigma2 <- function(
    df,
    theta_init,
    y_levels,
    sigma2.lower = NULL,
    sigma2.quantile = 0.01
  ) {
  X_cols <- names(theta_init)[!startsWith(names(theta_init), "theta_Y")]
  
  count_df <- compute_counts(df, y_levels = y_levels, X_cols = X_cols)
  
  if (length(X_cols) == 0) {
    df <- df |>
      bind_cols(count_df) |>
      bind_cols(theta_init)
  } else {
    df <- df |>
      left_join(count_df, by = X_cols) |>
      left_join(theta_init, by = X_cols)
  }

  # Split levels: j = all except K
  yK <- .get_yK_value(y_levels, theta_init)
  yj_levels <- setdiff(y_levels, yK)
  
  # Collect needed column names
  theta_cols <- paste0("theta_Y", yj_levels)
  pred_yj_cols <- paste0("pred_Y", yj_levels)
  count_yj_cols <- paste0("count_Y", yj_levels, "_So")
  count_yK_col  <- paste0("count_Y", yK, "_So")
  pred_yK_col <- paste0("pred_Y", yK)

  stopifnot(all(theta_cols %in% names(df)))
  stopifnot(all(pred_yj_cols %in% names(df)))
  stopifnot(all(count_yj_cols %in% names(df)))
  stopifnot(count_yK_col %in% names(df))
  stopifnot(pred_yK_col %in% names(df))

  # Matrices for vectorized rowwise sums
  theta_init <- as.matrix(df[, theta_cols, drop = FALSE])  # n x (K-1)
  pred_yj <- as.matrix(df[, pred_yj_cols, drop = FALSE])
  count_yj <- as.matrix(df[, count_yj_cols, drop = FALSE])
  pred_yK <- df[[pred_yK_col]]
  count_yK <- df[[count_yK_col]]

  # Compute variance components
  term1 <- (df$pred_D / (df$count_D1_Se^2) + (1 - df$pred_D) / (df$count_D0_Se^2)) * df$pred_S_e
  term2 <- rowSums(theta_init^2 * pred_yj / count_yj^2, na.rm = TRUE) * df$pred_S_o
  term3 <- rowSums(theta_init, na.rm = TRUE)^2 * pred_yK / count_yK^2 * df$pred_S_o
  sigma2 <- term1 + term2 + term3
  
  # Stabilization: lower bound sigma2
  if (is.null(sigma2.lower)) {
    sigma2.lower <- quantile(sigma2, sigma2.quantile, na.rm = TRUE)
  }
  sigma2 <- pmax(sigma2, sigma2.lower)

  return(sigma2)
}


#' Compute H matrix (efficient weights)
#'
#' Computes the H matrix used in efficient weighting for the RSV estimator.
#' H = delta_o / sigma^2.
#'
#' @param df Data frame with observations and predictions
#' @param theta_init Data frame with theta estimates
#' @param sigma2.lower Lower bound for sigma-squared stabilization (optional)
#' @param sigma2.quantile Quantile for sigma-squared lower bound
#'
#' @return Data frame with H columns (one per non-reference Y category)
#'
#' @keywords internal
compute_H <- function(
    df,
    theta_init,
    y_levels,
    sigma2.lower = NULL,
    sigma2.quantile = 0.01,
    pred_y_only = FALSE
    ) {

  yK <- .get_yK_value(y_levels, theta_init)

  if (pred_y_only) {
    yj_levels    <- setdiff(y_levels, yK)
    pred_yj_cols <- paste0("pred_Y", yj_levels)
    pred_yj      <- as.matrix(df[, pred_yj_cols, drop = FALSE])
    out          <- as.data.frame(pred_yj)
    names(out)   <- paste0("H_Y", yj_levels)
    return(out)
  }

  X_cols <- names(theta_init)[!startsWith(names(theta_init), "theta_Y")]
  delta_o <- compute_deltas(df, y_levels = y_levels, yK = yK, X_cols = X_cols, use_pred = TRUE) |>
    select(starts_with("delta_o_Y"))

  sigma2 <- compute_sigma2(df, theta_init, y_levels = y_levels, sigma2.lower = sigma2.lower,
                           sigma2.quantile = sigma2.quantile)

  out <- delta_o / sigma2
  names(out) <- gsub("^delta_o_", "H_", names(out))
  out
}


#' Compute final theta estimate
#'
#' Computes the final efficient theta estimate by solving the estimating
#' equations with efficient weights. Stratifies by X if provided.
#'
#' @param df Data frame with observations and predictions
#' @param theta_init Data frame with initial theta estimates
#' @param sigma2.lower Lower bound for sigma-squared stabilization (optional)
#' @param sigma2.quantile Quantile for sigma-squared lower bound
#'
#' @return List with:
#' \describe{
#'   \item{theta}{Final theta estimates (by X stratum if applicable)}
#'   \item{numerator}{Numerator of estimating equation}
#'   \item{denominator}{Denominator of estimating equation}
#'   \item{H}{Efficient weights}
#' }
#'
#' @details
#' Solves the weighted estimating equation:
#' sum_i H_i * (delta_e_i - delta_o_i' * theta) = 0
#'
#' which gives theta = \eqn{(\sum_i H_i \delta_{oi} \delta_{oi}')^{-1} (\sum_i H_i \delta_{oi} \delta_{ei})}
#'
#' @keywords internal
compute_theta <- function(
    df,
    theta_init,
    y_levels,
    sigma2.lower = NULL,
    sigma2.quantile = 0.01,
    pred_y_only = FALSE,
    weights = NULL
  ) {

  # 1) Compute H (same rows as df)
  H_df <- compute_H(df, theta_init, y_levels = y_levels, sigma2.lower = sigma2.lower,
                    sigma2.quantile = sigma2.quantile, pred_y_only = pred_y_only)
  
  # 2) Compute deltas (same rows as df)
  yK <- .get_yK_value(y_levels, theta_init)
  X_cols <- names(theta_init)[!startsWith(names(theta_init), "theta_Y")]
  delta_df <- compute_deltas(df, y_levels = y_levels, yK = yK, X_cols = X_cols, use_pred = FALSE)
  merged_df <- dplyr::bind_cols(delta_df, H_df)
  if (!is.null(weights)) merged_df[[".__bb_w__"]] <- weights
  
  delta_o_cols <- names(merged_df)[grepl("^delta_o_Y", names(merged_df))]
  H_cols       <- names(merged_df)[grepl("^H_Y", names(merged_df))]

  stopifnot(length(delta_o_cols) == length(H_cols))
  stopifnot("delta_e" %in% names(merged_df))

  solve_group <- function(d) {
    B <- as.matrix(d[, H_cols, drop = FALSE])        # n x (K-1)
    O <- as.matrix(d[, delta_o_cols, drop = FALSE])  # n x (K-1)
    e <- as.matrix(d[, "delta_e", drop = FALSE])     # n
    n_d <- nrow(d)
    if (".__bb_w__" %in% names(d)) {
      w     <- d[[".__bb_w__"]]
      denom <- crossprod(B * w, O)
      num   <- crossprod(B * w, e)
    } else {
      denom <- crossprod(B, O) / n_d  # (K-1) x (K-1)
      num   <- crossprod(B, e) / n_d  # (K-1) x 1
    }
    
    theta <- as.numeric(solve_ls(num, denom))
    names(theta) <- gsub("^delta_o_", "theta_", delta_o_cols)
    
    list(
      theta = tibble::as_tibble(as.list(theta)),
      numerator = tibble::tibble(numerator = list(num)),
      denominator = tibble::tibble(denominator = list(denom))
    )
  }
  
  perX <- merged_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(X_cols))) |>
    dplyr::group_modify(~{
      res <- solve_group(.x)
      dplyr::bind_cols(res$theta, res$numerator, res$denominator)
    }) |>
    dplyr::ungroup()
  
  list(
    theta = perX |> dplyr::select(dplyr::all_of(X_cols), dplyr::starts_with("theta_")),
    numerator = perX |> dplyr::select(dplyr::all_of(X_cols), "numerator") ,
    denominator = perX |>  dplyr::select(dplyr::all_of(X_cols), "denominator"),
    H = H_df
  )
}


#' Compute average treatment effect (ATE) from stratified theta
#'
#' Averages theta estimates across X strata with weights proportional to
#' stratum sizes. If no X, returns theta directly.
#'
#' @param df Original data frame (used to compute stratum weights)
#' @param theta Data frame with theta estimates by X stratum
#'
#' @return One-row data frame with averaged theta estimates
#'
#' @details
#' Computes theta_ATE = sum_x w_x * theta_x where w_x = n_x / n
#'
#' @keywords internal
compute_theta_ate <- function(df, theta, stratum_weights = NULL) {
  X_cols <- names(theta)[!startsWith(names(theta), "theta_Y")]

  if (length(X_cols)==0) {
    theta
  } else {
    wX <- if (!is.null(stratum_weights)) {
      dplyr::rename(stratum_weights, w = ".__w_sum__")
    } else {
      df |>
        count(across(all_of(X_cols)), name = "n") |>
        mutate(w = n / sum(n)) |>
        select(-n)
    }

    theta |>
      left_join(wX, by = X_cols) |>
      summarise(across(starts_with("theta_"), ~ sum(.x * .data$w, na.rm = TRUE)))
  }
}


#' Compute final treatment effect theta0
#'
#' Converts the theta vector (treatment effects for non-reference categories)
#' into a single scalar treatment effect estimate theta0.
#'
#' @param theta_ate One-row data frame with averaged theta estimates
#' @param y_levels Vector of outcome levels for Y
#' @param yK Reference category for Y
#'
#' @return Scalar theta0 (final treatment effect estimate)
#'
#' @details
#' Computes theta0 = sum_j y_j * theta_j + y_K * (1 - sum_j theta_j)
#'
#' This gives the expected value of Y under treatment.
#'
#' @keywords internal
compute_theta0 <- function(theta_ate, y_levels, yK) {
  yj_levels <- setdiff(y_levels, yK)

  # Extract theta values in order
  theta_vec <- as.numeric(theta_ate[1, paste0("theta_Y", yj_levels), drop = TRUE])

  theta0 <- sum((yj_levels - yK) * theta_vec, na.rm = TRUE)
  theta0
}


#' Compute RSV treatment effect estimate
#'
#' Main function that computes the RSV treatment effect estimate from a
#' combined data frame containing observations and predictions.
#'
#' @param df Data frame with Y, D, S_e, S_o, X (optional), and pred_* columns
#' @param theta_init Initial theta estimate (from compute_theta_init)
#' @param sigma2.lower Lower bound for sigma-squared stabilization (optional)
#' @param sigma2.quantile Quantile for sigma-squared lower bound (default 0.01)
#'
#' @return List with:
#' \describe{
#'   \item{coefficients}{Final treatment effect estimate (scalar, named "D")}
#'   \item{numerator}{Numerator of estimating equation}
#'   \item{denominator}{Denominator of estimating equation}
#'   \item{weights}{Efficient weights (H matrix)}
#'   \item{theta_init}{Initial theta estimate passed in}
#'   \item{n.obs}{Sample size in observational sample}
#'   \item{n.exp}{Sample size in experimental sample}
#'   \item{n.both}{Sample size in both samples}
#' }
#'
#' @details
#' This is the main computational function. It:
#' \enumerate{
#'   \item Computes final theta using efficient weights
#'   \item Averages across X strata if applicable
#'   \item Converts theta vector to scalar theta0
#'   \item Returns estimate with metadata
#' }
#'
#' @keywords internal
compute_estimate <- function(df, theta_init, y_levels, sigma2.lower = NULL, sigma2.quantile = 0.01,
                            pred_y_only = FALSE, weights = NULL) {
  # load("~/Documents/remoteoutcome_v2/data/pred_real_Ycons.rda")
  # df <- pred_real_Ycons |> mutate(pred_Y1 = pred_Y, pred_Y0 = 1-pred_Y, X1 = runif(n())>0.5)
  # theta_init <- compute_theta_init(df, X_cols = c("X1"), yK=0)
  # sigma2.lower = NULL
  # sigma2.quantile = 0.01
  
  # Convert logical to numeric for S_e and S_o
  df <- .ensure_numeric_indicators(df)

  # Sample sizes
  n.exp <- sum(df$S_e == 1)
  n.obs <- sum(df$S_o == 1)
  n.both <- sum(df$S_e == 1 & df$S_o == 1)

  
  # Compute final theta
  theta_df <- compute_theta(df, theta_init, y_levels = y_levels, sigma2.lower = sigma2.lower,
                            sigma2.quantile = sigma2.quantile, pred_y_only = pred_y_only,
                            weights = weights)
  theta <- theta_df$theta

  # Average across X strata; optional observation weights (e.g. clustered bootstrap)
  X_cols_est <- names(theta_init)[!startsWith(names(theta_init), "theta_Y")]
  strat_w <- if (!is.null(weights) && length(X_cols_est) > 0) {
    w_df <- df[, X_cols_est, drop = FALSE]
    w_df[[".__bb_w__"]] <- weights
    dplyr::group_by(w_df, dplyr::across(dplyr::all_of(X_cols_est))) |>
      dplyr::summarise(.__w_sum__ = sum(.data$`.__bb_w__`), .groups = "drop")
  } else NULL
  theta_ate <- compute_theta_ate(df, theta, stratum_weights = strat_w)

  # Convert to scalar
  yK <- .get_yK_value(y_levels, theta_init)
  theta0 <- compute_theta0(theta_ate, y_levels = y_levels, yK = yK)

  # ATE-specific relevance: (lambda' J^{-1} lambda)^{-1} pooled with pi weights
  denom_list <- theta_df$denominator$denominator
  pi_weights <- if (length(X_cols_est) > 0 && length(denom_list) > 1) {
    counts   <- dplyr::count(df, dplyr::across(dplyr::all_of(X_cols_est)))
    denom_joined <- dplyr::left_join(theta_df$denominator, counts, by = X_cols_est)
    denom_joined$n / sum(denom_joined$n)
  } else NULL
  relevance_ate <- compute_relevance_ate(denom_list, y_levels, yK, pi_weights)

  list(
    coefficients = structure(theta0, names = "D"),
    numerator = theta_df$numerator,
    denominator = theta_df$denominator,
    relevance_ate = relevance_ate,
    weights = theta_df$H,
    theta_init = theta_init,
    theta = theta,
    theta_ate = theta_ate,
    theta0 = theta0,
    n.obs = n.obs,
    n.exp = n.exp,
    n.both = n.both
  )
}