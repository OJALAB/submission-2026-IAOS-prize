library(data.table)

weighted_mean_zt <- function(x, w) sum(w * x) / sum(w)

weighted_var_zt <- function(x, w) {
  mu <- weighted_mean_zt(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

safe_mu_zt <- function(eta, eps = 1e-12) pmax(exp(pmin(eta, 20)), eps)

size_from_sigma_zt <- function(sigma, eps = 1e-12) 1 / pmax(sigma, eps)

p0_count <- function(mu, sigma = NULL, family = "NB", eps = 1e-8) {
  if (family == "NB") {
    size <- size_from_sigma_zt(sigma, eps = eps)
    return(pmin(dnbinom(0, mu = mu, size = size), 1 - eps))
  }
  pmin(dpois(0, lambda = mu), 1 - eps)
}

dzt_count <- function(y, mu, sigma = NULL, family = "NB", log = FALSE, eps = 1e-8) {
  p0 <- p0_count(mu = mu, sigma = sigma, family = family, eps = eps)
  if (family == "NB") {
    size <- size_from_sigma_zt(sigma, eps = eps)
    out <- dnbinom(y, mu = mu, size = size, log = FALSE) / pmax(1 - p0, eps)
  } else {
    out <- dpois(y, lambda = mu, log = FALSE) / pmax(1 - p0, eps)
  }
  out <- pmax(out, eps)
  if (log) return(log(out))
  out
}

pzt_ge_count <- function(C, mu, sigma = NULL, family = "NB", log = FALSE, eps = 1e-8) {
  p0 <- p0_count(mu = mu, sigma = sigma, family = family, eps = eps)
  if (family == "NB") {
    size <- size_from_sigma_zt(sigma, eps = eps)
    tail_prob <- pmax(1 - pnbinom(C - 1, mu = mu, size = size), eps)
  } else {
    tail_prob <- pmax(1 - ppois(C - 1, lambda = mu), eps)
  }
  out <- tail_prob / pmax(1 - p0, eps)
  out <- pmax(out, eps)
  if (log) return(log(out))
  out
}

ezt_mean <- function(mu, sigma = NULL, family = "NB", eps = 1e-8) {
  p0 <- p0_count(mu = mu, sigma = sigma, family = family, eps = eps)
  mu / pmax(1 - p0, eps)
}

ezt_tail_mean <- function(C, mu, sigma = NULL, family = "NB", eps = 1e-8) {
  if (family == "NB") {
    size <- size_from_sigma_zt(sigma, eps = eps)
    tail_prob <- pmax(1 - pnbinom(C - 1, mu = mu, size = size), eps)
    if (C > 1) {
      expected_less_C <- Reduce(`+`, lapply(seq_len(C - 1), function(v) v * dnbinom(v, mu = mu, size = size)))
    } else {
      expected_less_C <- rep(0, length(mu))
    }
    return(pmax(mu - expected_less_C, 0) / tail_prob)
  }
  tail_prob <- pmax(1 - ppois(C - 1, lambda = mu), eps)
  if (C > 1) {
    expected_less_C <- Reduce(`+`, lapply(seq_len(C - 1), function(v) v * dpois(v, lambda = mu)))
  } else {
    expected_less_C <- rep(0, length(mu))
  }
  pmax(mu - expected_less_C, 0) / tail_prob
}

fit_zt_cens <- function(formula,
                        data,
                        weights,
                        family = "NB",
                        start = NULL,
                        control = list()) {
  dep_vars <- all.vars(formula[[2]])
  y_name <- dep_vars[1]
  status_name <- dep_vars[2]

  data <- as.data.table(copy(data))
  data <- data[!is.na(data[[y_name]]) & !is.na(data[[status_name]])]

  y <- as.numeric(data[[y_name]])
  status <- as.numeric(data[[status_name]])

  expl_formula <- formula
  expl_formula[[2]] <- NULL

  expl_vars <- all.vars(expl_formula)
  factor_vars <- intersect(expl_vars, names(data)[vapply(data, is.factor, logical(1))])
  if (length(factor_vars) > 0) {
    xlevels <- lapply(data[, factor_vars, with = FALSE], levels)
  } else {
    xlevels <- list()
  }

  mf <- model.frame(expl_formula, data = data, xlev = xlevels, drop.unused.levels = FALSE)
  X_mat <- model.matrix(expl_formula, mf)

  if (missing(weights)) {
    w_vec <- rep(1, NROW(data))
  } else {
    w_vec <- as.numeric(data[[weights]])
  }
  w_vec[!is.finite(w_vec)] <- 0
  w_vec <- pmax(w_vec, 0)

  p <- NCOL(X_mat)
  eps <- 1e-8
  big_nll <- 1e100

  nll <- function(params) {
    beta <- params[seq_len(p)]
    eta <- as.numeric(X_mat %*% beta)
    if (any(!is.finite(eta))) return(big_nll)
    mu <- safe_mu_zt(eta, eps = eps)

    if (family == "NB") {
      log_sigma <- params[p + 1L]
      sigma <- exp(log_sigma)
      if (!is.finite(sigma) || sigma <= eps) return(big_nll)
      size <- size_from_sigma_zt(sigma, eps = eps)
      p_0 <- pmin(dnbinom(0, mu = mu, size = size), 1 - eps)

      ll_exact <- 0
      idx_exact <- which(status == 1)
      if (length(idx_exact) > 0) {
        log_p <- dnbinom(y[idx_exact], mu = mu[idx_exact], size = size, log = TRUE)
        ll_exact <- sum(w_vec[idx_exact] * (log_p - log(pmax(1 - p_0[idx_exact], eps))))
      }

      ll_cens <- 0
      idx_cens <- which(status == 0)
      if (length(idx_cens) > 0) {
        p_ge <- pmax(1 - pnbinom(y[idx_cens] - 1, mu = mu[idx_cens], size = size), eps)
        ll_cens <- sum(w_vec[idx_cens] * (log(p_ge) - log(pmax(1 - p_0[idx_cens], eps))))
      }

      out <- -(ll_exact + ll_cens)
      if (!is.finite(out)) return(big_nll)
      return(out)
    }

    p_0 <- pmin(dpois(0, lambda = mu), 1 - eps)

    ll_exact <- 0
    idx_exact <- which(status == 1)
    if (length(idx_exact) > 0) {
      log_p <- dpois(y[idx_exact], lambda = mu[idx_exact], log = TRUE)
      ll_exact <- sum(w_vec[idx_exact] * (log_p - log(pmax(1 - p_0[idx_exact], eps))))
    }

    ll_cens <- 0
    idx_cens <- which(status == 0)
    if (length(idx_cens) > 0) {
      p_ge <- pmax(1 - ppois(y[idx_cens] - 1, lambda = mu[idx_cens]), eps)
      ll_cens <- sum(w_vec[idx_cens] * (log(p_ge) - log(pmax(1 - p_0[idx_cens], eps))))
    }

    out <- -(ll_exact + ll_cens)
    if (!is.finite(out)) return(big_nll)
    out
  }

  idx_init <- which(status == 1 & is.finite(y) & w_vec > 0)
  if (length(idx_init) >= p) {
    init_beta <- tryCatch({
      init_model <- suppressWarnings(
        stats::glm.fit(
          x = X_mat[idx_init, , drop = FALSE],
          y = y[idx_init],
          weights = w_vec[idx_init],
          family = stats::poisson(link = "log")
        )
      )
      b <- as.numeric(init_model$coefficients)
      b[is.na(b) | !is.finite(b)] <- 0
      b
    }, error = function(e) rep(0, p))
  } else {
    init_beta <- rep(0, p)
  }
  names(init_beta) <- colnames(X_mat)

  if (family == "NB") {
    if (length(idx_init) >= 2) {
      mu0 <- weighted_mean_zt(y[idx_init], w_vec[idx_init])
      var0 <- weighted_var_zt(y[idx_init], w_vec[idx_init])
      sigma0 <- max((var0 - mu0) / pmax(mu0^2, eps), 1e-4)
    } else {
      sigma0 <- 0.25
    }
    default_init <- c(init_beta, log_sigma = log(sigma0))
  } else {
    default_init <- init_beta
  }

  init_params <- default_init
  if (!is.null(start)) {
    start_num <- as.numeric(start)
    start_names <- names(start)
    if (!is.null(start_names)) {
      common <- intersect(names(default_init), start_names)
      if (length(common) > 0) init_params[common] <- start_num[match(common, start_names)]
    } else if (length(start_num) == length(default_init)) {
      init_params <- start_num
      names(init_params) <- names(default_init)
    }
  }
  init_params <- as.numeric(init_params)
  names(init_params) <- names(default_init)

  default_control <- list(rel.tol = 1e-10, x.tol = 1e-10, eval.max = 2000, iter.max = 2000)
  default_control[names(control)] <- control

  lower_bounds <- rep(-20, length(init_params))
  upper_bounds <- rep(20, length(init_params))

  if (family == "NB") {
    upper_bounds[length(init_params)] <- 5
  }

  opt <- nlminb(
    start = init_params,
    objective = nll,
    control = default_control,
    lower = lower_bounds,
    upper = upper_bounds
  )

  beta_hat <- opt$par[seq_len(p)]
  names(beta_hat) <- colnames(X_mat)

  res <- list(
    coefficients = beta_hat,
    family = family,
    terms = expl_formula,
    xlevels = xlevels,
    convergence = opt$convergence,
    objective = opt$objective,
    message = opt$message,
    iterations = opt$iterations,
    nobs = NROW(data),
    weights_sum = sum(w_vec)
  )
  if (family == "NB") res$sigma <- exp(opt$par[p + 1L])

  class(res) <- "zt_cens"
  res
}

coef.zt_cens <- function(object, ...) {
  if (identical(object$family, "NB")) {
    c(object$coefficients, log_sigma = log(object$sigma))
  } else {
    object$coefficients
  }
}

predict.zt_cens <- function(object, newdata, what = c("mu", "sigma")) {
  what <- what[1]
  if (what == "sigma") {
    if (identical(object$family, "NB")) return(rep(object$sigma, NROW(newdata)))
    return(rep(NA_real_, NROW(newdata)))
  }
  mf <- model.frame(object$terms, newdata, xlev = object$xlevels, drop.unused.levels = FALSE)
  X_mat <- model.matrix(object$terms, mf)
  beta <- object$coefficients
  eta <- as.numeric(X_mat %*% beta)
  safe_mu_zt(eta)
}

generate_labels <- function(df, Pi) {
  df[, code_ml := {
    g <- as.character(code[1])
    as.character(
      sample(
        x = colnames(Pi),
        size = .N,
        replace = TRUE,
        prob = Pi[g, ]
      )
    )
  }, by = .(code)]

  df[, code_ml := factor(code_ml, levels = levels(df[["code"]]))]

  df
}

estimate_Pi <- function(df) {
  count_table <- table(df[["code"]], df[["code_ml"]])
  Pi_hat <- prop.table(count_table, margin = 1)
  as.matrix(Pi_hat)
}

baffour_em <- function(df, latent_classes = 2, tol = 1e-4, max_iter = 1000) {
  grid <- CJ(
    code = unique(df[["code"]]),
    I1 = factor(c("0", "1"), levels = c("0", "1")),
    I2 = factor(c("0", "1"), levels = c("0", "1")),
    I3 = factor(c("0", "1"), levels = c("0", "1"))
  )

  data_full <- merge(
    grid,
    df[, .(code, I1, I2, I3, count)],
    by = c("code", "I1", "I2", "I3"),
    all.x = TRUE
  )

  data_full[is.na(count), count := 0]
  data_full[, is_missing_cell := (I1 == "0" & I2 == "0" & I3 == "0")]

  data_em <- data_full[rep(1:.N, each = latent_classes)]
  data_em[, X := factor(rep(1:latent_classes, times = nrow(data_full)))]

  random_split <- runif(nrow(data_em) / latent_classes, 0.45, 0.55)
  data_em[X == "1", em_count := count * random_split]
  data_em[X == "2", em_count := count * (1 - random_split)]

  data_em[is_missing_cell == TRUE, em_count := 1]

  last_coef <- NULL
  converged <- FALSE
  iter <- 0

  while (!converged && iter < max_iter) {
    iter <- iter + 1

    model <- glm(
      em_count ~ (I1 + I2 + I3) * X + code * X + I1:I2 + I1:I3,
      family = poisson(),
      data = data_em
    )

    data_em[, mu_pred := fitted(model)]
    current_coef <- coef(model)

    if (!is.null(last_coef)) {
      diff <- max(abs(current_coef - last_coef))
      if (diff < tol) converged <- TRUE
    }
    last_coef <- current_coef

    data_em[, mu_total := sum(mu_pred), by = .(code, I1, I2, I3)]

    data_em[is_missing_cell == FALSE, em_count := count / mu_total * mu_pred]
    data_em[is_missing_cell == TRUE, em_count := mu_pred]
  }

  missing_preds <- data_em[
    is_missing_cell == TRUE,
    .(n_pred_lca = sum(em_count)),
    by = .(code)
  ]

  missing_preds
}

baffour_misclass_full <- function(df_obs,
                                  Pi,
                                  latent_classes = 2,
                                  C = 2,
                                  vac_family = "NB",
                                  tol = 1e-5,
                                  max_iter = 100,
                                  tau_floor = 1e-10,
                                  trace = FALSE) {
  obs_levels <- levels(df_obs[["code"]])

  Pi_df <- as.data.table(as.table(Pi))
  setnames(Pi_df, c("code", "code_obs", "pi_prob"))
  Pi_df[, code := factor(code, levels = obs_levels)]
  Pi_df[, code_obs := factor(code_obs, levels = obs_levels)]

  df_obs <- as.data.table(copy(df_obs))
  df_obs[, id := .I]
  df_obs[, `:=`(
    I1 = factor(I1, levels = c(0, 1)),
    I2 = factor(I2, levels = c(0, 1)),
    I3 = factor(I3, levels = c(0, 1)),
    code_obs = factor(code_obs, levels = obs_levels)
  )]

  true_grid <- CJ(
    id = df_obs[["id"]],
    code = factor(levels(Pi_df[["code"]]), levels = levels(Pi_df[["code"]])),
    X = factor(seq_len(latent_classes))
  )

  E_matrix <- merge(
    true_grid,
    df_obs[, .(id, I1, I2, I3, code_obs, vac_obs, v_status)],
    by = "id",
    all.x = TRUE,
    sort = FALSE
  )
  E_matrix <- merge(E_matrix, Pi_df, by = c("code", "code_obs"), all.x = TRUE, sort = FALSE)
  E_matrix[is.na(pi_prob), pi_prob := 0]
  E_matrix[I1 == "1", pi_prob := fifelse(code == code_obs, 1, 0)]
  E_matrix[, row_id := .I]

  split_dt <- unique(E_matrix[, .(id, code)])
  split_dt[, split_x := runif(.N, 0.35, 0.65)]
  E_matrix <- merge(E_matrix, split_dt, by = c("id", "code"), all.x = TRUE, sort = FALSE)
  E_matrix[, tau := pi_prob * fifelse(X == levels(X)[1], split_x, 1 - split_x)]
  E_matrix[!is.finite(tau), tau := 0]
  E_matrix[, tau := { s <- sum(tau); if (s > 0) tau / s else rep(1 / .N, .N) }, by = id]

  data_em <- CJ(
    code = factor(levels(Pi_df[["code"]]), levels = levels(Pi_df[["code"]])),
    I1 = factor(c("0", "1"), levels = c("0", "1")),
    I2 = factor(c("0", "1"), levels = c("0", "1")),
    I3 = factor(c("0", "1"), levels = c("0", "1")),
    X = factor(seq_len(latent_classes))
  )
  data_em[, is_missing_cell := (I1 == "0" & I2 == "0" & I3 == "0")]

  agg_initial <- E_matrix[, .(em_count = sum(tau)), by = .(I1, I2, I3, code, X)]
  data_em <- merge(data_em, agg_initial, by = c("I1", "I2", "I3", "code", "X"), all.x = TRUE, sort = FALSE)
  data_em[is.na(em_count), em_count := 0]
  data_em[is_missing_cell == TRUE, em_count := 1]

  E_matrix[, status_num := fifelse(v_status == "censored", 0, 1)]

  converged <- FALSE
  iter <- 0L
  tau_old <- NULL
  mu_crc_old <- NULL
  mu_vac_old <- NULL
  last_vac_start <- NULL
  last_delta <- Inf

  while (!converged && iter < max_iter) {
    iter <- iter + 1L
    if (isTRUE(trace)) message(sprintf("baffour_misclass_full: iter = %d", iter))

    model_crc <- tryCatch(
      suppressWarnings(glm(
        em_count ~ (I1 + I2 + I3) * X + code * X + I1:I2 + I1:I3,
        family = poisson(),
        data = data_em,
        method = brglm2::brglmFit
      )),
      error = function(e) suppressWarnings(glm(
        em_count ~ (I1 + I2 + I3) * X + code * X + I1:I2 + I1:I3,
        family = poisson(),
        data = data_em
      ))
    )
    data_em[, mu_pred_crc := pmax(as.numeric(fitted(model_crc)), tau_floor)]

    E_fit <- E_matrix[v_status != "missing", .(vac_obs, status_num, code, w = pmax(tau, 0))]
    E_fit[, code := factor(code, levels = obs_levels)]
    E_fit <- E_fit[, .(w = sum(w)), by = .(vac_obs, status_num, code)]
    E_fit[, code := factor(code, levels = obs_levels)]

    model_vac <- fit_zt_cens(
      formula = Surv(vac_obs, status_num) ~ code,
      data = E_fit,
      weights = "w",
      family = vac_family,
      start = last_vac_start
    )
    last_vac_start <- coef(model_vac)

    E_pred <- copy(E_matrix)
    E_pred[, code := factor(code, levels = obs_levels)]
    E_matrix[, mu_vac := predict.zt_cens(model_vac, newdata = E_pred, what = "mu")]
    E_matrix[, sigma_vac := predict.zt_cens(model_vac, newdata = E_pred, what = "sigma")]

    mu_to_join <- data_em[, .(I1, I2, I3, code, X, mu_pred_crc)]
    setkey(mu_to_join, I1, I2, I3, code, X)
    setkey(E_matrix, I1, I2, I3, code, X)
    E_matrix[mu_to_join, mu_pred_crc := i.mu_pred_crc]
    setorder(E_matrix, row_id)

    E_matrix[v_status == "exact", L_V := dzt_count(y = vac_obs, mu = mu_vac, sigma = sigma_vac, family = vac_family)]
    E_matrix[v_status == "censored", L_V := pzt_ge_count(C = C, mu = mu_vac, sigma = sigma_vac, family = vac_family)]
    E_matrix[v_status == "missing", L_V := 1]

    E_matrix[, log_num := log(pmax(mu_pred_crc, 1e-300)) + log(pmax(pi_prob, 1e-300)) + log(pmax(L_V, 1e-300))]
    E_matrix[, log_num_max := max(log_num), by = id]
    E_matrix[, num_stab := exp(log_num - log_num_max)]
    E_matrix[, den_stab := pmax(sum(num_stab), 1e-300), by = id]
    E_matrix[, tau := num_stab / den_stab]

    agg_n_hat <- E_matrix[, .(new_em_count = sum(tau)), by = .(I1, I2, I3, code, X)]
    data_em_no_count <- copy(data_em)
    data_em_no_count[, em_count := NULL]
    data_em <- merge(
      data_em_no_count,
      agg_n_hat,
      by = c("I1", "I2", "I3", "code", "X"),
      all.x = TRUE,
      sort = FALSE
    )
    setnames(data_em, "new_em_count", "em_count")
    data_em[is.na(em_count), em_count := 0]
    data_em[is_missing_cell == TRUE, em_count := mu_pred_crc]

    delta_tau <- if (is.null(tau_old)) Inf else max(abs(E_matrix$tau - tau_old))
    delta_crc <- if (is.null(mu_crc_old)) Inf else max(abs(log(pmax(data_em$mu_pred_crc, tau_floor)) - log(pmax(mu_crc_old, tau_floor))))
    delta_vac <- if (is.null(mu_vac_old)) Inf else max(abs(log(pmax(E_matrix$mu_vac, tau_floor)) - log(pmax(mu_vac_old, tau_floor))))
    last_delta <- max(delta_tau, delta_crc, delta_vac)
    converged <- is.finite(last_delta) && (last_delta < tol)

    tau_old <- copy(E_matrix$tau)
    mu_crc_old <- copy(data_em$mu_pred_crc)
    mu_vac_old <- copy(E_matrix$mu_vac)
  }

  E_matrix[, nu_hat := NA_real_]
  E_matrix[v_status == "missing", nu_hat := ezt_mean(mu = mu_vac, sigma = sigma_vac, family = vac_family)]
  E_matrix[v_status == "exact", nu_hat := vac_obs]
  E_matrix[v_status == "censored", nu_hat := ezt_tail_mean(C = C, mu = mu_vac, sigma = sigma_vac, family = vac_family)]

  n_preds <- data_em[, .(n_pred_lca = sum(em_count)), by = .(code)]
  vac_preds_obs <- E_matrix[, .(vac_obs = sum(tau * nu_hat)), by = .(code)]

  hidden_pop <- data_em[is_missing_cell == TRUE, .(code, X, em_count)]
  hidden_newdata_df <- data.frame(
    code = factor(hidden_pop[["code"]], levels = obs_levels),
    vac_obs = 1,
    status_num = 1
  )
  pred_mu_hidden <- predict.zt_cens(model_vac, newdata = hidden_newdata_df, what = "mu")
  pred_sigma_hidden <- predict.zt_cens(model_vac, newdata = hidden_newdata_df, what = "sigma")
  hidden_pop[, mu_vac_hidden := ezt_mean(mu = pred_mu_hidden, sigma = pred_sigma_hidden, family = vac_family)]

  vac_preds_hidden <- hidden_pop[, .(vac_hidden = sum(em_count * mu_vac_hidden)), by = .(code)]
  vac_preds <- merge(vac_preds_obs, vac_preds_hidden, by = "code", all = TRUE, sort = FALSE)
  vac_preds[is.na(vac_obs), vac_obs := 0]
  vac_preds[is.na(vac_hidden), vac_hidden := 0]
  vac_preds[, total_vac_pred := vac_obs + vac_hidden]

  list(
    n_preds = n_preds,
    vac_preds = vac_preds,
    iter = iter,
    converged = converged,
    last_delta = last_delta,
    vac_family = vac_family
  )
}

sim_vac <- function(df,
                    p_easy,
                    p_hard,
                    prob_hard_vec,
                    cens_frac,
                    missing_frac,
                    Pi,
                    val_sample_size,
                    C = 2,
                    vac_family = "NB",
                    use_estimated_Pi = FALSE,
                    tol_em = 1e-5,
                    max_iter_em = 100,
                    trace_em = FALSE) {
  true_total_vac <- sum(df[["vac"]])
  true_n <- NROW(df)

  df <- as.data.table(copy(df))
  df[, code := relevel(code, ref = "2")]

  true_vac_by_code <- df[, .(true_vac_code = sum(vac)), by = .(code)]

  df <- generate_labels(df, Pi)

  Pi_est <- Pi
  if (isTRUE(use_estimated_Pi)) {
    code_dist <- prop.table(table(df[["code"]]))
    df_val <- data.table(
      code = factor(
        sample(
          x = names(code_dist),
          size = val_sample_size,
          replace = TRUE,
          prob = as.numeric(code_dist)
        ),
        levels = levels(df[["code"]])
      )
    )
    df_val <- generate_labels(df_val, Pi)
    Pi_est <- estimate_Pi(df_val)
  }

  params_probs <- data.table(code = factor(1:9, levels = levels(df[["code"]])), prob_is_hard = prob_hard_vec)
  df <- merge(df, params_probs, by = "code", sort = FALSE)
  df[, is_hard := rbinom(.N, 1, prob_is_hard)]
  df[, `:=`(
    prob_1 = fifelse(is_hard == 1, p_hard[1], p_easy[1]),
    prob_2 = fifelse(is_hard == 1, p_hard[2], p_easy[2]),
    prob_3 = fifelse(is_hard == 1, p_hard[3], p_easy[3])
  )]
  df[, I1 := rbinom(.N, 1, prob_1)]
  df[, prob_2 := pmax(prob_2 - (I1 * 0.05), 0)]
  df[, I2 := rbinom(.N, 1, prob_2)]
  df[, prob_3 := pmax(prob_3 - (I1 * 0.05), 0)]
  df[, I3 := rbinom(.N, 1, prob_3)]
  df[, c("prob_is_hard", "is_hard", "prob_1", "prob_2", "prob_3") := NULL]

  df_obs <- copy(df[I1 == 1 | I2 == 1 | I3 == 1])
  df_obs[, code_obs := code_ml]
  df_obs[I1 == 1, code_obs := code]
  df_obs[, code_obs := factor(code_obs, levels = levels(df[["code"]]))]

  df_obs[, vac := as.numeric(vac)]
  df_obs[, vac_true := vac]
  df_obs[, vac_obs := vac]
  df_obs[, v_status := "exact"]

  df_obs[I1 == 0, rand_val := runif(.N)]
  df_obs[I1 == 0 & rand_val <= missing_frac, v_status := "missing"]
  df_obs[v_status == "missing", vac_obs := NA_real_]
  df_obs[I1 == 0 & vac_true >= C & rand_val > missing_frac & rand_val <= (missing_frac + cens_frac), v_status := "censored"]
  df_obs[v_status == "censored", vac_obs := as.numeric(C)]
  df_obs[, rand_val := NULL]

  all_combinations_naive <- CJ(
    code_obs = factor(levels(df[["code"]]), levels = levels(df[["code"]])),
    I1 = factor(c("0", "1"), levels = c("0", "1")),
    I2 = factor(c("0", "1"), levels = c("0", "1")),
    I3 = factor(c("0", "1"), levels = c("0", "1"))
  )
  all_combinations_naive <- all_combinations_naive[!(I1 == "0" & I2 == "0" & I3 == "0")]

  df_agg_naive <- df_obs[, .(count = .N), by = .(
    code_obs,
    I1 = factor(as.character(I1), levels = c("0", "1")),
    I2 = factor(as.character(I2), levels = c("0", "1")),
    I3 = factor(as.character(I3), levels = c("0", "1"))
  )]
  df_agg_naive <- df_agg_naive[all_combinations_naive, on = .(code_obs, I1, I2, I3)]
  df_agg_naive[is.na(count), count := 0]
  setnames(df_agg_naive, "code_obs", "code")

  df_hidden <- CJ(
    code = factor(levels(df[["code"]]), levels = levels(df[["code"]])),
    I1 = factor("0", levels = c("0", "1")),
    I2 = factor("0", levels = c("0", "1")),
    I3 = factor("0", levels = c("0", "1"))
  )

  model_lca <- baffour_em(df_agg_naive)

  model_misclass_base <- baffour_misclass_full(
    df_obs = df_obs,
    Pi = Pi,
    C = C,
    vac_family = vac_family,
    tol = tol_em,
    max_iter = max_iter_em,
    trace = trace_em
  )

  if (isTRUE(use_estimated_Pi)) {
    model_misclass <- baffour_misclass_full(
      df_obs = df_obs,
      Pi = Pi_est,
      C = C,
      vac_family = vac_family,
      tol = tol_em,
      max_iter = max_iter_em,
      trace = trace_em
    )
  } else {
    model_misclass <- model_misclass_base
  }

  df_obs[, status := 1L]
  df_obs[v_status == "censored", status := 0L]
  df_obs_to_reg <- df_obs[v_status != "missing"]
  df_obs_to_reg[, code_obs := factor(code_obs, levels = levels(df[["code"]]))]

  model_reg <- fit_zt_cens(
    formula = Surv(vac_obs, status) ~ code_obs,
    data = df_obs_to_reg,
    family = vac_family
  )

  idx_cens <- which(df_obs$status == 0L)
  if (length(idx_cens) > 0) {
    tmp <- copy(df_obs[idx_cens])
    tmp[, code_obs := factor(code_obs, levels = levels(df[["code"]]))]
    mu_pred_cens <- predict.zt_cens(model_reg, newdata = tmp, what = "mu")
    sigma_pred_cens <- predict.zt_cens(model_reg, newdata = tmp, what = "sigma")
    df_obs[idx_cens, vac_obs := ezt_tail_mean(C = C, mu = mu_pred_cens, sigma = sigma_pred_cens, family = vac_family)]
  }

  idx_missing <- which(df_obs$v_status == "missing")
  if (length(idx_missing) > 0) {
    tmp <- copy(df_obs[idx_missing])
    tmp[, code_obs := factor(code_obs, levels = levels(df[["code"]]))]
    mu_pred_missing <- predict.zt_cens(model_reg, newdata = tmp, what = "mu")
    sigma_pred_missing <- predict.zt_cens(model_reg, newdata = tmp, what = "sigma")
    df_obs[idx_missing, vac_obs := ezt_mean(mu = mu_pred_missing, sigma = sigma_pred_missing, family = vac_family)]
  }

  predicted_means <- df_obs[, .(mean_pred = mean(vac_obs, na.rm = TRUE)), by = .(code_obs)]
  fallback_mean <- mean(df_obs[["vac_obs"]], na.rm = TRUE)
  predicted_means[is.nan(mean_pred), mean_pred := fallback_mean]
  setnames(predicted_means, "code_obs", "code")

  df_hidden <- merge(df_hidden, model_lca, by = "code", all.x = TRUE, sort = FALSE)
  df_hidden <- merge(df_hidden, predicted_means, by = "code", all.x = TRUE, sort = FALSE)
  df_hidden[is.na(mean_pred), mean_pred := fallback_mean]
  df_hidden[, vac_pred := mean_pred * n_pred_lca]

  est_n_naive <- NROW(df_obs) + sum(df_hidden[["n_pred_lca"]], na.rm = TRUE)
  est_vac_naive <- sum(df_obs[["vac_obs"]], na.rm = TRUE) + sum(df_hidden[["vac_pred"]], na.rm = TRUE)
  est_n_base <- sum(model_misclass_base$n_preds[["n_pred_lca"]], na.rm = TRUE)
  est_vac_base <- sum(model_misclass_base$vac_preds[["total_vac_pred"]], na.rm = TRUE)
  est_n <- sum(model_misclass$n_preds[["n_pred_lca"]], na.rm = TRUE)
  est_vac <- sum(model_misclass$vac_preds[["total_vac_pred"]], na.rm = TRUE)

  vac_obs_naive_by_code <- df_obs[, .(vac_obs_naive = sum(vac_obs, na.rm = TRUE)), by = .(code_obs)]
  setnames(vac_obs_naive_by_code, "code_obs", "code")

  vac_hidden_naive_by_code <- df_hidden[, .(vac_hidden_naive = sum(vac_pred, na.rm = TRUE)), by = .(code)]

  est_vac_naive_by_code <- merge(
    vac_obs_naive_by_code,
    vac_hidden_naive_by_code,
    by = "code",
    all = TRUE,
    sort = FALSE
  )
  est_vac_naive_by_code[is.na(vac_obs_naive), vac_obs_naive := 0]
  est_vac_naive_by_code[is.na(vac_hidden_naive), vac_hidden_naive := 0]
  est_vac_naive_by_code[, est_vac_naive_code := vac_obs_naive + vac_hidden_naive]

  est_vac_base_by_code <- copy(
    model_misclass_base$vac_preds[, .(code, est_vac_base_code = total_vac_pred)]
  )

  est_vac_by_code <- copy(
    model_misclass$vac_preds[, .(code, est_vac_code = total_vac_pred)]
  )

  by_code <- Reduce(
    function(x, y) merge(x, y, by = "code", all = TRUE, sort = FALSE),
    list(
      true_vac_by_code,
      est_vac_naive_by_code[, .(code, est_vac_naive_code)],
      est_vac_base_by_code,
      est_vac_by_code
    )
  )

  num_cols <- setdiff(names(by_code), "code")
  for (nm in num_cols) {
    by_code[is.na(get(nm)), (nm) := 0]
  }

  by_code[, code_num := as.integer(as.character(code))]
  setorder(by_code, code_num)
  by_code[, code_num := NULL]

  code_chr <- as.character(by_code$code)
  by_code_vec <- c(
    setNames(by_code$true_vac_code, paste0("true_vac_code_", code_chr)),
    setNames(by_code$est_vac_naive_code, paste0("est_vac_naive_code_", code_chr)),
    setNames(by_code$est_vac_base_code, paste0("est_vac_base_code_", code_chr)),
    setNames(by_code$est_vac_code, paste0("est_vac_code_", code_chr))
  )

  c(
    true_n = true_n,
    true_total_vac = true_total_vac,
    est_n_naive = est_n_naive,
    est_vac_naive = est_vac_naive,
    est_n_base = est_n_base,
    est_vac_base = est_vac_base,
    est_n = est_n,
    est_vac = est_vac,
    by_code_vec
  )
}

eval <- function(df_results) {
  df_results <- as.data.table(df_results)
  method_map <- c(est_vac_naive = "M1", est_vac_base = "M2", est_vac = "M3")

  metrics_list <- list()
  for (est in c("est_vac_naive", "est_vac_base", "est_vac")) {
    error <- df_results[[est]] - df_results[["true_total_vac"]]
    metrics_list[[est]] <- data.table(
      method = method_map[[est]],
      mean_value = mean(df_results[[est]]),
      mean_bias = mean(error),
      mean_rel_bias = 100 * mean(error / df_results[["true_total_vac"]]),
      rmse = sqrt(mean(error^2))
    )
  }
  metrics_vac <- rbindlist(metrics_list)

  codes <- 1:9

  metrics_list <- list()
  for (code in codes) {
    true_col <- paste0("true_vac_code_", code)

    for (est in c("est_vac_naive_code_", "est_vac_base_code_", "est_vac_code_")) {
      est_col <- paste0(est, code)
      error <- df_results[[est_col]] - df_results[[true_col]]
      est_name <- substr(est, 1, nchar(est) - 6)

      metrics_list[[est_col]] <- data.table(
        code = code,
        method = method_map[[est_name]],
        mean_value = mean(df_results[[est_col]]),
        mean_bias = mean(error),
        mean_rel_bias = 100 * mean(error / df_results[[true_col]]),
        rmse = sqrt(mean(error^2))
      )
    }
  }
  metrics_vac_code <- rbindlist(metrics_list)

  list(metrics_vac, metrics_vac_code)
}
