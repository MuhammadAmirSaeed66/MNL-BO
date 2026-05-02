############################################################
# MNL-based Bayesian Optimization for a Categorical Benchmark
#
# Search space: 5 categorical variables, 4 levels each
# Objective: main effects, pairwise interactions, and Gaussian noise
# Surrogate: multinomial logit / conditional logit via mlogit
# Uncertainty: asymptotic covariance with delta method
# Acquisition: probability-weighted EI and UCB
############################################################

required_packages <- c(
  "mlogit", "data.table", "ggplot2", "plotly", "viridis",
  "Rtsne", "igraph", "visNetwork"
)

install_missing_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) install.packages(missing_packages)
}

install_missing_packages(required_packages)

suppressPackageStartupMessages({
  library(mlogit)
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(viridis)
  library(Rtsne)
  library(igraph)
  library(visNetwork)
})

set.seed(123)

levels4 <- 1:4
n_variables <- 5
variable_names <- paste0("c", seq_len(n_variables))

X_all <- expand.grid(
  c1 = levels4,
  c2 = levels4,
  c3 = levels4,
  c4 = levels4,
  c5 = levels4
)

for (j in seq_len(n_variables)) {
  X_all[[j]] <- factor(X_all[[j]], levels = levels4)
}

n_configs <- nrow(X_all)
interaction_pairs <- combn(seq_len(n_variables), 2, simplify = FALSE)

X_mm <- model.matrix(~ (c1 + c2 + c3 + c4 + c5)^2 - 1, data = X_all)
n_features <- ncol(X_mm)

make_theta <- function(n_variables = 5, n_levels = 4) {
  theta_main <- lapply(seq_len(n_variables), function(j) {
    rnorm(n_levels, mean = 0, sd = 1)
  })
  
  pairs <- combn(seq_len(n_variables), 2, simplify = FALSE)
  
  theta_pair <- setNames(
    vector("list", length(pairs)),
    vapply(pairs, function(pk) paste0("c", pk[1], "_c", pk[2]), character(1))
  )
  
  for (nm in names(theta_pair)) {
    theta_pair[[nm]] <- matrix(rnorm(n_levels^2, mean = 0, sd = 0.6), nrow = n_levels)
  }
  
  list(theta_main = theta_main, theta_pair = theta_pair)
}

boost_landscape <- function(theta, c_star, boost_best = 6.0, boost_near = 4.5) {
  for (j in seq_along(c_star)) {
    theta$theta_main[[j]][c_star[j]] <-
      theta$theta_main[[j]][c_star[j]] + boost_best / length(c_star)
  }
  
  pairs <- combn(seq_along(c_star), 2, simplify = FALSE)
  
  for (pk in pairs) {
    nm <- paste0("c", pk[1], "_c", pk[2])
    theta$theta_pair[[nm]][c_star[pk[1]], c_star[pk[2]]] <-
      theta$theta_pair[[nm]][c_star[pk[1]], c_star[pk[2]]] + boost_best / length(pairs)
  }
  
  near_configs <- list(
    c(c_star[1], c_star[2], c_star[3], c_star[4], 2),
    c(c_star[1], c_star[2], 4,         c_star[4], c_star[5]),
    c(2,         c_star[2], c_star[3], c_star[4], c_star[5])
  )
  
  for (cn in near_configs) {
    for (j in seq_along(cn)) {
      theta$theta_main[[j]][cn[j]] <-
        theta$theta_main[[j]][cn[j]] + boost_near / length(cn)
    }
    
    for (pk in pairs) {
      nm <- paste0("c", pk[1], "_c", pk[2])
      theta$theta_pair[[nm]][cn[pk[1]], cn[pk[2]]] <-
        theta$theta_pair[[nm]][cn[pk[1]], cn[pk[2]]] + boost_near / (length(pairs) + 2)
    }
  }
  
  theta
}

f_obj <- function(c_row, theta, noise_sd = 0.05) {
  c_int <- as.integer(c_row)
  value <- 0
  
  for (j in seq_along(c_int)) {
    value <- value + theta$theta_main[[j]][c_int[j]]
  }
  
  pairs <- combn(seq_along(c_int), 2, simplify = FALSE)
  
  for (pk in pairs) {
    nm <- paste0("c", pk[1], "_c", pk[2])
    value <- value + theta$theta_pair[[nm]][c_int[pk[1]], c_int[pk[2]]]
  }
  
  value + rnorm(1, mean = 0, sd = noise_sd)
}

theta <- make_theta(n_variables = n_variables, n_levels = length(levels4))

c_star <- c(1, 2, 3, 4, 1)
names(c_star) <- variable_names

theta <- boost_landscape(theta, c_star)

y_true <- vapply(
  seq_len(n_configs),
  function(i) f_obj(X_all[i, ], theta, noise_sd = 0),
  numeric(1)
)

oracle_best_idx <- which.max(y_true)
oracle_best_cfg <- X_all[oracle_best_idx, , drop = FALSE]
oracle_best_val <- y_true[oracle_best_idx]

cat("Oracle best index:", oracle_best_idx, "\n")
cat("Oracle best configuration:\n")
print(oracle_best_cfg)
cat("Oracle best value:", oracle_best_val, "\n\n")

make_choice_long <- function(eval_ids, y_obs, chid) {
  chosen_id <- eval_ids[which.max(y_obs[as.character(eval_ids)])]
  
  choice_df <- data.frame(
    chid = chid,
    alt = as.character(eval_ids),
    choice = as.integer(eval_ids == chosen_id),
    stringsAsFactors = FALSE
  )
  
  cbind(choice_df, X_mm[eval_ids, , drop = FALSE])
}

fit_mnl_surrogate <- function(choice_long_df) {
  mdat <- mlogit.data(
    choice_long_df,
    choice = "choice",
    shape = "long",
    chid.var = "chid",
    alt.var = "alt"
  )
  
  formula_mnl <- as.formula(
    paste0("choice ~ 0 + ", paste(colnames(X_mm), collapse = " + "))
  )
  
  tryCatch(
    mlogit(formula_mnl, data = mdat),
    error = function(e) NULL
  )
}

predict_mnl_all <- function(model) {
  if (is.null(model)) {
    return(list(
      mu = rep(0, n_configs),
      sigma = rep(1, n_configs),
      pi = rep(1 / n_configs, n_configs)
    ))
  }
  
  beta <- stats::coef(model)
  covariance <- tryCatch(stats::vcov(model), error = function(e) NULL)
  
  beta_vec <- rep(0, n_features)
  names(beta_vec) <- colnames(X_mm)
  beta_vec[names(beta)] <- beta
  
  mu <- as.numeric(X_mm %*% beta_vec)
  
  if (is.null(covariance)) {
    sigma <- rep(1, n_configs)
  } else {
    covariance_full <- matrix(
      0,
      nrow = n_features,
      ncol = n_features,
      dimnames = list(colnames(X_mm), colnames(X_mm))
    )
    
    covariance_full[rownames(covariance), colnames(covariance)] <- covariance
    
    XV <- X_mm %*% covariance_full
    variance <- rowSums(XV * X_mm)
    sigma <- sqrt(pmax(variance, 1e-12))
  }
  
  exp_mu <- exp(mu - max(mu))
  pi <- exp_mu / sum(exp_mu)
  
  list(mu = mu, sigma = sigma, pi = pi)
}

acq_ei_mnl <- function(mu, pi, best_y) {
  pmax(mu - best_y, 0) * pi
}

acq_ucb_mnl <- function(mu, sigma, pi, kappa = 2.0) {
  mu + kappa * sigma * pi
}

run_mnl_bo <- function(
    T_budget = 60,
    n0 = 10,
    noise_sd = 0.05,
    kappa = 2.0,
    acq = c("EI", "UCB")
) {
  acq <- match.arg(acq)
  
  eval_ids <- sample(seq_len(n_configs), size = n0, replace = FALSE)
  y_obs <- vapply(
    eval_ids,
    function(i) f_obj(X_all[i, ], theta, noise_sd = noise_sd),
    numeric(1)
  )
  names(y_obs) <- as.character(eval_ids)
  
  history <- data.frame(
    t = seq_len(n0),
    chosen_id = eval_ids,
    y = y_obs,
    best_y = cummax(y_obs)
  )
  
  choice_long <- NULL
  
  for (t in seq.int(n0 + 1, T_budget)) {
    choice_long_t <- make_choice_long(eval_ids, y_obs, chid = t - n0)
    choice_long <- if (is.null(choice_long)) {
      choice_long_t
    } else {
      rbind(choice_long, choice_long_t)
    }
    
    model <- fit_mnl_surrogate(choice_long)
    pred <- predict_mnl_all(model)
    
    unevaluated_ids <- setdiff(seq_len(n_configs), eval_ids)
    if (length(unevaluated_ids) == 0) break
    
    acquisition <- switch(
      acq,
      EI = acq_ei_mnl(pred$mu, pred$pi, max(y_obs)),
      UCB = acq_ucb_mnl(pred$mu, pred$sigma, pred$pi, kappa = kappa)
    )
    
    next_id <- unevaluated_ids[which.max(acquisition[unevaluated_ids])]
    y_next <- f_obj(X_all[next_id, ], theta, noise_sd = noise_sd)
    
    eval_ids <- c(eval_ids, next_id)
    y_obs[as.character(next_id)] <- y_next
    
    history <- rbind(history, data.frame(
      t = t,
      chosen_id = next_id,
      y = y_next,
      best_y = max(y_obs)
    ))
  }
  
  best_id <- as.integer(names(y_obs)[which.max(y_obs)])
  
  list(
    history = history,
    best_id = best_id,
    best_cfg = X_all[best_id, , drop = FALSE],
    best_y = max(y_obs),
    oracle_best_id = oracle_best_idx,
    oracle_best_cfg = oracle_best_cfg,
    oracle_best_y = oracle_best_val
  )
}

res_EI <- run_mnl_bo(T_budget = 60, n0 = 10, noise_sd = 0.05, acq = "EI", kappa = 2.0)
res_UCB <- run_mnl_bo(T_budget = 60, n0 = 10, noise_sd = 0.05, acq = "UCB", kappa = 2.0)

print_run_result <- function(result, label) {
  cat("\n===", label, "===\n")
  cat("Best found objective:", result$best_y, "\n")
  cat("Best found configuration:\n")
  print(result$best_cfg)
  cat("Oracle best objective:", result$oracle_best_y, "\n")
  cat("Oracle best configuration:\n")
  print(result$oracle_best_cfg)
}

print_run_result(res_EI, "EI-MNL BO")
print_run_result(res_UCB, "UCB-MNL BO")

cfg_key <- function(df_row) paste(as.integer(df_row), collapse = "-")

safe_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

make_hist_dt <- function(result, method_name) {
  h <- as.data.table(result$history)
  h[, method := method_name]
  h[, best_y_lag := shift(best_y, type = "lag"), by = method]
  h[is.na(best_y_lag), best_y_lag := best_y[1], by = method]
  h[, improvement := best_y - best_y_lag]
  h[]
}

runs <- list(EI = res_EI, UCB = res_UCB)

hist_dt <- rbindlist(lapply(names(runs), function(m) {
  make_hist_dt(runs[[m]], m)
}))

setorder(hist_dt, method, t)

convergence_summary <- hist_dt[, .(
  T = max(t),
  init_best = best_y[t == min(t)],
  final_best = best_y[t == max(t)],
  total_improvement = best_y[t == max(t)] - best_y[t == min(t)],
  n_improvements = sum(improvement > 0, na.rm = TRUE)
), by = method]

best_cfg_dt <- rbindlist(lapply(names(runs), function(m) {
  data.table(
    method = m,
    best_y = runs[[m]]$best_y,
    best_id = runs[[m]]$best_id,
    cfg = cfg_key(runs[[m]]$best_cfg)
  )
}))

cat("\n=== Convergence summary ===\n")
print(convergence_summary)

cat("\n=== Best configuration by method ===\n")
print(best_cfg_dt)

p_best <- ggplot(hist_dt, aes(x = t, y = best_y, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = oracle_best_val, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Iteration", y = "Best-so-far objective") +
  safe_theme() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    plot.title = element_blank()
  )

print(p_best)

p_selected <- ggplot(hist_dt, aes(x = t, y = y, color = method)) +
  geom_line(alpha = 0.8) +
  geom_point(size = 1.2, alpha = 0.8) +
  labs(title = "Selected Objective Values", x = "Iteration", y = "Selected objective") +
  safe_theme()

print(p_selected)

hist_dt[, y_roll := frollmean(y, n = 10, align = "right"), by = method]
hist_dt[, regret := oracle_best_val - best_y]
hist_dt[, regret_eps := pmax(regret, 1e-8)]

p_rolling <- ggplot(hist_dt, aes(x = t, y = y_roll, color = method)) +
  geom_line(linewidth = 1) +
  labs(title = "Selected Objective Rolling Mean", x = "Iteration", y = "Rolling mean") +
  safe_theme()

print(p_rolling)

p_regret <- ggplot(hist_dt, aes(x = t, y = regret, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2, alpha = 0.8) +
  labs(title = "Simple Regret", x = "Iteration", y = "Oracle - best-so-far") +
  safe_theme()

print(p_regret)

p_log_regret <- ggplot(hist_dt, aes(x = t, y = log10(regret_eps), color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2, alpha = 0.8) +
  labs(title = "Log10 Regret", x = "Iteration", y = "log10(regret)") +
  safe_theme()

print(p_log_regret)

plotly_convergence <- plot_ly(hist_dt, x = ~t, color = ~method) %>%
  add_lines(y = ~best_y, name = "Best-so-far") %>%
  add_markers(y = ~best_y, name = "Best-so-far points", marker = list(size = 5)) %>%
  layout(
    title = "Convergence",
    xaxis = list(title = "Iteration"),
    yaxis = list(title = "Best-so-far objective")
  )

print(plotly_convergence)

p_histogram <- ggplot(hist_dt, aes(x = y, fill = method)) +
  geom_histogram(bins = 25, alpha = 0.65, position = "identity") +
  labs(title = "Selected Objective Distribution", x = "Selected objective", y = "Count") +
  safe_theme()

print(p_histogram)

chosen_cfgs <- rbindlist(lapply(names(runs), function(m) {
  ids <- runs[[m]]$history$chosen_id
  dt <- as.data.table(X_all[ids, ])
  dt[, method := m]
  dt
}))

for (j in seq_len(n_variables)) {
  chosen_cfgs[[paste0("c", j)]] <- as.integer(chosen_cfgs[[paste0("c", j)]])
}

level_counts <- rbindlist(lapply(seq_len(n_variables), function(j) {
  chosen_cfgs[, .N, by = .(method, level = get(paste0("c", j)))][
    , variable := paste0("c", j)
  ]
}))

p_level_counts <- ggplot(level_counts, aes(x = factor(level), y = N, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~variable, ncol = 3) +
  labs(title = "Categorical Level Selection Frequency", x = "Level", y = "Count") +
  safe_theme()

print(p_level_counts)

land_dt <- as.data.table(X_all)

for (j in seq_len(n_variables)) {
  land_dt[[paste0("c", j)]] <- as.integer(land_dt[[paste0("c", j)]])
}

land_dt[, y_true := y_true]
land_dt[, cfg_id := .I]
land_dt[, cfg := apply(.SD, 1, paste, collapse = "-"), .SDcols = variable_names]

slice_dt <- copy(land_dt)[c4 == 1 & c5 == 1]

p_slice <- ggplot(slice_dt, aes(x = factor(c1), y = factor(c2), fill = y_true)) +
  geom_tile() +
  facet_wrap(~c3, ncol = 4) +
  scale_fill_viridis_c() +
  labs(
    title = "Landscape Slices with c4 = 1 and c5 = 1",
    x = "c1 level",
    y = "c2 level",
    fill = "Objective"
  ) +
  safe_theme()

print(p_slice)

violin_dt <- rbindlist(lapply(seq_len(n_variables), function(j) {
  dt <- land_dt[, .(level = get(paste0("c", j)), y_true)]
  dt[, variable := paste0("c", j)]
  dt
}))

p_violin <- ggplot(violin_dt, aes(x = factor(level), y = y_true)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.2) +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Objective Distribution by Level", x = "Level", y = "Objective") +
  safe_theme()

print(p_violin)

hamming_dist <- function(A) {
  n <- nrow(A)
  D <- matrix(0, n, n)
  
  for (i in seq_len(n)) {
    D[i, ] <- rowSums(t(t(A) != A[i, ]))
  }
  
  as.dist(D)
}

A <- as.matrix(land_dt[, ..variable_names])
Dham <- hamming_dist(A)

hc <- hclust(Dham, method = "average")
plot(hc, labels = FALSE, main = "Hierarchical Clustering with Hamming Distance")

mds_xy <- cmdscale(Dham, k = 2)

mds_dt <- data.table(
  x = mds_xy[, 1],
  y = mds_xy[, 2],
  y_true = land_dt$y_true,
  cfg_id = land_dt$cfg_id
)

p_mds <- ggplot(mds_dt, aes(x = x, y = y, color = y_true)) +
  geom_point(size = 1.3, alpha = 0.85) +
  scale_color_viridis_c() +
  labs(
    title = "MDS Embedding of the Categorical Space",
    x = "MDS 1",
    y = "MDS 2",
    color = "Objective"
  ) +
  safe_theme()

print(p_mds)

set.seed(42)

idx_tsne <- sample(seq_len(nrow(A)), 400)
Dsub <- as.matrix(Dham)[idx_tsne, idx_tsne]
tsne <- Rtsne(Dsub, is_distance = TRUE, perplexity = 15)

tsne_dt <- data.table(
  x = tsne$Y[, 1],
  y = tsne$Y[, 2],
  y_true = land_dt$y_true[idx_tsne]
)

p_tsne <- ggplot(tsne_dt, aes(x = x, y = y, color = y_true)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_viridis_c() +
  labs(
    title = "t-SNE Embedding with Hamming Distance",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Objective"
  ) +
  safe_theme()

print(p_tsne)

make_choice_dataset_from_result <- function(result) {
  eval_ids <- result$history$chosen_id
  y_obs <- result$history$y
  names(y_obs) <- as.character(eval_ids)
  
  choice_long <- NULL
  
  for (t in seq.int(2, length(eval_ids))) {
    ids_t <- eval_ids[seq_len(t)]
    y_t <- y_obs[as.character(ids_t)]
    cl <- make_choice_long(ids_t, y_t, chid = t - 1)
    
    choice_long <- if (is.null(choice_long)) {
      cl
    } else {
      rbind(choice_long, cl)
    }
  }
  
  choice_long
}

build_surrogate_diag <- function(result) {
  choice_long <- make_choice_dataset_from_result(result)
  model <- fit_mnl_surrogate(choice_long)
  pred <- predict_mnl_all(model)
  
  data.table(
    cfg_id = seq_along(pred$mu),
    mu = pred$mu,
    sigma = pred$sigma,
    pi = pred$pi,
    y_true = y_true,
    evaluated = seq_along(pred$mu) %in% result$history$chosen_id
  )
}

edges <- rbindlist(lapply(interaction_pairs, function(pk) {
  nm <- paste0("c", pk[1], "_c", pk[2])
  data.table(
    from = paste0("c", pk[1]),
    to = paste0("c", pk[2]),
    weight = mean(abs(theta$theta_pair[[nm]]))
  )
}))

g <- graph_from_data_frame(edges, directed = FALSE)

plot(
  g,
  vertex.size = 28,
  vertex.label.cex = 1.1,
  edge.width = 2 + 10 * (E(g)$weight / max(E(g)$weight)),
  main = "Variable Interaction Graph"
)

vis_nodes <- data.frame(id = variable_names, label = variable_names)

vis_edges <- data.frame(
  from = edges$from,
  to = edges$to,
  value = edges$weight,
  title = paste0("mean absolute interaction = ", signif(edges$weight, 3))
)

visNetwork(vis_nodes, vis_edges, main = "Interaction Graph") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

main_dt <- rbindlist(lapply(seq_len(n_variables), function(j) {
  data.table(
    variable = paste0("c", j),
    level = levels4,
    theta = theta$theta_main[[j]]
  )
}))

best_level <- main_dt[, .(best_level = level[which.max(theta)]), by = variable]

cat("\n=== Best main-effect level by variable ===\n")
print(best_level)

set.seed(999)

seeds <- seq_len(30)
T_budget <- 60
n0 <- 10
noise_sd <- 0.05
kappa <- 2.0
methods <- c("EI", "UCB")

all_runs <- list()

for (s in seeds) {
  set.seed(s)
  
  for (m in methods) {
    result <- run_mnl_bo(
      T_budget = T_budget,
      n0 = n0,
      noise_sd = noise_sd,
      kappa = kappa,
      acq = m
    )
    
    dt <- as.data.table(result$history)
    dt[, seed := s]
    dt[, method := m]
    
    all_runs[[length(all_runs) + 1]] <- dt
  }
}

bo_dt <- rbindlist(all_runs)
setorder(bo_dt, method, seed, t)

summary_table <- bo_dt[, .(
  mean_init_best = mean(best_y[t == min(t)]),
  sd_init_best = sd(best_y[t == min(t)]),
  mean_final_best = mean(best_y[t == max(t)]),
  sd_final_best = sd(best_y[t == max(t)]),
  mean_improvement = mean(best_y[t == max(t)] - best_y[t == min(t)])
), by = method]

cat("\n=== Multi-seed final performance summary ===\n")
print(summary_table)

curve_dt <- bo_dt[, .(
  mean_best = mean(best_y),
  sd_best = sd(best_y)
), by = .(method, t)]

curve_dt[, upper := mean_best + sd_best]
curve_dt[, lower := mean_best - sd_best]

p_mean_curve <- ggplot(curve_dt, aes(x = t, y = mean_best, color = method)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = method),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Iteration", y = "Best-so-far objective") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    plot.title = element_blank()
  )

print(p_mean_curve)

make_pairwise_long <- function(eval_ids, y_obs, n_pairs = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  eval_ids <- as.integer(eval_ids)
  if (length(eval_ids) < 2) return(NULL)
  
  sampled_pairs <- t(replicate(n_pairs, sample(eval_ids, 2, replace = FALSE)))
  dt_list <- vector("list", nrow(sampled_pairs))
  
  for (k in seq_len(nrow(sampled_pairs))) {
    i <- sampled_pairs[k, 1]
    j <- sampled_pairs[k, 2]
    
    yi <- y_obs[as.character(i)]
    yj <- y_obs[as.character(j)]
    
    if (is.na(yi) || is.na(yj) || yi == yj) next
    
    winner <- if (yi > yj) i else j
    loser <- if (yi > yj) j else i
    
    dt_list[[k]] <- data.table(
      chid = k,
      alt = as.character(c(winner, loser)),
      choice = c(1L, 0L)
    )
  }
  
  pair_dt <- rbindlist(dt_list, fill = TRUE)
  pair_dt[!is.na(chid)]
}

fit_mnl_pairs <- function(pair_long_dt) {
  if (is.null(pair_long_dt) || nrow(pair_long_dt) < 10) return(NULL)
  
  alt_id <- as.integer(pair_long_dt$alt)
  pair_df <- cbind(as.data.frame(pair_long_dt), X_mm[alt_id, , drop = FALSE])
  
  mdat <- mlogit.data(
    pair_df,
    choice = "choice",
    shape = "long",
    chid.var = "chid",
    alt.var = "alt"
  )
  
  formula_mnl <- as.formula(
    paste0("choice ~ 0 + ", paste(colnames(X_mm), collapse = " + "))
  )
  
  tryCatch(mlogit(formula_mnl, data = mdat), error = function(e) NULL)
}

prob_mass_topK_over_time <- function(result, K = 10, n_pairs = 400) {
  y_obs <- result$history$y
  eval_ids <- result$history$chosen_id
  names(y_obs) <- as.character(eval_ids)
  
  top_ids <- order(-y_true)[seq_len(K)]
  out <- vector("list", nrow(result$history) - 2)
  idx <- 1
  
  for (t in seq.int(3, nrow(result$history))) {
    ids_t <- eval_ids[seq_len(t)]
    y_t <- y_obs[as.character(ids_t)]
    
    pair_long <- make_pairwise_long(ids_t, y_t, n_pairs = n_pairs, seed = 100 + t)
    model <- fit_mnl_pairs(pair_long)
    pred <- predict_mnl_all(model)
    
    out[[idx]] <- data.table(
      iter = t,
      prob_topK = sum(pred$pi[top_ids]),
      fit_ok = !is.null(model)
    )
    
    idx <- idx + 1
  }
  
  rbindlist(out, fill = TRUE)
}

fit_final_pairwise_diag <- function(result, n_pairs = 2000) {
  eval_ids <- result$history$chosen_id
  y_obs <- result$history$y
  names(y_obs) <- as.character(eval_ids)
  
  pair_long <- make_pairwise_long(eval_ids, y_obs, n_pairs = n_pairs, seed = 777)
  model <- fit_mnl_pairs(pair_long)
  pred <- predict_mnl_all(model)
  
  data.table(
    cfg_id = seq_along(pred$mu),
    mu = pred$mu,
    pi = pred$pi,
    sigma = pred$sigma,
    y_true = y_true,
    evaluated = seq_along(pred$mu) %in% eval_ids
  )
}

diag_EI <- fit_final_pairwise_diag(res_EI)

diag_EI[, mu_rank := frank(mu, ties.method = "average") / .N]
diag_EI[, mu_bin := cut(mu_rank, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)]

calib_dt <- diag_EI[, .(
  mean_true = mean(y_true),
  sd_true = sd(y_true),
  n = .N
), by = mu_bin]

p_calibration <- ggplot(calib_dt, aes(x = mu_bin, y = mean_true)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_true - sd_true, ymax = mean_true + sd_true),
    width = 0.2
  ) +
  labs(
    title = "Surrogate Calibration with Pairwise-Preference MNL",
    x = "Predicted utility quantile bin",
    y = "Mean noise-free objective"
  ) +
  safe_theme() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_calibration)

cat("\n=== Completed MNL-BO categorical benchmark analysis ===\n")
