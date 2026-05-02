############################################################
# MNL-based Bayesian Optimization for TSP
#
# Problem:
#   Single-objective TSP minimization with a fixed distance matrix
#
# Representation:
#   Continuous keys in [0, 1]^n define a tour by sorting cities
#
# Surrogates:
#   MNL-BO: conditional logit trained from preference choice sets
#   RF-SMAC: random forest surrogate with tree-wise uncertainty
#
# Baselines:
#   Random Search
#   2-opt Local Search
#   Simulated Annealing
#   Restarted Hill Climbing
#
# Outputs:
#   Single-run best route table
#   Single-run progress table
#   Multi-run summary table
#   Single-run convergence plot
#   Best route MDS plot
#   Multi-run convergence plot
############################################################

required_packages <- c("data.table", "mlogit", "ggplot2", "ranger")

install_missing_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  
  if (length(missing_packages) > 0) {
    install.packages(missing_packages)
  }
}

install_missing_packages(required_packages)

suppressPackageStartupMessages({
  library(data.table)
  library(mlogit)
  library(ggplot2)
  library(ranger)
})

set.seed(123)

############################################################
# 1. Cities and distance matrix
############################################################

cities <- c(
  "Atlanta", "Boston", "Chicago", "Dallas", "Denver", "Houston",
  "Las Vegas", "Los Angeles", "Miami", "New Orleans", "New York",
  "Phoenix", "San Francisco", "Seattle", "Washington D.C."
)

distance_matrix <- matrix(
  c(
    0, 1095, 715, 805, 1437, 844, 1920, 2230, 675, 499, 884, 1832, 2537, 2730, 657,
    1095, 0, 983, 1815, 1991, 1886, 2500, 3036, 1539, 1541, 213, 2664, 3179, 3043, 44,
    715, 983, 0, 931, 1050, 1092, 1500, 2112, 1390, 947, 840, 1729, 2212, 2052, 695,
    805, 1815, 931, 0, 801, 242, 1150, 1425, 1332, 504, 1604, 1027, 1765, 2122, 1372,
    1437, 1991, 1050, 801, 0, 1032, 885, 1174, 2094, 1305, 1780, 836, 1266, 1373, 1635,
    844, 1886, 1092, 242, 1032, 0, 1525, 1556, 1237, 365, 1675, 1158, 1958, 2348, 1443,
    1920, 2500, 1500, 1150, 885, 1525, 0, 289, 2640, 1805, 2486, 294, 573, 1188, 2568,
    2230, 3036, 2112, 1425, 1174, 1556, 289, 0, 2757, 1921, 2825, 398, 403, 1150, 2680,
    675, 1539, 1390, 1332, 2094, 1237, 2640, 2757, 0, 892, 1328, 2359, 3097, 3389, 1101,
    499, 1541, 947, 504, 1305, 365, 1805, 1921, 892, 0, 1330, 1523, 2269, 2626, 1098,
    884, 213, 840, 1604, 1780, 1675, 2486, 2825, 1328, 1330, 0, 2442, 3036, 2900, 229,
    1832, 2664, 1729, 1027, 836, 1158, 294, 398, 2359, 1523, 2442, 0, 800, 1482, 2278,
    2537, 3179, 2212, 1765, 1266, 1958, 573, 403, 3097, 2269, 3036, 800, 0, 817, 2864,
    2730, 3043, 2052, 2122, 1373, 2348, 1188, 1150, 3389, 2626, 2900, 1482, 817, 0, 2755,
    657, 44, 695, 1372, 1635, 1443, 2568, 2680, 1101, 1098, 229, 2278, 2864, 2755, 0
  ),
  nrow = length(cities),
  byrow = TRUE
)

dimnames(distance_matrix) <- list(cities, cities)
n_cities <- length(cities)

############################################################
# 2. Core TSP utilities
############################################################

calculate_distance <- function(path) {
  total <- 0
  
  for (i in seq_len(length(path) - 1)) {
    total <- total + distance_matrix[path[i], path[i + 1]]
  }
  
  total + distance_matrix[path[length(path)], path[1]]
}

keys_to_path <- function(keys) {
  cities[order(keys)]
}

tour_key <- function(path) {
  paste(path, collapse = "|")
}

tour_string <- function(path) {
  paste(c(path, path[1]), collapse = " \u2192 ")
}

############################################################
# 3. Tour features
############################################################

cities_safe <- make.names(cities)
city_map <- setNames(cities_safe, cities)

make_features <- function(path) {
  feature_names <- as.vector(
    outer(paste0("pos", seq_len(n_cities), "_"), cities_safe, paste0)
  )
  
  z <- setNames(rep(0, length(feature_names)), feature_names)
  
  for (pos in seq_len(n_cities)) {
    feature_name <- paste0("pos", pos, "_", city_map[[path[pos]]])
    z[[feature_name]] <- 1
  }
  
  z
}

tsp_features_df <- function(paths) {
  template <- make_features(paths[[1]])
  
  Z <- t(vapply(
    paths,
    function(path) make_features(path),
    FUN.VALUE = template
  ))
  
  as.data.frame(Z, check.names = FALSE)
}

############################################################
# 4. Preference data for MNL
############################################################

build_choice_set <- function(N, K) {
  sample.int(N, size = min(K, N), replace = FALSE)
}

make_mlogit_data <- function(choice_sets, paths, dists) {
  rows <- vector("list", length = 0)
  set_id <- 1
  
  for (idx in choice_sets) {
    winner <- idx[which.min(dists[idx])]
    
    for (j in idx) {
      z <- make_features(paths[[j]])
      
      rows[[length(rows) + 1]] <- data.frame(
        set_id = set_id,
        alt_id = paste0("a", j),
        choice = as.integer(j == winner),
        t(z),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
    
    set_id <- set_id + 1
  }
  
  df <- do.call(rbind, rows)
  
  mlogit.data(
    df,
    choice = "choice",
    shape = "long",
    chid.var = "set_id",
    alt.var = "alt_id"
  )
}

fit_mnl_surrogate <- function(mdata) {
  cols <- setdiff(colnames(mdata), c("choice", "idx"))
  
  formula_mnl <- as.formula(
    paste("choice ~ 0 +", paste(cols, collapse = " + "))
  )
  
  tryCatch(
    mlogit(formula_mnl, data = mdata),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(mlogit(formula_mnl, data = mdata))
  )
}

predict_mu_sigma <- function(model, paths_new) {
  beta <- coef(model)
  covariance <- vcov(model)
  
  template <- make_features(paths_new[[1]])
  
  Z <- t(vapply(
    paths_new,
    function(path) make_features(path),
    FUN.VALUE = template
  ))
  
  Zdf <- as.data.frame(Z, check.names = FALSE)
  
  for (nm in setdiff(names(beta), names(Zdf))) {
    Zdf[[nm]] <- 0
  }
  
  Zmat <- as.matrix(Zdf[, names(beta), drop = FALSE])
  
  mu <- as.numeric(Zmat %*% beta)
  sigma2 <- rowSums((Zmat %*% covariance) * Zmat)
  sigma <- sqrt(pmax(sigma2, 1e-12))
  
  list(mu = mu, sigma = sigma)
}

############################################################
# 5. Acquisition functions
############################################################

ei_max <- function(mu, sigma, best) {
  sigma <- pmax(sigma, 1e-12)
  z <- (mu - best) / sigma
  (mu - best) * pnorm(z) + sigma * dnorm(z)
}

ucb_max <- function(mu, sigma, kappa = 2) {
  mu + kappa * sigma
}

############################################################
# 6. Local search moves and candidate generation
############################################################

swap_move <- function(path) {
  n <- length(path)
  i <- sample.int(n, 1)
  j <- sample(setdiff(seq_len(n), i), 1)
  
  out <- path
  temp <- out[i]
  out[i] <- out[j]
  out[j] <- temp
  
  out
}

insertion_move <- function(path) {
  n <- length(path)
  i <- sample.int(n, 1)
  j <- sample(setdiff(seq_len(n), i), 1)
  
  out <- path
  city <- out[i]
  out <- out[-i]
  append(out, city, after = j - (j > i))
}

two_opt_move <- function(path) {
  n <- length(path)
  i <- sample.int(n - 1, 1)
  j <- sample(seq.int(i + 1, n), 1)
  
  out <- path
  out[i:j] <- rev(out[i:j])
  
  out
}

random_neighbor <- function(path, ops = c("2opt", "swap", "insert")) {
  op <- sample(ops, 1)
  
  if (op == "2opt") return(two_opt_move(path))
  if (op == "swap") return(swap_move(path))
  if (op == "insert") return(insertion_move(path))
  
  stop("Unknown move type.")
}

generate_candidates <- function(
    best_path,
    n_local = 90,
    n_global = 30,
    op = c("2opt", "swap", "insert"),
    cities
) {
  op <- match.arg(op)
  
  local_paths <- vector("list", n_local)
  
  for (k in seq_len(n_local)) {
    if (op == "2opt") local_paths[[k]] <- two_opt_move(best_path)
    if (op == "swap") local_paths[[k]] <- swap_move(best_path)
    if (op == "insert") local_paths[[k]] <- insertion_move(best_path)
  }
  
  global_paths <- replicate(
    n_global,
    sample(cities, length(cities), replace = FALSE),
    simplify = FALSE
  )
  
  all_paths <- c(local_paths, global_paths)
  keys <- vapply(all_paths, tour_key, character(1))
  
  all_paths[!duplicated(keys)]
}

############################################################
# 7. MNL-BO for TSP
############################################################

run_mnl_bo_tsp <- function(
    n_init = 30,
    n_iter = 150,
    n_choice_sets = 40,
    K_choice = 5,
    n_cand = 120,
    kappa = 2,
    acq_mode = c("UCB", "EI"),
    op = c("2opt", "swap", "insert"),
    seed = 1
) {
  set.seed(seed)
  
  acq_mode <- match.arg(acq_mode)
  op <- match.arg(op)
  
  keys_eval <- replicate(n_init, runif(n_cities), simplify = FALSE)
  paths_eval <- lapply(keys_eval, keys_to_path)
  dists_eval <- vapply(paths_eval, calculate_distance, numeric(1))
  utils_eval <- -dists_eval
  
  best_idx <- which.min(dists_eval)
  best_path <- paths_eval[[best_idx]]
  best_dist <- dists_eval[best_idx]
  best_utility <- max(utils_eval)
  
  N_total <- n_init + n_iter
  
  log <- data.table(
    step = seq_len(N_total),
    best_dist = NA_real_,
    selected_dist = NA_real_
  )
  
  log$selected_dist[seq_len(n_init)] <- dists_eval
  log$best_dist[seq_len(n_init)] <- cummin(dists_eval)
  
  seen <- new.env(parent = emptyenv())
  
  for (i in seq_along(paths_eval)) {
    seen[[tour_key(paths_eval[[i]])]] <- dists_eval[i]
  }
  
  for (iter in seq_len(n_iter)) {
    choice_sets <- replicate(
      n_choice_sets,
      build_choice_set(length(paths_eval), K_choice),
      simplify = FALSE
    )
    
    mdata <- make_mlogit_data(choice_sets, paths_eval, dists_eval)
    model <- fit_mnl_surrogate(mdata)
    
    cand_paths <- generate_candidates(
      best_path = best_path,
      n_local = 90,
      n_global = 30,
      op = op,
      cities = cities
    )
    
    if (length(cand_paths) > n_cand) {
      cand_paths <- cand_paths[seq_len(n_cand)]
    }
    
    if (length(cand_paths) < 2) {
      cand_paths <- replicate(
        10,
        sample(cities, n_cities, replace = FALSE),
        simplify = FALSE
      )
    }
    
    if (is.null(model)) {
      pick <- sample.int(length(cand_paths), 1)
    } else {
      pred <- predict_mu_sigma(model, cand_paths)
      
      acquisition <- if (acq_mode == "UCB") {
        ucb_max(pred$mu, pred$sigma, kappa = kappa)
      } else {
        ei_max(pred$mu, pred$sigma, best = best_utility)
      }
      
      pick <- which.max(acquisition)
    }
    
    next_path <- cand_paths[[pick]]
    key <- tour_key(next_path)
    
    if (!is.null(seen[[key]])) {
      next_dist <- seen[[key]]
    } else {
      next_dist <- calculate_distance(next_path)
      seen[[key]] <- next_dist
    }
    
    next_utility <- -next_dist
    
    paths_eval <- c(paths_eval, list(next_path))
    dists_eval <- c(dists_eval, next_dist)
    utils_eval <- c(utils_eval, next_utility)
    
    if (next_dist < best_dist) {
      best_dist <- next_dist
      best_path <- next_path
    }
    
    best_utility <- max(best_utility, next_utility)
    
    idx <- n_init + iter
    log$selected_dist[idx] <- next_dist
    log$best_dist[idx] <- best_dist
  }
  
  list(
    best_dist = best_dist,
    best_path = best_path,
    log = log,
    paths_eval = paths_eval,
    dists_eval = dists_eval
  )
}

############################################################
# 8. Baseline methods
############################################################

run_random_search_tsp <- function(n_init = 30, n_iter = 150, seed = 1) {
  set.seed(seed)
  
  N_total <- n_init + n_iter
  best_dist <- Inf
  best_path <- NULL
  
  log <- data.table(
    step = seq_len(N_total),
    best_dist = NA_real_,
    selected_dist = NA_real_
  )
  
  seen <- new.env(parent = emptyenv())
  
  for (t in seq_len(N_total)) {
    path <- sample(cities, n_cities, replace = FALSE)
    key <- tour_key(path)
    
    if (!is.null(seen[[key]])) {
      dist <- seen[[key]]
    } else {
      dist <- calculate_distance(path)
      seen[[key]] <- dist
    }
    
    if (dist < best_dist) {
      best_dist <- dist
      best_path <- path
    }
    
    log$selected_dist[t] <- dist
    log$best_dist[t] <- best_dist
  }
  
  list(best_dist = best_dist, best_path = best_path, log = log)
}

run_twoopt_local_search_tsp <- function(n_init = 30, n_iter = 150, seed = 1) {
  set.seed(seed)
  
  N_total <- n_init + n_iter
  
  init_paths <- replicate(
    n_init,
    sample(cities, n_cities, replace = FALSE),
    simplify = FALSE
  )
  
  init_dists <- vapply(init_paths, calculate_distance, numeric(1))
  
  best_idx <- which.min(init_dists)
  best_path <- init_paths[[best_idx]]
  best_dist <- init_dists[best_idx]
  
  log <- data.table(
    step = seq_len(N_total),
    best_dist = NA_real_,
    selected_dist = NA_real_
  )
  
  log$selected_dist[seq_len(n_init)] <- init_dists
  log$best_dist[seq_len(n_init)] <- cummin(init_dists)
  
  seen <- new.env(parent = emptyenv())
  
  for (i in seq_along(init_paths)) {
    seen[[tour_key(init_paths[[i]])]] <- init_dists[i]
  }
  
  for (t in seq.int(n_init + 1, N_total)) {
    candidate <- two_opt_move(best_path)
    key <- tour_key(candidate)
    
    if (!is.null(seen[[key]])) {
      dist <- seen[[key]]
    } else {
      dist <- calculate_distance(candidate)
      seen[[key]] <- dist
    }
    
    if (dist < best_dist) {
      best_dist <- dist
      best_path <- candidate
    }
    
    log$selected_dist[t] <- dist
    log$best_dist[t] <- best_dist
  }
  
  list(best_dist = best_dist, best_path = best_path, log = log)
}

run_simulated_annealing_tsp <- function(
    n_init = 30,
    n_iter = 150,
    seed = 1,
    temp0 = 2000,
    cooling = 0.98
) {
  set.seed(seed)
  
  N_total <- n_init + n_iter
  
  init_paths <- replicate(
    n_init,
    sample(cities, n_cities, replace = FALSE),
    simplify = FALSE
  )
  
  init_dists <- vapply(init_paths, calculate_distance, numeric(1))
  
  best_idx <- which.min(init_dists)
  current_path <- init_paths[[best_idx]]
  current_dist <- init_dists[best_idx]
  best_path <- current_path
  best_dist <- current_dist
  
  log <- data.table(
    step = seq_len(N_total),
    best_dist = NA_real_,
    selected_dist = NA_real_
  )
  
  log$selected_dist[seq_len(n_init)] <- init_dists
  log$best_dist[seq_len(n_init)] <- cummin(init_dists)
  
  seen <- new.env(parent = emptyenv())
  
  for (i in seq_along(init_paths)) {
    seen[[tour_key(init_paths[[i]])]] <- init_dists[i]
  }
  
  temperature <- temp0
  
  for (iter in seq_len(n_iter)) {
    candidate <- random_neighbor(current_path)
    key <- tour_key(candidate)
    
    if (!is.null(seen[[key]])) {
      next_dist <- seen[[key]]
    } else {
      next_dist <- calculate_distance(candidate)
      seen[[key]] <- next_dist
    }
    
    delta <- next_dist - current_dist
    
    if (delta < 0 || runif(1) < exp(-delta / temperature)) {
      current_path <- candidate
      current_dist <- next_dist
    }
    
    if (current_dist < best_dist) {
      best_dist <- current_dist
      best_path <- current_path
    }
    
    idx <- n_init + iter
    log$selected_dist[idx] <- next_dist
    log$best_dist[idx] <- best_dist
    
    temperature <- temperature * cooling
  }
  
  list(best_dist = best_dist, best_path = best_path, log = log)
}

run_restarted_hillclimb_tsp <- function(
    n_init = 30,
    n_iter = 150,
    seed = 1,
    restart_patience = 15
) {
  set.seed(seed)
  
  N_total <- n_init + n_iter
  
  init_paths <- replicate(
    n_init,
    sample(cities, n_cities, replace = FALSE),
    simplify = FALSE
  )
  
  init_dists <- vapply(init_paths, calculate_distance, numeric(1))
  
  best_idx <- which.min(init_dists)
  current_path <- init_paths[[best_idx]]
  current_dist <- init_dists[best_idx]
  best_path <- current_path
  best_dist <- current_dist
  
  log <- data.table(
    step = seq_len(N_total),
    best_dist = NA_real_,
    selected_dist = NA_real_
  )
  
  log$selected_dist[seq_len(n_init)] <- init_dists
  log$best_dist[seq_len(n_init)] <- cummin(init_dists)
  
  seen <- new.env(parent = emptyenv())
  
  for (i in seq_along(init_paths)) {
    seen[[tour_key(init_paths[[i]])]] <- init_dists[i]
  }
  
  no_improve <- 0
  
  for (iter in seq_len(n_iter)) {
    candidate <- random_neighbor(current_path)
    key <- tour_key(candidate)
    
    if (!is.null(seen[[key]])) {
      next_dist <- seen[[key]]
    } else {
      next_dist <- calculate_distance(candidate)
      seen[[key]] <- next_dist
    }
    
    if (next_dist < current_dist) {
      current_path <- candidate
      current_dist <- next_dist
      no_improve <- 0
    } else {
      no_improve <- no_improve + 1
    }
    
    if (no_improve >= restart_patience) {
      restart_path <- sample(cities, n_cities, replace = FALSE)
      restart_key <- tour_key(restart_path)
      
      if (!is.null(seen[[restart_key]])) {
        restart_dist <- seen[[restart_key]]
      } else {
        restart_dist <- calculate_distance(restart_path)
        seen[[restart_key]] <- restart_dist
      }
      
      current_path <- restart_path
      current_dist <- restart_dist
      next_dist <- restart_dist
      no_improve <- 0
    }
    
    if (current_dist < best_dist) {
      best_dist <- current_dist
      best_path <- current_path
    }
    
    idx <- n_init + iter
    log$selected_dist[idx] <- next_dist
    log$best_dist[idx] <- best_dist
  }
  
  list(best_dist = best_dist, best_path = best_path, log = log)
}

############################################################
# 9. RF-SMAC baseline
############################################################

rf_predict_mean_sd <- function(rf_model, X_new) {
  predictions <- predict(rf_model, data = X_new, predict.all = TRUE)$predictions
  
  list(
    mu = rowMeans(predictions),
    sigma = apply(predictions, 1, sd)
  )
}

run_rf_smac_tsp <- function(
    n_init = 30,
    n_iter = 150,
    n_cand = 120,
    kappa = 2,
    acq_mode = c("EI", "UCB"),
    op = c("2opt", "swap", "insert"),
    seed = 1,
    num_trees = 500,
    min_node = 5
) {
  set.seed(seed)
  
  acq_mode <- match.arg(acq_mode)
  op <- match.arg(op)
  
  paths_eval <- replicate(
    n_init,
    sample(cities, n_cities, replace = FALSE),
    simplify = FALSE
  )
  
  dists_eval <- vapply(paths_eval, calculate_distance, numeric(1))
  utils_eval <- -dists_eval
  
  best_idx <- which.max(utils_eval)
  best_path <- paths_eval[[best_idx]]
  best_utility <- utils_eval[best_idx]
  
  N_total <- n_init + n_iter
  
  log <- data.table(
    step = seq_len(N_total),
    best_dist = NA_real_,
    selected_dist = NA_real_
  )
  
  log$selected_dist[seq_len(n_init)] <- dists_eval
  log$best_dist[seq_len(n_init)] <- cummin(dists_eval)
  
  seen <- new.env(parent = emptyenv())
  
  for (i in seq_along(paths_eval)) {
    seen[[tour_key(paths_eval[[i]])]] <- dists_eval[i]
  }
  
  for (iter in seq_len(n_iter)) {
    X_train <- tsp_features_df(paths_eval)
    
    rf_model <- ranger(
      x = X_train,
      y = utils_eval,
      num.trees = num_trees,
      min.node.size = min_node,
      seed = seed + iter
    )
    
    cand_paths <- generate_candidates(
      best_path = best_path,
      n_local = 90,
      n_global = 30,
      op = op,
      cities = cities
    )
    
    if (length(cand_paths) > n_cand) {
      cand_paths <- cand_paths[seq_len(n_cand)]
    }
    
    if (length(cand_paths) < 2) {
      cand_paths <- replicate(
        10,
        sample(cities, n_cities, replace = FALSE),
        simplify = FALSE
      )
    }
    
    X_cand <- tsp_features_df(cand_paths)
    pred <- rf_predict_mean_sd(rf_model, X_cand)
    
    acquisition <- if (acq_mode == "EI") {
      ei_max(pred$mu, pred$sigma, best = best_utility)
    } else {
      ucb_max(pred$mu, pred$sigma, kappa = kappa)
    }
    
    pick <- which.max(acquisition)
    next_path <- cand_paths[[pick]]
    key <- tour_key(next_path)
    
    if (!is.null(seen[[key]])) {
      next_dist <- seen[[key]]
    } else {
      next_dist <- calculate_distance(next_path)
      seen[[key]] <- next_dist
    }
    
    next_utility <- -next_dist
    
    paths_eval <- c(paths_eval, list(next_path))
    dists_eval <- c(dists_eval, next_dist)
    utils_eval <- c(utils_eval, next_utility)
    
    if (next_utility > best_utility) {
      best_utility <- next_utility
      best_path <- next_path
    }
    
    idx <- n_init + iter
    log$selected_dist[idx] <- next_dist
    log$best_dist[idx] <- min(log$best_dist[idx - 1], next_dist)
  }
  
  list(
    best_dist = min(dists_eval),
    best_path = best_path,
    log = log
  )
}

############################################################
# 10. Single-run study
############################################################

single_run <- run_mnl_bo_tsp(
  n_init = 30,
  n_iter = 150,
  n_choice_sets = 40,
  K_choice = 5,
  n_cand = 120,
  kappa = 2,
  acq_mode = "UCB",
  op = "2opt",
  seed = 123
)

bo_log <- copy(single_run$log)
best_path <- single_run$best_path
best_dist <- single_run$best_dist
paths_eval <- single_run$paths_eval
dists_eval <- single_run$dists_eval

bo_log[, improvement := shift(best_dist, type = "lag") - best_dist]
bo_log[is.na(improvement), improvement := 0]

tab_best_route <- data.table(
  best_distance = best_dist,
  tour = tour_string(best_path)
)

top_tours <- data.table(
  tour_id = seq_along(paths_eval),
  tour_key = vapply(paths_eval, tour_key, character(1)),
  distance = dists_eval,
  tour = vapply(paths_eval, tour_string, character(1))
)

setorder(top_tours, distance)
top_tours <- top_tours[!duplicated(tour_key)]

progress_checkpoints <- unique(c(1, 5, 10, 25, 50, 75, 100, 125, 150))
progress_checkpoints <- progress_checkpoints[
  progress_checkpoints <= max(bo_log$step)
]

best_tour_up_to <- function(step) {
  n_init <- length(dists_eval) - max(bo_log$step) + 30
  n_seen <- min(n_init + step, length(dists_eval))
  best_idx <- which.min(dists_eval[seq_len(n_seen)])
  
  list(
    best_dist = dists_eval[best_idx],
    best_tour = paths_eval[[best_idx]]
  )
}

tab_progress <- rbindlist(lapply(progress_checkpoints, function(step) {
  incumbent <- best_tour_up_to(step)
  
  data.table(
    step = step,
    selected_dist = bo_log[step == !!step, selected_dist],
    best_so_far_dist = incumbent$best_dist,
    best_tour = tour_string(incumbent$best_tour)
  )
}), fill = TRUE)

tab_progress[, delta_best := best_so_far_dist - shift(best_so_far_dist)]
tab_progress[is.na(delta_best), delta_best := 0]

cat("\n=== Best route table ===\n")
print(tab_best_route)

cat("\n=== Top 10 unique tours ===\n")
print(top_tours[seq_len(min(10, .N)), .(rank = .I, distance, tour)])

cat("\n=== Progress table ===\n")
print(tab_progress)

p_single_convergence <- ggplot(bo_log, aes(x = step, y = best_dist)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.2, alpha = 0.8) +
  labs(
    x = "Iteration",
    y = "Best-so-far distance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    plot.title = element_blank()
  )

print(p_single_convergence)

############################################################
# 11. Best-route MDS visualization
############################################################

make_route_plot_data <- function(path, distance_matrix) {
  xy <- cmdscale(as.dist(distance_matrix), k = 2)
  
  pos_df <- data.table(
    city = rownames(xy),
    x = xy[, 1],
    y = xy[, 2]
  )
  
  seg_df <- pos_df[match(path, city)]
  seg_df[, next_city := shift(city, type = "lead", fill = city[1])]
  seg_df[, next_x := shift(x, type = "lead", fill = x[1])]
  seg_df[, next_y := shift(y, type = "lead", fill = y[1])]
  
  seg_df
}

route_df <- make_route_plot_data(best_path, distance_matrix)

p_best_route <- ggplot(route_df, aes(x = x, y = y)) +
  geom_segment(
    aes(xend = next_x, yend = next_y),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.7
  ) +
  geom_point(size = 2.4) +
  geom_text(aes(label = city), vjust = -0.7, size = 3.2) +
  labs(
    x = "MDS dimension 1",
    y = "MDS dimension 2"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    plot.title = element_blank()
  )

print(p_best_route)

############################################################
# 12. Multi-run study
############################################################

run_multirun_tsp_study <- function(
    R = 20,
    seeds = 1:20,
    n_init = 30,
    n_iter = 100,
    n_choice_sets = 40,
    K_choice = 5,
    n_cand = 120,
    kappa = 2,
    op = "2opt"
) {
  stopifnot(length(seeds) >= R)
  
  N_total <- n_init + n_iter
  
  methods <- c(
    "MNL-BO (EI)",
    "MNL-BO (UCB)",
    "RF-SMAC (EI)",
    "RF-SMAC (UCB)",
    "Simulated Annealing",
    "Restarted Hill Climbing",
    "2-opt Local Search",
    "Random Search"
  )
  
  logs_best <- lapply(methods, function(.) {
    matrix(NA_real_, nrow = R, ncol = N_total)
  })
  names(logs_best) <- methods
  
  finals <- data.table(
    method = rep(methods, each = R),
    run = rep(seq_len(R), times = length(methods)),
    final_best = NA_real_
  )
  
  for (r in seq_len(R)) {
    seed <- seeds[r]
    
    res_mnl_ei <- run_mnl_bo_tsp(
      n_init = n_init,
      n_iter = n_iter,
      n_choice_sets = n_choice_sets,
      K_choice = K_choice,
      n_cand = n_cand,
      kappa = kappa,
      acq_mode = "EI",
      op = op,
      seed = seed
    )
    logs_best[["MNL-BO (EI)"]][r, ] <- res_mnl_ei$log$best_dist
    finals[method == "MNL-BO (EI)" & run == r, final_best := tail(res_mnl_ei$log$best_dist, 1)]
    
    res_mnl_ucb <- run_mnl_bo_tsp(
      n_init = n_init,
      n_iter = n_iter,
      n_choice_sets = n_choice_sets,
      K_choice = K_choice,
      n_cand = n_cand,
      kappa = kappa,
      acq_mode = "UCB",
      op = op,
      seed = seed + 10000
    )
    logs_best[["MNL-BO (UCB)"]][r, ] <- res_mnl_ucb$log$best_dist
    finals[method == "MNL-BO (UCB)" & run == r, final_best := tail(res_mnl_ucb$log$best_dist, 1)]
    
    res_rf_ei <- run_rf_smac_tsp(
      n_init = n_init,
      n_iter = n_iter,
      n_cand = n_cand,
      kappa = kappa,
      acq_mode = "EI",
      op = op,
      seed = seed + 20000
    )
    logs_best[["RF-SMAC (EI)"]][r, ] <- res_rf_ei$log$best_dist
    finals[method == "RF-SMAC (EI)" & run == r, final_best := tail(res_rf_ei$log$best_dist, 1)]
    
    res_rf_ucb <- run_rf_smac_tsp(
      n_init = n_init,
      n_iter = n_iter,
      n_cand = n_cand,
      kappa = kappa,
      acq_mode = "UCB",
      op = op,
      seed = seed + 30000
    )
    logs_best[["RF-SMAC (UCB)"]][r, ] <- res_rf_ucb$log$best_dist
    finals[method == "RF-SMAC (UCB)" & run == r, final_best := tail(res_rf_ucb$log$best_dist, 1)]
    
    res_sa <- run_simulated_annealing_tsp(
      n_init = n_init,
      n_iter = n_iter,
      seed = seed + 40000
    )
    logs_best[["Simulated Annealing"]][r, ] <- res_sa$log$best_dist
    finals[method == "Simulated Annealing" & run == r, final_best := tail(res_sa$log$best_dist, 1)]
    
    res_rhc <- run_restarted_hillclimb_tsp(
      n_init = n_init,
      n_iter = n_iter,
      seed = seed + 50000
    )
    logs_best[["Restarted Hill Climbing"]][r, ] <- res_rhc$log$best_dist
    finals[method == "Restarted Hill Climbing" & run == r, final_best := tail(res_rhc$log$best_dist, 1)]
    
    res_twoopt <- run_twoopt_local_search_tsp(
      n_init = n_init,
      n_iter = n_iter,
      seed = seed + 60000
    )
    logs_best[["2-opt Local Search"]][r, ] <- res_twoopt$log$best_dist
    finals[method == "2-opt Local Search" & run == r, final_best := tail(res_twoopt$log$best_dist, 1)]
    
    res_random <- run_random_search_tsp(
      n_init = n_init,
      n_iter = n_iter,
      seed = seed + 70000
    )
    logs_best[["Random Search"]][r, ] <- res_random$log$best_dist
    finals[method == "Random Search" & run == r, final_best := tail(res_random$log$best_dist, 1)]
    
    if (r %% 5 == 0) {
      cat(sprintf("Completed %d/%d multi-run replications\n", r, R))
    }
  }
  
  summary_table <- finals[, .(
    mean_final = mean(final_best, na.rm = TRUE),
    sd_final = sd(final_best, na.rm = TRUE),
    best_min = min(final_best, na.rm = TRUE)
  ), by = method]
  
  summary_table[, method := factor(method, levels = methods)]
  summary_table <- summary_table[order(method)]
  summary_table[, method := as.character(method)]
  
  trajectory_table <- rbindlist(lapply(methods, function(method_name) {
    mat <- logs_best[[method_name]]
    
    data.table(
      step = seq_len(N_total),
      method = method_name,
      mean_best = colMeans(mat, na.rm = TRUE),
      sd_best = apply(mat, 2, sd, na.rm = TRUE)
    )
  }))
  
  trajectory_table[, method := factor(method, levels = methods)]
  trajectory_table <- trajectory_table[order(method, step)]
  trajectory_table[, method := as.character(method)]
  trajectory_table[, ymin := mean_best - sd_best]
  trajectory_table[, ymax := mean_best + sd_best]
  
  list(
    summary = summary_table,
    trajectory = trajectory_table,
    finals = finals,
    logs_best = logs_best
  )
}

multi_run <- run_multirun_tsp_study(
  R = 20,
  seeds = 1:20,
  n_init = 30,
  n_iter = 100,
  n_choice_sets = 40,
  K_choice = 5,
  n_cand = 120,
  kappa = 2,
  op = "2opt"
)

cat("\n=== Multi-run summary table ===\n")
print(multi_run$summary)

############################################################
# 13. Multi-run convergence plot
############################################################

method_colors <- c(
  "MNL-BO (EI)" = "#1b9e77",
  "MNL-BO (UCB)" = "#d95f02",
  "RF-SMAC (EI)" = "#7570b3",
  "RF-SMAC (UCB)" = "#e7298a",
  "Simulated Annealing" = "#66a61e",
  "Restarted Hill Climbing" = "#e6ab02",
  "2-opt Local Search" = "#a6761d",
  "Random Search" = "#666666"
)

p_multirun <- ggplot(
  multi_run$trajectory,
  aes(x = step, y = mean_best, color = method, fill = method)
) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "Function evaluations",
    y = "Best-so-far tour length"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_blank()
  )

print(p_multirun)

############################################################
# 15. More plots from TSP study
############################################################

extra_packages <- c(
  "plotly", "leaflet", "leaflet.extras2", "pheatmap",
  "reshape2", "viridis", "htmlwidgets"
)

install_missing_packages(extra_packages)

suppressPackageStartupMessages({
  library(plotly)
  library(leaflet)
  library(leaflet.extras2)
  library(pheatmap)
  library(reshape2)
  library(viridis)
  library(htmlwidgets)
})

############################################################
# 15.1 Interactive BO convergence plot
############################################################

p_interactive_convergence <- plot_ly(bo_log, x = ~step) %>%
  add_lines(
    y = ~best_dist,
    name = "Best distance",
    line = list(width = 3)
  ) %>%
  add_lines(
    y = ~selected_dist,
    name = "Selected distance",
    line = list(width = 2)
  ) %>%
  layout(
    title = "BO Convergence (Interactive)",
    xaxis = list(title = "Iteration"),
    yaxis = list(title = "Distance")
  )

print(p_interactive_convergence)

############################################################
# 15.2 Single-run convergence plot with title
############################################################

p_single_convergence_titled <- ggplot(bo_log, aes(x = step, y = best_dist)) +
  geom_line(linewidth = 1.2, color = "black") +
  labs(
    title = "TSP: Best-so-far Distance (Single Run)",
    x = "Iteration",
    y = "Best-so-far distance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  )

print(p_single_convergence_titled)

# ggsave(
#   "plots/tsp_single_run_best_so_far.png",
#   p_single_convergence_titled,
#   width = 8,
#   height = 5,
#   dpi = 300
# )

############################################################
# 15.3 Leaflet route map with manual coordinates
############################################################

city_coordinates <- data.frame(
  city = c(
    "Atlanta", "Boston", "Chicago", "Dallas", "Denver", "Houston",
    "Las Vegas", "Los Angeles", "Miami", "New Orleans", "New York",
    "Phoenix", "San Francisco", "Seattle", "Washington D.C."
  ),
  lat = c(
    33.7490, 42.3601, 41.8781, 32.7767, 39.7392, 29.7604,
    36.1699, 34.0522, 25.7617, 29.9511, 40.7128,
    33.4484, 37.7749, 47.6062, 38.9072
  ),
  lng = c(
    -84.3880, -71.0589, -87.6298, -96.7970, -104.9903, -95.3698,
    -115.1398, -118.2437, -80.1918, -90.0715, -74.0060,
    -112.0740, -122.4194, -122.3321, -77.0369
  )
)

route_coords <- city_coordinates[match(best_path, city_coordinates$city), ]
route_coords_closed <- rbind(route_coords, route_coords[1, ])

leaflet_route <- leaflet() %>%
  addTiles(urlTemplate = "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png") %>%
  addMarkers(
    data = city_coordinates,
    lng = ~lng,
    lat = ~lat,
    popup = ~city,
    label = ~city
  ) %>%
  setView(
    lng = mean(city_coordinates$lng),
    lat = mean(city_coordinates$lat),
    zoom = 4
  ) %>%
  addPolylines(
    data = route_coords_closed,
    lng = ~lng,
    lat = ~lat,
    weight = 3,
    opacity = 0.8,
    color = "blue"
  )

for (i in seq_len(nrow(route_coords_closed) - 1)) {
  segment_i <- route_coords_closed[c(i, i + 1), ]
  
  leaflet_route <- leaflet_route %>%
    addArrowhead(
      data = segment_i,
      lng = ~lng,
      lat = ~lat,
      color = "blue",
      weight = 2,
      opacity = 0.7,
      options = arrowheadOptions(size = 0.12)
    )
}

leaflet_route <- leaflet_route %>%
  addAwesomeMarkers(
    lng = route_coords$lng[1],
    lat = route_coords$lat[1],
    popup = paste0("START: ", route_coords$city[1]),
    label = paste0("START: ", route_coords$city[1])
  )

leaflet_route

# htmlwidgets::saveWidget(
#   leaflet_route,
#   "plots/tsp_best_route_leaflet.html",
#   selfcontained = TRUE
# )

############################################################
# 15.4 Distance matrix heatmap with clustering
############################################################

pheatmap(
  as.matrix(distance_matrix),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  main = "",
  fontsize = 14,
  fontsize_row = 14,
  fontsize_col = 14,
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
  angle_col = 45,
  border_color = "grey70"
)

distance_melted <- melt(distance_matrix)
colnames(distance_melted) <- c("City1", "City2", "Distance")

p_distance_heatmap <- ggplot(distance_melted, aes(x = City1, y = City2, fill = Distance)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL, fill = "Distance") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    plot.title = element_blank()
  )

print(p_distance_heatmap)

# ggsave(
#   "plots/tsp_distance_matrix_heatmap.png",
#   p_distance_heatmap,
#   width = 10,
#   height = 7,
#   dpi = 300
# )

############################################################
# 15.5 Base-R multi-run plot: four-method version
############################################################

plot_four_method_multirun <- function(traj_dt) {
  dt <- copy(as.data.table(traj_dt))
  
  method_order <- c(
    "MNL-BO (EI)",
    "MNL-BO (UCB)",
    "Random Search",
    "2-opt Local Search"
  )
  
  dt <- dt[method %in% method_order]
  dt[, method := factor(method, levels = method_order)]
  dt[, ymin := mean_best - sd_best]
  dt[, ymax := mean_best + sd_best]
  
  colors <- c(
    "MNL-BO (EI)" = "#1b9e77",
    "MNL-BO (UCB)" = "#d95f02",
    "Random Search" = "#7570b3",
    "2-opt Local Search" = "#e7298a"
  )
  
  x <- sort(unique(dt$step))
  y_min <- min(dt$ymin, na.rm = TRUE)
  y_max <- max(dt$ymax, na.rm = TRUE)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mar = c(5.2, 5.2, 1.2, 1.2))
  
  plot(
    NA,
    xlim = range(x),
    ylim = c(y_min, y_max),
    xlab = "Function evaluations",
    ylab = "Best-so-far tour length",
    main = ""
  )
  
  for (m in method_order) {
    d <- dt[method == m][order(step)]
    
    polygon(
      c(d$step, rev(d$step)),
      c(d$ymin, rev(d$ymax)),
      border = NA,
      col = adjustcolor(colors[m], alpha.f = 0.15)
    )
    
    lines(
      d$step,
      d$mean_best,
      lwd = 2.8,
      col = colors[m]
    )
  }
  
  legend(
    "topright",
    legend = method_order,
    col = colors[method_order],
    lwd = 3,
    bty = "n",
    cex = 0.95
  )
}

plot_four_method_multirun(multi_run$trajectory)

############################################################
# 15.6 Base-R multi-run plot: all-method version
############################################################

plot_all_method_multirun_base <- function(traj_dt) {
  dt <- copy(as.data.table(traj_dt))
  
  method_order <- c(
    "MNL-BO (EI)",
    "MNL-BO (UCB)",
    "RF-SMAC (EI)",
    "RF-SMAC (UCB)",
    "Simulated Annealing",
    "Restarted Hill Climbing",
    "2-opt Local Search",
    "Random Search"
  )
  
  method_order <- intersect(method_order, unique(dt$method))
  dt <- dt[method %in% method_order]
  dt[, method := factor(method, levels = method_order)]
  dt[, ymin := mean_best - sd_best]
  dt[, ymax := mean_best + sd_best]
  
  colors <- c(
    "MNL-BO (EI)" = "#1b9e77",
    "MNL-BO (UCB)" = "#d95f02",
    "RF-SMAC (EI)" = "#7570b3",
    "RF-SMAC (UCB)" = "#e7298a",
    "Simulated Annealing" = "#66a61e",
    "Restarted Hill Climbing" = "#e6ab02",
    "2-opt Local Search" = "#a6761d",
    "Random Search" = "#666666"
  )
  
  x <- sort(unique(dt$step))
  y_min <- min(dt$ymin, na.rm = TRUE)
  y_max <- max(dt$ymax, na.rm = TRUE)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(
    mar = c(5.5, 6, 1.5, 1.5),   # more space for bigger labels
    mgp = c(3.5, 1.2, 0)
  )
  
  plot(
    NA,
    xlim = range(x),
    ylim = c(y_min, y_max),
    xlab = "Function evaluations",
    ylab = "Best-so-far tour length",
    main = "",
    yaxt = "n",
    cex.lab = 2.2,  
    cex.axis = 2.0 
  )
  
  # Custom Y axis with commas
  y_ticks <- axTicks(2)
  axis(
    side = 2,
    at = y_ticks,
    labels = ifelse(
      abs(y_ticks) >= 10000,
      format(y_ticks, big.mark = ",", scientific = FALSE, trim = TRUE),
      as.character(y_ticks)
    ),
    cex.axis = 1.4
  )
  
  for (m in method_order) {
    d <- dt[method == m][order(step)]
    
    polygon(
      c(d$step, rev(d$step)),
      c(d$ymin, rev(d$ymax)),
      border = NA,
      col = adjustcolor(colors[m], alpha.f = 0.15)
    )
    
    lines(
      d$step,
      d$mean_best,
      lwd = 3.2,   # thicker for print
      col = colors[m]
    )
  }
  
  legend(
    "topright",
    legend = method_order,
    col = colors[method_order],
    lwd = 4,
    bty = "n",
    cex = 1.4   
  )
}

plot_all_method_multirun_base(multi_run$trajectory)

############################################################
# 15.7 ggplot multi-run plot: all methods
############################################################

p_multirun_ggplot_extended <- ggplot(
  multi_run$trajectory,
  aes(x = step, y = mean_best, color = method, fill = method)
) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  labs(
    x = "Function evaluations",
    y = "Best-so-far tour length"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.title = element_blank()
  )

print(p_multirun_ggplot_extended)

# ggsave(
#   "plots/tsp_multirun_all_methods.png",
#   p_multirun_ggplot_extended,
#   width = 12,
#   height = 6,
#   dpi = 300
# )

############################################################
# 15.8 Export all final tables
############################################################

cat("\n=== Extended output tables ===\n")

cat("\nBest route table:\n")
print(tab_best_route)

cat("\nTop 10 unique tours:\n")
print(top_tours[seq_len(min(10, .N)), .(rank = .I, distance, tour)])

cat("\nProgress table:\n")
print(tab_progress)

cat("\nMulti-run summary table:\n")
print(multi_run$summary)

cat("\nExtended plots completed.\n")

