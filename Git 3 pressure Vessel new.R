#===============================================================================
# MNL-Based Bayesian Optimization for Pressure Vessel Design
# Mixed discrete, continuous, and categorical variables
#
# Main method:
#   Pairwise preference Bayesian optimization using MNL / Bradley-Terry choice model
#
# Baseline:
#   Random-forest SMAC-style optimization using ranger
#
# Objective:
#   Minimize total pressure-vessel cost with constraint penalties
#===============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mlogit)
  library(ggplot2)
})

set.seed(123)

#-------------------------
# 1) Materials + properties
#-------------------------
materials <- c("Steel A", "Steel B", "Steel C")

mat_tbl <- data.table(
  material = materials,
  rho = c(0.2830, 0.2836, 0.2844),   # lb/in^3
  Cs  = c(0.35,   0.30,   0.40),     # $/lb rolled
  Ch  = c(1.00,   1.05,   0.95),     # $/lb forged
  Cw  = c(8.00,   8.50,   7.50)      # $/lb welding (shell material only)
)

# convenient map
mat_props <- function(y) mat_tbl[match(y, material)]

#-------------------------
# 2) Decision encoding + bounds
#    n_s, n_h are integer indices -> delta = 0.0625 * n
#-------------------------
delta_step <- 0.0625
ns_min <- 1
ns_max <- as.integer(round(2 / delta_step)) 
nh_min <- 1
nh_max <- as.integer(round(2 / delta_step))

R_min <- 1e-6
R_max <- 100

L_min <- 1e-6
L_max <- 240

#-------------------------
# 3) Objective (minimization): cost + constraints
#-------------------------
pressure_vessel_cost <- function(ns, nh, R, L, ys, yh, penalty = 1e6) {
  # map to thicknesses
  ds <- delta_step * ns
  dh <- delta_step * nh
  
  # material props
  shell <- mat_props(ys)
  head  <- mat_props(yh)
  
  rho_s <- shell$rho; Cs <- shell$Cs; Cw <- shell$Cw
  rho_h <- head$rho;  Ch <- head$Ch
  
  # Constraints (as in your formulation)
  g1 <- 0.0193 * R - ds
  g2 <- 0.00954 * R - dh
  g3 <- 1.296e6 - pi * R^2 * L - (4/3) * pi * R^3
  
  feasible <- (g1 <= 0) && (g2 <= 0) && (g3 <= 0) &&
    (ns >= ns_min) && (ns <= ns_max) &&
    (nh >= nh_min) && (nh <= nh_max) &&
    (R > 0) && (R <= R_max) &&
    (L > 0) && (L <= L_max) &&
    (ys %in% materials) && (yh %in% materials)
  cost <- 2*pi*R*(rho_s*Cs*ds*L + rho_h*Ch*dh*R) +
    (4*pi*ds^2*rho_s*Cw/9) * (L + 2*pi*R)
  
  if (!is.finite(cost)) return(penalty)
  if (!feasible) return(penalty + abs(g1) + abs(g2) + abs(g3))
  
  cost
}

# Utility for MNL (maximize utility)
utility <- function(cost) -cost

#-------------------------
# 4) Helpers: represent a design, clamp/repair, and pretty printing
#-------------------------
make_design <- function(ns, nh, R, L, ys, yh) {
  list(ns = as.integer(ns), nh = as.integer(nh),
       R = as.numeric(R), L = as.numeric(L),
       ys = as.character(ys), yh = as.character(yh))
}

repair_design <- function(x) {
  x$ns <- max(ns_min, min(ns_max, as.integer(round(x$ns))))
  x$nh <- max(nh_min, min(nh_max, as.integer(round(x$nh))))
  x$R  <- max(R_min,  min(R_max,  as.numeric(x$R)))
  x$L  <- max(L_min,  min(L_max,  as.numeric(x$L)))
  if (!(x$ys %in% materials)) x$ys <- sample(materials, 1)
  if (!(x$yh %in% materials)) x$yh <- sample(materials, 1)
  x
}

design_key <- function(x) {
  paste(x$ns, x$nh, sprintf("%.6f", x$R), sprintf("%.6f", x$L), x$ys, x$yh, sep="|")
}

design_to_row <- function(x) {
  data.frame(
    ns = x$ns,
    nh = x$nh,
    R  = x$R,
    L  = x$L,
    ys = factor(x$ys, levels = materials),
    yh = factor(x$yh, levels = materials),
    stringsAsFactors = FALSE
  )
}

#-------------------------
# 5) Feature encoding for MNL (mixed variables)
#    - numeric: ns, nh, R, L (with mild scaling)
#    - categorical: ys, yh (one-hot)
#    - interaction between ys and yh can be added
#-------------------------
make_feature_matrix <- function(df) {
  df2 <- copy(as.data.table(df))
  df2[, ns_sc := ns / ns_max]
  df2[, nh_sc := nh / nh_max]
  df2[, R_sc  := R / R_max]
  df2[, L_sc  := L / L_max]
  
  mm <- model.matrix(
    ~ 0 + ns_sc + nh_sc + R_sc + L_sc + ys + yh + ys:yh,
    data = df2
  )
  colnames(mm) <- make.names(colnames(mm), unique = TRUE)
  
  mm
}

#-------------------------
# 6) Pairwise preference dataset
#    Each pair is a choice situation with 2 alternatives:
#      winner = lower cost (higher utility), loser = higher cost
#    We fit mlogit as conditional logit:
#      choice ~ 0 + features
#-------------------------
make_pairwise_long <- function(ids_a, ids_b, X_mm, costs, chid_start = 1) {
  # ids_a, ids_b are equal-length integer vectors of indices into archive
  stopifnot(length(ids_a) == length(ids_b))
  m <- length(ids_a)
  
  rows <- vector("list", 2*m)
  k <- 1
  
  for (i in seq_len(m)) {
    a <- ids_a[i]
    b <- ids_b[i]
    # winner has smaller cost
    win <- if (costs[a] <= costs[b]) a else b
    los <- if (win == a) b else a
    
    # two-alternative choice set
    chid <- chid_start + i - 1
    
    # winner row (choice=1)
    rows[[k]] <- data.frame(
      chid = chid,
      alt  = paste0("x", win),
      choice = 1L,
      X_mm[win, , drop = FALSE],
      check.names = FALSE
    ); k <- k + 1
    
    # loser row (choice=0)
    rows[[k]] <- data.frame(
      chid = chid,
      alt  = paste0("x", los),
      choice = 0L,
      X_mm[los, , drop = FALSE],
      check.names = FALSE
    ); k <- k + 1
  }
  
  long_df <- rbindlist(rows, fill = TRUE)
  long_df[]
}

fit_mnl_pairs <- function(long_df) {
  
  # ✅ ensure feature column names are safe
  feat_cols <- setdiff(colnames(long_df), c("chid","alt","choice"))
  safe_cols <- make.names(feat_cols, unique = TRUE)
  
  # rename in the data.frame (mlogit uses names internally)
  setnames(long_df, feat_cols, safe_cols)
  
  mdat <- mlogit.data(
    long_df,
    choice   = "choice",
    shape    = "long",
    chid.var = "chid",
    alt.var  = "alt"
  )
  
  fml <- as.formula(paste("choice ~ 0 +", paste(safe_cols, collapse = " + ")))
  
  tryCatch(mlogit(fml, data = mdat), error = function(e) NULL)
}

# Predict mu and sigma on ANY set of designs using the same feature template
predict_mu_sigma <- function(mod, X_mm_new) {
  n <- nrow(X_mm_new)
  if (is.null(mod)) {
    return(list(mu = rep(0, n), sigma = rep(1, n)))
  }
  
  beta <- coef(mod)
  V <- tryCatch(vcov(mod), error = function(e) NULL)
  
  # Align columns
  cols <- colnames(X_mm_new)
  beta_full <- rep(0, length(cols)); names(beta_full) <- cols
  beta_full[names(beta)] <- beta
  
  mu <- as.numeric(X_mm_new %*% beta_full)
  
  if (is.null(V)) {
    sigma <- rep(1, n)
  } else {
    V_full <- matrix(0, nrow = length(cols), ncol = length(cols),
                     dimnames = list(cols, cols))
    V_full[rownames(V), colnames(V)] <- V
    XV <- X_mm_new %*% V_full
    var <- rowSums(XV * X_mm_new)
    sigma <- sqrt(pmax(var, 1e-12))
  }
  
  list(mu = mu, sigma = sigma)
}

#-------------------------
# 7) Candidate generation (local + global)
#    - local: around incumbent best
#    - global: random within bounds
#-------------------------
mutate_design_local <- function(x, p_mat = 0.25, p_thick = 0.40,
                                sd_R = 6, sd_L = 15) {
  y <- x
  
  # jitter continuous
  y$R <- y$R + rnorm(1, 0, sd_R)
  y$L <- y$L + rnorm(1, 0, sd_L)
  
  # occasionally change thickness indices
  if (runif(1) < p_thick) y$ns <- y$ns + sample(c(-2,-1,1,2), 1)
  if (runif(1) < p_thick) y$nh <- y$nh + sample(c(-2,-1,1,2), 1)
  
  # occasionally flip materials
  if (runif(1) < p_mat) y$ys <- sample(materials, 1)
  if (runif(1) < p_mat) y$yh <- sample(materials, 1)
  
  repair_design(y)
}

random_design_global <- function() {
  make_design(
    ns = sample(ns_min:ns_max, 1),
    nh = sample(nh_min:nh_max, 1),
    R  = runif(1, 1, R_max),
    L  = runif(1, 1, L_max),
    ys = sample(materials, 1),
    yh = sample(materials, 1)
  )
}

generate_candidates <- function(best_x, n_local = 80, n_global = 40) {
  locals <- replicate(n_local, mutate_design_local(best_x), simplify = FALSE)
  globals <- replicate(n_global, random_design_global(), simplify = FALSE)
  
  all <- c(locals, globals)
  
  # deduplicate
  keys <- vapply(all, design_key, FUN.VALUE = character(1))
  all[!duplicated(keys)]
}

#-------------------------
# 8) Acquisition functions (utility scale, maximize)
#-------------------------
acq_ucb <- function(mu, sigma, kappa = 2.0) mu + kappa * sigma
acq_ei  <- function(mu, best_mu) pmax(mu - best_mu, 0)

#===============================================================================
# 9) BO loop (pairwise MNL)
#===============================================================================
run_pv_mnl_bo <- function(
    n_init = 50,
    n_iter = 150,
    n_pairs = 80,       # number of pairwise comparisons per refit
    n_cand = 120,
    acq_mode = c("UCB","EI"),
    kappa = 2.0,
    seed = 1
) {
  set.seed(seed)
  acq_mode <- match.arg(acq_mode)
  
  # Initial archive
  X_list <- replicate(n_init, random_design_global(), simplify = FALSE)
  X_df <- rbindlist(lapply(X_list, design_to_row), fill = TRUE)
  costs <- vapply(
    seq_len(nrow(X_df)),
    function(i) {
      pressure_vessel_cost(
        ns = X_df$ns[i],
        nh = X_df$nh[i],
        R  = X_df$R[i],
        L  = X_df$L[i],
        ys = as.character(X_df$ys[i]),
        yh = as.character(X_df$yh[i])
      )
    },
    numeric(1)
  )
  
  utils <- utility(costs)
  
  # incumbent (min cost)
  best_idx <- which.min(costs)
  best_cost <- costs[best_idx]
  best_x <- X_list[[best_idx]]
  best_mu <- max(utils)
  
  # log
  bo_log <- data.table(iter = integer(), selected_cost = numeric(), best_cost = numeric())
  
  # cache (optional)
  seen <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(X_df))) {
    seen[[design_key(X_list[[i]])]] <- costs[i]
  }
  
  for (t in seq_len(n_iter)) {
    
    #------------------------------------------------------------
    # (A) Build pairwise training data from the current archive
    #------------------------------------------------------------
    N <- nrow(X_df)
    if (N < 4) {
      mod <- NULL
    } else {
      # Sample random pairs from archive
      pairs_a <- sample.int(N, size = n_pairs, replace = TRUE)
      pairs_b <- sample.int(N, size = n_pairs, replace = TRUE)
      # avoid identical pairs
      same <- which(pairs_a == pairs_b)
      if (length(same) > 0) {
        pairs_b[same] <- ((pairs_b[same] + sample(1:(N-1), length(same), replace=TRUE) - 1) %% N) + 1
      }
      
      # Feature matrix on archive
      X_mm <- make_feature_matrix(X_df)
      
      # Pairwise long data
      long_df <- make_pairwise_long(pairs_a, pairs_b, X_mm, costs, chid_start = 1L)
      
      # Fit MNL on pairs
      mod <- fit_mnl_pairs(long_df)
    }
    
    #------------------------------------------------------------
    # (B) Generate candidates (local + global)
    #------------------------------------------------------------
    cand_list <- generate_candidates(best_x, n_local = 90, n_global = 30)
    if (length(cand_list) > n_cand) cand_list <- cand_list[1:n_cand]
    cand_df <- rbindlist(lapply(cand_list, design_to_row), fill = TRUE)
    
    #------------------------------------------------------------
    # (C) Select next point by acquisition (fallback random if fail)
    #------------------------------------------------------------
    if (is.null(mod)) {
      pick <- sample.int(nrow(cand_df), 1)
    } else {
      X_cand_mm <- make_feature_matrix(cand_df)
      pred <- predict_mu_sigma(mod, X_cand_mm)
      mu <- pred$mu
      sigma <- pred$sigma
      
      if (acq_mode == "UCB") {
        a <- acq_ucb(mu, sigma, kappa = kappa)
      } else {
        a <- acq_ei(mu, best_mu)
      }
      pick <- which.max(a)
    }
    
    x_next <- cand_list[[pick]]
    k <- design_key(x_next)
    
    #------------------------------------------------------------
    # (D) Evaluate objective (with cache)
    #------------------------------------------------------------
    if (!is.null(seen[[k]])) {
      c_next <- seen[[k]]
    } else {
      c_next <- pressure_vessel_cost(x_next$ns, x_next$nh, x_next$R, x_next$L, x_next$ys, x_next$yh)
      seen[[k]] <- c_next
    }
    u_next <- utility(c_next)
    
    #------------------------------------------------------------
    # (E) Update archive
    #------------------------------------------------------------
    X_list <- c(X_list, list(x_next))
    X_df <- rbind(X_df, as.data.table(design_to_row(x_next)), fill = TRUE)
    costs <- c(costs, c_next)
    utils <- c(utils, u_next)
    
    # incumbent update
    if (c_next < best_cost) {
      best_cost <- c_next
      best_x <- x_next
    }
    best_mu <- max(best_mu, u_next)
    
    bo_log <- rbind(bo_log, data.table(iter = t, selected_cost = c_next, best_cost = best_cost))
    
    if (t %% 20 == 0) {
      cat(sprintf("Iter %d | Best cost: %.2f\n", t, best_cost))
    }
  }
  
  list(
    best_x = best_x,
    best_cost = best_cost,
    history = bo_log,
    archive = X_df,
    archive_cost = costs
  )
}

#===============================================================================
# 10) Run a single experiment
#===============================================================================
res <- run_pv_mnl_bo(
  n_init = 50,
  n_iter = 150,
  n_pairs = 80,
  n_cand = 120,
  acq_mode = "UCB",   # "UCB" or "EI"
  kappa = 2.0,
  seed = 7
)

cat("\n=== BEST DESIGN FOUND ===\n")
print(res$best_x)
cat(sprintf("Best cost: %.2f\n", res$best_cost))

#===============================================================================
# 11) outputs: tables + figures
#===============================================================================

bo_log <- as.data.table(res$history)
setorder(bo_log, iter)

# ---- Table A: progress checkpoints (cost should decrease; lower is better)
iters_show <- unique(c(1, 5, 10, 25, 50, 75, 100, 125, 150))
iters_show <- iters_show[iters_show <= max(bo_log$iter)]

tab_progress <- bo_log[iter %in% iters_show,
                       .(iter, selected_cost, best_so_far_cost = best_cost)]

cat("\n=== TABLE: Pressure Vessel Progress Over Iterations ===\n")
print(tab_progress)

# ---- Table B: best design (single run)
best_row <- res$best_x
tab_best <- data.table(
  best_cost = res$best_cost,
  ns = best_row$ns, nh = best_row$nh,
  delta_s = delta_step * best_row$ns,
  delta_h = delta_step * best_row$nh,
  R = best_row$R, L = best_row$L,
  ys = best_row$ys, yh = best_row$yh
)

cat("\n=== TABLE: Best Pressure Vessel Design (single run) ===\n")
print(tab_best)

# ---- Figure 1: single-run convergence (best-so-far cost)
p_best_single <- ggplot(bo_log, aes(x = iter, y = best_cost)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Pressure Vessel: Best-so-far Cost (Single Run)",
    x = "Iteration",
    y = "Best-so-far total cost"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
print(p_best_single)

# ---- Figure 2: selected cost vs iteration (shows exploration)
p_sel <- ggplot(bo_log, aes(x = iter, y = selected_cost)) +
  geom_line(alpha = 0.85) +
  labs(
    title = "Pressure Vessel: Selected Cost per Iteration",
    x = "Iteration",
    y = "Selected total cost"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
print(p_sel)

# 1) Run 10 independent runs (same settings)
# -------------------------
# MULTI-RUN EXPERIMENTS
# -------------------------
seeds <- 1:10

runs_list <- lapply(seeds, function(s) {
  out <- run_pv_mnl_bo(
    n_init = 50,
    n_iter = 150,
    n_pairs = 80,
    n_cand = 120,
    acq_mode = "UCB",
    kappa = 2.0,
    seed = s
  )
  dt <- as.data.table(out$history)
  dt[, seed := s]
  list(res = out, log = dt)
})

logs_all <- rbindlist(lapply(runs_list, `[[`, "log"), fill = TRUE)
setorder(logs_all, seed, iter)

# 2) Table — Summary over runs (mean ± sd)
# initial best and final best per run
init_dt <- logs_all[iter == min(iter), .(init_best = best_cost[1]), by = seed]
final_dt <- logs_all[, .SD[which.max(iter)], by = seed]
setnames(final_dt, "best_cost", "final_best")

tab_summary <- merge(init_dt, final_dt[, .(seed, final_best)], by = "seed")
tab_summary[, improvement := final_best - init_best]  # negative = improvement

tab_summary_paper <- tab_summary[, .(
  mean_init_best  = mean(init_best),
  sd_init_best    = sd(init_best),
  mean_final_best = mean(final_best),
  sd_final_best   = sd(final_best),
  mean_improvement = mean(improvement),
  sd_improvement   = sd(improvement)
)]

cat("\n=== TABLE: Pressure Vessel summary over runs (mean ± sd) ===\n")
print(tab_summary_paper)

# 3) Figure — Mean best-so-far curve (mean ± 1 sd ribbon)
curve_dt <- logs_all[, .(
  mean_best = mean(best_cost),
  sd_best   = sd(best_cost)
), by = iter]

curve_dt[, `:=`(
  upper = mean_best + sd_best,
  lower = pmax(mean_best - sd_best, 0)
)]

p_mean_curve <- ggplot(curve_dt, aes(x = iter, y = mean_best)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  labs(
    title = "Pressure Vessel: Mean Best-so-far Cost Over Runs",
    x = "Iteration",
    y = "Best-so-far total cost"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_mean_curve)

# 4) Table — infeasible rate over the whole run
penalty <- 1e6

# using the archive costs from each run
tab_feas <- rbindlist(lapply(runs_list, function(x) {
  costs <- x$res$archive_cost
  data.table(
    infeasible_rate = mean(costs >= penalty),
    feasible_best   = min(costs[costs < penalty], na.rm = TRUE)
  )
}), fill = TRUE)

cat("\n=== TABLE (optional): infeasible rate over runs ===\n")
print(tab_feas[, .(
  mean_infeasible = mean(infeasible_rate),
  sd_infeasible   = sd(infeasible_rate),
  mean_best_feasible = mean(feasible_best, na.rm = TRUE),
  sd_best_feasible   = sd(feasible_best, na.rm = TRUE)
)])

# 5) Figure — cumulative infeasible fraction (single run)
penalty <- 1e6
bo_log <- as.data.table(res$history)
bo_log[, infeasible := selected_cost >= penalty]
bo_log[, cum_infeasible := cumsum(infeasible) / .I]

p_infeas <- ggplot(bo_log, aes(x = iter, y = cum_infeasible)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Pressure Vessel: Cumulative Infeasible Fraction (Single Run)",
    x = "Iteration",
    y = "Fraction infeasible so far"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_infeas)

#####################
# 0) Setup (themes + palettes)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

theme_paper <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

## A) Single-run plots (strong for “method behavior”)
# 1) Best-so-far cost (single run)
bo_log <- as.data.table(res$history)
setorder(bo_log, iter)

p_best_single <- ggplot(bo_log, aes(iter, best_cost)) +
  geom_line(linewidth = 1.2, color = "#1f77b4") +
  labs(
    title = "Pressure Vessel: Best-so-far Cost (Single Run)",
    x = "Iteration", y = "Best-so-far total cost"
  ) +
  theme_paper()

print(p_best_single)

## B) Multi-run plots
# 4) Mean best-so-far ± 1 sd ribbon (robustness)
curve_dt <- logs_all[, .(
  mean_best = mean(best_cost),
  sd_best   = sd(best_cost)
), by = iter]
curve_dt[, `:=`(
  upper = mean_best + sd_best,
  lower = pmax(mean_best - sd_best, 0)
)]

p_mean_curve <- ggplot(curve_dt, aes(iter, mean_best)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#1f77b4", alpha = 0.2) +
  geom_line(linewidth = 1.2, color = "#1f77b4") +
  labs(
    title = "Pressure Vessel: Mean Best-so-far Cost Over Runs",
    x = "Iteration", y = "Best-so-far total cost"
  ) +
  theme_paper()

print(p_mean_curve)

# 5) “Spaghetti plot” of best-so-far across runs (with light lines)
p_spaghetti <- ggplot(logs_all, aes(iter, best_cost, group = seed)) +
  geom_line(color = "#1f77b4", alpha = 0.18, linewidth = 0.9) +
  geom_line(data = curve_dt, aes(iter, mean_best),
            inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  labs(
    title = "Pressure Vessel: Best-so-far Trajectories (10 Runs)",
    x = "Iteration", y = "Best-so-far total cost"
  ) +
  theme_paper()

print(p_spaghetti)

comma_5plus <- function(x) {
  ifelse(abs(x) >= 10000,
         format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
         as.character(round(x)))
}

p_spaghetti <- ggplot(logs_all, aes(iter, best_cost, group = seed)) +
  geom_line(color = "#377EB8", alpha = 0.18, linewidth = 0.9) +
  geom_line(
    data = curve_dt,
    aes(iter, mean_best),
    inherit.aes = FALSE,
    color = "#E41A1C",
    linewidth = 1.1
  ) +
  scale_y_continuous(labels = comma_5plus) +
  labs(
    x = "Iteration",
    y = "Best-so-far total cost"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text  = element_text(size = 15),
    plot.title = element_blank()
  )

print(p_spaghetti)

# 6) Distribution of final best costs (violin + box) — very paper-friendly
final_dt <- logs_all[, .SD[which.max(iter)], by = seed]

p_final_dist <- ggplot(final_dt, aes(x = "Final", y = best_cost)) +
  geom_violin(fill = "#9467bd", alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.5, fill = "#9467bd", alpha = 0.7) +
  geom_jitter(width = 0.06, alpha = 0.6, size = 1.7, color = "black") +
  labs(
    title = "Pressure Vessel: Final Best Cost Distribution (10 Runs)",
    x = "", y = "Final best cost"
  ) +
  theme_paper()

print(p_final_dist)

comma_5plus <- function(x) {
  ifelse(abs(x) >= 10000,
         format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
         as.character(round(x)))
}

final_dt <- logs_all[, .SD[which.max(iter)], by = seed]

p_final_dist <- ggplot(final_dt, aes(x = "Final", y = best_cost)) +
  geom_violin(
    fill = "#984EA3",
    alpha = 0.35,
    color = NA
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.alpha = 0.5,
    fill = "#984EA3",
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.06,
    alpha = 0.6,
    size = 1.7,
    color = "black"
  ) +
  scale_y_continuous(labels = comma_5plus) +
  labs(
    x = "",
    y = "Final best cost"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text  = element_text(size = 15),
    plot.title = element_blank()
  )

print(p_final_dist)

# 7) Cumulative infeasible fraction (single run)
penalty <- 1e6
bo_log[, infeasible := selected_cost >= penalty]
bo_log[, cum_infeasible := cumsum(infeasible) / .I]

p_infeas <- ggplot(bo_log, aes(iter, cum_infeasible)) +
  geom_line(linewidth = 1.2, color = "#d62728") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Pressure Vessel: Cumulative Infeasible Fraction (Single Run)",
    x = "Iteration", y = "Fraction infeasible so far"
  ) +
  theme_paper()

print(p_infeas)

# 8) Selected costs with penalty threshold line (single run)
p_penalty_line <- ggplot(bo_log, aes(iter, selected_cost)) +
  geom_point(alpha = 0.65, size = 1.6, color = "#ff7f0e") +
  geom_hline(yintercept = penalty, linetype = "dashed", color = "#d62728") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Pressure Vessel: Selected Costs and Penalty Threshold",
    x = "Iteration", y = "Selected cost (penalty dashed line)"
  ) +
  theme_paper()

print(p_penalty_line)

## D) “Design insight” plots (excellent for mixed categorical + continuous)
# 9) Material pair frequency among best designs (across runs)
# best design of each run
best_runs <- rbindlist(lapply(seq_along(runs_list), function(i) {
  out <- runs_list[[i]]$res
  x <- out$best_x
  data.table(seed = seeds[i], ys = x$ys, yh = x$yh, best_cost = out$best_cost)
}), fill = TRUE)

# 10) Scatter of evaluated designs: (R, L) colored by feasibility (single run archive)
archive_dt <- as.data.table(res$archive)
archive_dt[, cost := res$archive_cost]
archive_dt[, feasible := cost < penalty]

p_RL <- ggplot(archive_dt, aes(R, L)) +
  geom_point(aes(color = feasible), alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"),
                     labels = c("TRUE" = "Feasible", "FALSE" = "Penalized")) +
  labs(
    title = "Pressure Vessel: Evaluated Designs in (R, L) Space",
    x = "Radius R", y = "Length L", color = ""
  ) +
  theme_paper()

print(p_RL)

# 11) Cost vs R (or L) colored by shell material (single run, feasible only)
archive_feas <- archive_dt[feasible == TRUE]

p_cost_R <- ggplot(archive_feas, aes(R, cost, color = ys)) +
  geom_point(alpha = 0.75, size = 1.7) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Pressure Vessel: Feasible Costs vs Radius (colored by shell material)",
    x = "Radius R", y = "Total cost", color = "Shell material"
  ) +
  theme_paper()

print(p_cost_R)

# 12) Thickness indices / thickness values distribution (single run, feasible)
archive_feas[, delta_s := delta_step * ns]
archive_feas[, delta_h := delta_step * nh]

p_thickness <- ggplot(melt(archive_feas[, .(delta_s, delta_h)]), aes(x = value, fill = variable)) +
  geom_histogram(bins = 18, alpha = 0.55, position = "identity") +
  labs(
    title = "Pressure Vessel: Thickness Distributions (Feasible Designs)",
    x = "Thickness (in)", y = "Count", fill = ""
  ) +
  theme_paper()

print(p_thickness)

#===============================================================================
# PATCH B: SMAC-like Random Forest (RF) baseline for Pressure Vessel
# - Surrogate: Random Forest via ranger
# - Uncertainty: tree-wise sd
# - Acquisition: EI (maximize) on cost minimization, or LCB
# - Multi-run: mean±sd + LaTeX table
#===============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
})

#-------------------------
# RF utilities
#-------------------------
rf_predict_mean_sd <- function(rf_model, X_new) {
  pr <- predict(rf_model, data = X_new, predict.all = TRUE)
  all <- pr$predictions
  list(mu = rowMeans(all), sd = apply(all, 1, sd))
}

ei_min <- function(mu, sd, best) {
  sd <- pmax(sd, 1e-12)
  z <- (best - mu) / sd
  (best - mu) * pnorm(z) + sd * dnorm(z)
}

lcb_min <- function(mu, sd, kappa = 2) mu - kappa * sd

#-------------------------
# Single-run RF-SMAC PV
#-------------------------
run_rf_smac_pressure_vessel <- function(
    n_init = 50,
    n_iter = 150,
    n_cand = 120,
    kappa = 2,
    acq = c("EI","LCB"),
    seed = 1,
    num_trees = 500,
    min_node = 5
) {
  set.seed(seed)
  acq <- match.arg(acq)
  
  # init archive
  X_list <- replicate(n_init, random_design_global(), simplify = FALSE)
  X_df <- rbindlist(lapply(X_list, design_to_row), fill = TRUE)
  
  costs <- vapply(seq_len(nrow(X_df)), function(i) {
    pressure_vessel_cost(
      ns = X_df$ns[i], nh = X_df$nh[i],
      R = X_df$R[i], L = X_df$L[i],
      ys = as.character(X_df$ys[i]),
      yh = as.character(X_df$yh[i])
    )
  }, numeric(1))
  
  best_idx <- which.min(costs)
  best_cost <- costs[best_idx]
  best_x <- X_list[[best_idx]]
  
  N_total <- n_init + n_iter
  log <- data.table(step = 1:N_total, best_cost = NA_real_, selected_cost = NA_real_)
  log$selected_cost[1:n_init] <- costs
  log$best_cost[1:n_init] <- cummin(costs)
  
  seen <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(X_df))) seen[[design_key(X_list[[i]])]] <- costs[i]
  
  for (t in 1:n_iter) {
    # fit RF on archive
    X_mm <- make_feature_matrix(X_df)
    rf <- ranger(
      x = as.data.frame(X_mm),
      y = costs,
      num.trees = num_trees,
      min.node.size = min_node,
      seed = seed + t
    )
    
    # candidates
    cand_list <- generate_candidates(best_x, n_local = 90, n_global = 30)
    if (length(cand_list) > n_cand) cand_list <- cand_list[1:n_cand]
    cand_df <- rbindlist(lapply(cand_list, design_to_row), fill = TRUE)
    X_cand_mm <- make_feature_matrix(cand_df)
    
    pr <- rf_predict_mean_sd(rf, as.data.frame(X_cand_mm))
    mu <- pr$mu
    sd <- pr$sd
    
    if (acq == "EI") {
      a <- ei_min(mu, sd, best = best_cost)
      pick <- which.max(a)  # EI: maximize improvement
    } else {
      a <- lcb_min(mu, sd, kappa = kappa)
      pick <- which.min(a)  # LCB: minimize
    }
    
    x_next <- cand_list[[pick]]
    k <- design_key(x_next)
    
    if (!is.null(seen[[k]])) {
      c_next <- seen[[k]]
    } else {
      c_next <- pressure_vessel_cost(x_next$ns, x_next$nh, x_next$R, x_next$L, x_next$ys, x_next$yh)
      seen[[k]] <- c_next
    }
    
    # update archive
    X_list <- c(X_list, list(x_next))
    X_df <- rbind(X_df, as.data.table(design_to_row(x_next)), fill = TRUE)
    costs <- c(costs, c_next)
    
    if (c_next < best_cost) {
      best_cost <- c_next
      best_x <- x_next
    }
    
    idx <- n_init + t
    log$selected_cost[idx] <- c_next
    log$best_cost[idx] <- min(log$best_cost[idx-1], c_next)
  }
  
  list(best_cost = best_cost, best_x = best_x, log = log)
}

#-------------------------
# Multi-run RF-SMAC PV
#-------------------------
run_multirun_rf_pv <- function(
    R = 10,
    seeds = 1:10,
    n_init = 50,
    n_iter = 150,
    n_cand = 120,
    kappa = 2,
    acq = c("EI","LCB")
) {
  acq <- match.arg(acq)
  stopifnot(length(seeds) >= R)
  
  N_total <- n_init + n_iter
  best_mat <- matrix(NA_real_, nrow = R, ncol = N_total)
  finals <- numeric(R)
  
  for (r in 1:R) {
    res <- run_rf_smac_pressure_vessel(
      n_init = n_init,
      n_iter = n_iter,
      n_cand = n_cand,
      kappa = kappa,
      acq = acq,
      seed = seeds[r]
    )
    best_mat[r, ] <- res$log$best_cost
    finals[r] <- tail(res$log$best_cost, 1)
    if (r %% 5 == 0) cat(sprintf("PV RF-%s completed %d/%d runs\n", acq, r, R))
  }
  
  traj <- data.table(
    step = 1:N_total,
    mean_best = colMeans(best_mat),
    sd_best = apply(best_mat, 2, sd)
  )
  
  summary <- data.table(
    method = sprintf("RF-SMAC (%s)", acq),
    mean_final = mean(finals),
    sd_final = sd(finals),
    best_min = min(finals)
  )
  
  list(summary = summary, traj = traj, finals = finals)
}

#-------------------------
# LaTeX table helper for PV
#-------------------------
latex_table_pv <- function(summary_dt, caption, label) {
  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lccc}",
    "\\hline",
    "Method & Mean final best & SD final best & Best (min) \\\\",
    "\\hline"
  )
  for (i in 1:nrow(summary_dt)) {
    lines <- c(lines, sprintf("%s & %.2f & %.2f & %.2f \\\\",
                              summary_dt$method[i], summary_dt$mean_final[i],
                              summary_dt$sd_final[i], summary_dt$best_min[i]))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  cat(paste(lines, collapse = "\n"), "\n")
}

# -------------------------
# Example usage (PV RF baselines, 10 runs)
# -------------------------
pv_rf_ei  <- run_multirun_rf_pv(R=10, seeds=1:10, n_init=50, n_iter=150, acq="EI")
pv_rf_lcb <- run_multirun_rf_pv(R=10, seeds=1:10, n_init=50, n_iter=150, acq="LCB")
print(pv_rf_ei$summary); print(pv_rf_lcb$summary)
latex_table_pv(rbind(pv_rf_ei$summary, pv_rf_lcb$summary),
               caption="Pressure vessel RF-SMAC baseline results over $R=10$ runs (lower is better).",
               label="tab:pv_rf_smac")

# ==========================================================
# PATCH: Combine MNL multi-run + RF-SMAC multi-run (Pressure Vessel)
# ==========================================================
library(data.table)
library(ggplot2)

# -------------------------
# 0) Build pv_out from your existing MNL multi-run objects
#     You already have: logs_all (columns: seed, iter, best_cost, selected_cost)
# -------------------------

# Safety check
stopifnot(exists("logs_all"))
stopifnot(all(c("seed","iter","best_cost") %in% names(logs_all)))

# MNL summary over runs (final best per seed)
pv_mnl_final <- logs_all[, .SD[which.max(iter)], by = seed]
pv_mnl_summary <- data.table(
  method = "MNL-BO (UCB)",
  mean_final = mean(pv_mnl_final$best_cost),
  sd_final   = sd(pv_mnl_final$best_cost),
  best_min   = min(pv_mnl_final$best_cost)
)

# MNL trajectory (mean ± sd over iterations)
pv_mnl_traj <- logs_all[, .(
  mean_best = mean(best_cost),
  sd_best   = sd(best_cost)
), by = iter]
setnames(pv_mnl_traj, "iter", "step")
pv_mnl_traj[, method := "MNL-BO (UCB)"]

# Assemble pv_out in the same “shape” as your TSP out object
pv_out <- list(
  summary = pv_mnl_summary,
  traj    = pv_mnl_traj
)

# -------------------------
# 1) Combine summaries (MNL + RF-SMAC)
#     You already have: pv_rf_ei$summary, pv_rf_lcb$summary
# -------------------------
stopifnot(exists("pv_rf_ei"))
stopifnot(exists("pv_rf_lcb"))

combined_pv <- rbindlist(
  list(as.data.table(pv_out$summary),
       as.data.table(pv_rf_ei$summary),
       as.data.table(pv_rf_lcb$summary)),
  fill = TRUE
)

print(combined_pv)

# -------------------------
# 2) LaTeX table printer (same format as TSP)
# -------------------------
latex_table_pv_combined <- function(summary_dt,
                                    caption = "Pressure vessel multi-run results over $R=10$ runs (lower is better).",
                                    label   = "tab:pv_multirun_all") {
  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lccc}",
    "\\hline",
    "Method & Mean final best & SD final best & Best (min) \\\\",
    "\\hline"
  )
  for (i in 1:nrow(summary_dt)) {
    lines <- c(lines, sprintf("%s & %.2f & %.2f & %.2f \\\\",
                              summary_dt$method[i],
                              summary_dt$mean_final[i],
                              summary_dt$sd_final[i],
                              summary_dt$best_min[i]))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  cat(paste(lines, collapse = "\n"), "\n")
}

latex_table_pv_combined(
  combined_pv,
  caption = "Pressure vessel multi-run results over $R=10$ runs (lower is better).",
  label   = "tab:pv_multirun_all"
)

# -------------------------
# 3) Build combined trajectories for plotting (MNL + RF-SMAC)
# -------------------------
pv_rf_ei_traj  <- copy(as.data.table(pv_rf_ei$traj))
pv_rf_lcb_traj <- copy(as.data.table(pv_rf_lcb$traj))

pv_rf_ei_traj[,  method := "RF-SMAC (EI)"]
pv_rf_lcb_traj[, method := "RF-SMAC (LCB)"]

# Ensure consistent column names
# (your RF traj uses step/mean_best/sd_best already)
stopifnot(all(c("step","mean_best","sd_best","method") %in% names(pv_rf_ei_traj)))
stopifnot(all(c("step","mean_best","sd_best","method") %in% names(pv_rf_lcb_traj)))
stopifnot(all(c("step","mean_best","sd_best","method") %in% names(pv_out$traj)))

pv_traj_all <- rbindlist(list(
  as.data.table(pv_out$traj),
  pv_rf_ei_traj,
  pv_rf_lcb_traj
), fill = TRUE)

pv_traj_all[, ymin := mean_best - sd_best]
pv_traj_all[, ymax := mean_best + sd_best]

# order legend nicely
pv_traj_all[, method := factor(method, levels = c("MNL-BO (UCB)", "RF-SMAC (EI)", "RF-SMAC (LCB)"))]
