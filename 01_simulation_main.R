##############################################################################
# 01_simulation_main.R
# Reproducible code for the oracle and main Monte Carlo analyses reported in
# FinalReport_JiayiWang.pdf.
#
# Report design:
#   Y = 0.5 X1 + 0.8 X2 + 0.3 Z + error
#   X1, X2 ~ N(0,1), Z ~ Bernoulli(0.5)
#   W1 = X1 + N(0,0.8^2), W2 = X2 + N(0,0.8^2)
#   complete-case probability = 0.05 + 0.35*H(0.5 Z)
#   Phase 1 prefers Z=1; Phase 2 prefers Z=0; each ~50%
#   oracle N=50,000
#   Monte Carlo N=1,000 and 2,000; 1,000 replications
#   probability modes: known ("true") and correctly estimated ("correct")
##############################################################################

suppressPackageStartupMessages({
  library(parallel)
})

source("00_mp_blipw_glm.R")

H       <- function(u) 1 / (1 + exp(-u))
SIGMA_E <- 0.8
GAMMA   <- 2.5
BETA    <- c(0, 0.5, 0.8, 0.3)
PARAMS  <- c("beta0", "beta_X1", "beta_X2", "beta_Z")

safe_solve <- function(A, B = NULL, ridge = 1e-8) {
  A <- (A + t(A)) / 2
  A <- A + ridge * diag(nrow(A))
  if (is.null(B)) solve(A) else solve(A, B)
}

make_true_pi_obj <- function(pi_true) {
  pi_true <- pmin(pmax(as.numeric(pi_true), 0.05), 0.999999)
  n <- length(pi_true)
  list(
    pi = pi_true,
    D = matrix(0, n, 1),
    fit = NULL,
    mu = pi_true,
    scores = matrix(0, n, 1),
    I_gam = matrix(1, 1, 1)
  )
}

oracle_comp <- function(N = 50000, seed = 42) {
  set.seed(seed)

  Z  <- rbinom(N, 1, 0.5)
  X1 <- rnorm(N)
  X2 <- rnorm(N)
  Y  <- BETA[1] + BETA[2] * X1 + BETA[3] * X2 + BETA[4] * Z + rnorm(N)
  W1 <- X1 + SIGMA_E * rnorm(N)
  W2 <- X2 + SIGMA_E * rnorm(N)

  pi_cc <- 0.05 + 0.35 * H(0.5 * Z)
  R <- rbinom(N, 1, pi_cc)

  a1 <- uniroot(function(a) mean(H(a + GAMMA * Z)) - 0.50, c(-8, 8))$root
  a2 <- uniroot(function(a) mean(H(a - GAMMA * Z)) - 0.50, c(-8, 8))$root

  q1 <- H(a1 + GAMMA * Z)
  q2 <- H(a2 - GAMMA * Z)

  R1 <- pmax(R, rbinom(N, 1, q1))
  R2 <- pmax(R, rbinom(N, 1, q2))

  pi1 <- ifelse(R == 1, 1, q1)
  pi2 <- ifelse(R == 1, 1, q2)

  ii <- which(R == 1)
  w  <- 1 / pmax(pi_cc[ii], 0.05)
  p  <- 4

  Dp <- cbind(1, X1[ii], X2[ii], Z[ii])
  G1 <- crossprod(Dp, Dp * w) / N
  b_s <- safe_solve(G1 * N, crossprod(Dp, Y[ii] * w))
  ep <- as.vector(Y[ii] - Dp %*% b_s)

  phi_cc <- Dp * (w * ep)
  phi_all <- matrix(0, N, p)
  phi_all[ii, ] <- phi_cc

  G1_inv <- safe_solve(G1)
  V_sipw <- G1_inv %*% (crossprod(phi_cc) / N) %*% t(G1_inv) / N

  do_phase <- function(Wk, Rk, pik) {
    ik <- which(Rk == 1)
    wk <- 1 / pmax(pik[ik], 0.05)

    Dw_cc <- cbind(1, Wk[ii])
    Dw_ph <- cbind(1, Wk[ik])

    G2 <- crossprod(Dw_cc, Dw_cc * w) / N
    th_cc <- safe_solve(G2 * N, crossprod(Dw_cc, Y[ii] * w))

    G3 <- crossprod(Dw_ph, Dw_ph * wk) / N
    th_ph <- safe_solve(G3 * N, crossprod(Dw_ph, Y[ik] * wk))

    psi_cc <- Dw_cc * (w * as.vector(Y[ii] - Dw_cc %*% th_cc))
    psi_ph <- Dw_ph * (wk * as.vector(Y[ik] - Dw_ph %*% th_ph))

    psi_cc_all <- matrix(0, N, 2)
    psi_ph_all <- matrix(0, N, 2)
    psi_cc_all[ii, ] <- psi_cc
    psi_ph_all[ik, ] <- psi_ph

    IC <- psi_cc_all - psi_ph_all
    B  <- G1_inv %*% (crossprod(phi_all, IC) / N)
    A  <- crossprod(IC) / N
    list(B = B, A = A, IC = IC)
  }

  r1 <- do_phase(W1, R1, pi1)
  r2 <- do_phase(W2, R2, pi2)

  V_b1 <- V_sipw - r1$B %*% safe_solve(r1$A) %*% t(r1$B) / N
  V_b2 <- V_sipw - r2$B %*% safe_solve(r2$A) %*% t(r2$B) / N

  A12 <- crossprod(r1$IC, r2$IC) / N
  A_full <- rbind(
    cbind(r1$A, A12),
    cbind(t(A12), r2$A)
  )
  A_full <- (A_full + t(A_full)) / 2
  B_full <- cbind(r1$B, r2$B)
  V_mp <- V_sipw - B_full %*% safe_solve(A_full) %*% t(B_full) / N

  list(
    V_sipw = V_sipw,
    V_b1 = V_b1,
    V_b2 = V_b2,
    V_mp = V_mp,
    miss = 1 - mean(R),
    frac1 = mean(R1),
    frac2 = mean(R2),
    pi1_range = range(q1),
    pi2_range = range(q2)
  )
}

one_rep <- function(n, seed, pi_mode = c("true", "correct")) {
  pi_mode <- match.arg(pi_mode)
  set.seed(seed)

  Z  <- rbinom(n, 1, 0.5)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y  <- BETA[1] + BETA[2] * X1 + BETA[3] * X2 + BETA[4] * Z + rnorm(n)
  W1 <- X1 + SIGMA_E * rnorm(n)
  W2 <- X2 + SIGMA_E * rnorm(n)

  pi_cc <- 0.05 + 0.35 * H(0.5 * Z)
  R <- rbinom(n, 1, pi_cc)
  if (mean(R) < 0.10 || mean(R) > 0.80) return(NULL)

  a1 <- uniroot(function(a) mean(H(a + GAMMA * Z)) - 0.50, c(-8, 8))$root
  a2 <- uniroot(function(a) mean(H(a - GAMMA * Z)) - 0.50, c(-8, 8))$root
  q1 <- H(a1 + GAMMA * Z)
  q2 <- H(a2 - GAMMA * Z)

  set.seed(seed + 1000000L)
  R1 <- pmax(R, rbinom(n, 1, q1))
  set.seed(seed + 2000000L)
  R2 <- pmax(R, rbinom(n, 1, q2))

  if (pi_mode == "true") {
    pi_obj <- make_true_pi_obj(pi_cc)
    pi1 <- ifelse(R == 1, 1, pmax(q1, 0.05))
    pi2 <- ifelse(R == 1, 1, pmax(q2, 0.05))
  } else {
    pi_obj <- tryCatch(
      fit_selection(R = R, D_sel = cbind(1, Z)),
      error = function(e) NULL
    )
    if (is.null(pi_obj)) return(NULL)

    # Avoid fragile model-variable naming by fitting with an explicit data frame.
    fit_phase <- function(Rk) {
      idx0 <- which(R == 0)
      d0 <- data.frame(Rk = Rk[idx0], Z = Z[idx0])
      fit <- tryCatch(glm(Rk ~ Z, data = d0, family = binomial()), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      qhat <- pmax(0.05, pmin(0.999999,
        as.vector(predict(fit, newdata = data.frame(Z = Z), type = "response"))
      ))
      ifelse(R == 1, 1, qhat)
    }

    pi1 <- fit_phase(R1)
    pi2 <- fit_phase(R2)
    if (is.null(pi1) || is.null(pi2)) return(NULL)
  }

  Xmat <- cbind(X1, X2)
  D1 <- cbind(1, W1)
  D2 <- cbind(1, W2)

  sip <- tryCatch(
    sipw_fit_glm(Y, Xmat, Z, R, gaussian(), pi_obj),
    error = function(e) NULL
  )
  bl1 <- tryCatch(
    mp_blipw_glm(Y, Xmat, Z, list(D1), pi_obj, list(pi1), R, list(R1), gaussian()),
    error = function(e) NULL
  )
  bl2 <- tryCatch(
    mp_blipw_glm(Y, Xmat, Z, list(D2), pi_obj, list(pi2), R, list(R2), gaussian()),
    error = function(e) NULL
  )
  mp2 <- tryCatch(
    mp_blipw_glm(Y, Xmat, Z, list(D1, D2), pi_obj, list(pi1, pi2),
                 R, list(R1, R2), gaussian()),
    error = function(e) NULL
  )

  if (any(vapply(list(sip, bl1, bl2, mp2), is.null, logical(1)))) return(NULL)
  list(sip = sip, bl1 = bl1, bl2 = bl2, mp2 = mp2, miss = 1 - mean(R))
}

run_reps <- function(n, pi_mode, n_rep = 1000, seed = 20260419,
                     cores = max(1L, parallel::detectCores() - 1L)) {
  seeds <- seed + seq_len(n_rep)

  worker <- function(s) one_rep(n, s, pi_mode)

  if (.Platform$OS.type == "windows" || cores <= 1L) {
    res <- lapply(seeds, worker)
  } else {
    res <- parallel::mclapply(seeds, worker, mc.cores = cores)
  }

  Filter(Negate(is.null), res)
}

summarise_mc <- function(res) {
  if (length(res) == 0) stop("No valid Monte Carlo replications.")
  z <- qnorm(0.975)
  keys <- c("sip", "bl1", "bl2", "mp2")
  labels <- c("SIPW", "BLIPW (W1)", "BLIPW (W2)", "MP-BLIPW")

  do.call(rbind, lapply(seq_along(BETA), function(j) {
    sip_est <- vapply(res, function(x) x$sip$coef[j], numeric(1))
    sip_var <- var(sip_est)

    do.call(rbind, lapply(seq_along(keys), function(m) {
      est <- vapply(res, function(x) x[[keys[m]]]$coef[j], numeric(1))
      se  <- vapply(res, function(x) x[[keys[m]]]$se[j], numeric(1))
      ok <- is.finite(est) & is.finite(se) & se > 0
      est <- est[ok]; se <- se[ok]

      data.frame(
        parameter = PARAMS[j],
        true = BETA[j],
        method = labels[m],
        bias = mean(est) - BETA[j],
        empirical_sd = sd(est),
        average_se = mean(se),
        coverage = mean(abs(est - BETA[j]) <= z * se),
        rmse = sqrt(mean((est - BETA[j])^2)),
        variance_reduction_pct =
          if (m == 1) 0 else 100 * (1 - var(est) / sip_var),
        n_valid = length(est),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

make_oracle_table <- function(ov) {
  se_sip <- sqrt(diag(ov$V_sipw))
  se_b1  <- sqrt(diag(ov$V_b1))
  se_b2  <- sqrt(diag(ov$V_b2))
  se_mp  <- sqrt(diag(ov$V_mp))

  data.frame(
    parameter = PARAMS,
    true = BETA,
    sipw_se = se_sip,
    blipw_w1_se = se_b1,
    blipw_w2_se = se_b2,
    mp_blipw_se = se_mp,
    blipw_w1_vr_pct = 100 * (1 - se_b1^2 / se_sip^2),
    blipw_w2_vr_pct = 100 * (1 - se_b2^2 / se_sip^2),
    mp_blipw_vr_pct = 100 * (1 - se_mp^2 / se_sip^2)
  )
}

make_figure_1 <- function(tab, filename) {
  png(filename, width = 1500, height = 700, res = 160)
  old <- par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
  on.exit({par(old); dev.off()}, add = TRUE)

  se_mat <- rbind(tab$sipw_se, tab$blipw_w1_se, tab$blipw_w2_se, tab$mp_blipw_se)
  colnames(se_mat) <- tab$parameter
  rownames(se_mat) <- c("SIPW", "BLIPW W1", "BLIPW W2", "MP-BLIPW")
  barplot(se_mat, beside = TRUE, las = 2, ylab = "Standard error",
          main = "Oracle standard errors", legend.text = rownames(se_mat),
          args.legend = list(x = "topright", cex = 0.8))

  vr_mat <- rbind(tab$blipw_w1_vr_pct, tab$blipw_w2_vr_pct, tab$mp_blipw_vr_pct)
  colnames(vr_mat) <- tab$parameter
  rownames(vr_mat) <- c("BLIPW W1", "BLIPW W2", "MP-BLIPW")
  barplot(vr_mat, beside = TRUE, las = 2, ylab = "Variance reduction (%)",
          main = "Variance reduction relative to SIPW",
          legend.text = rownames(vr_mat),
          args.legend = list(x = "topright", cex = 0.8))
}

make_figure_2 <- function(tab, filename) {
  dat <- tab[tab$pi_mode == "correct", ]
  Ns <- sort(unique(dat$N))

  png(filename, width = 1500, height = 1100, res = 160)
  old <- par(mfrow = c(2, 2), mar = c(6, 4, 3, 1))
  on.exit({par(old); dev.off()}, add = TRUE)

  metrics <- list(
    "Variance reduction (%)" = "variance_reduction_pct",
    "Empirical SD" = "empirical_sd",
    "Average estimated SE" = "average_se",
    "Coverage (%)" = "coverage"
  )

  for (nm in names(metrics)) {
    v <- metrics[[nm]]
    sub <- dat[dat$method == "MP-BLIPW", ]
    mat <- sapply(Ns, function(nn) {
      x <- sub[sub$N == nn, ]
      x <- x[match(PARAMS, x$parameter), ]
      val <- x[[v]]
      if (v == "coverage") val <- 100 * val
      val
    })
    rownames(mat) <- PARAMS
    colnames(mat) <- paste0("N=", Ns)
    barplot(t(mat), beside = TRUE, las = 2, ylab = nm, main = nm,
            legend.text = colnames(mat),
            args.legend = list(x = "topright", cex = 0.8))
    if (v == "coverage") abline(h = 95, lty = 2)
  }
}

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

n_rep <- as.integer(Sys.getenv("MPBLIPW_REPS", unset = "1000"))
cores <- as.integer(Sys.getenv(
  "MPBLIPW_CORES",
  unset = as.character(max(1L, parallel::detectCores() - 1L))
))

message("Running oracle analysis...")
oracle <- oracle_comp(N = 50000, seed = 42)
oracle_table <- make_oracle_table(oracle)
write.csv(oracle_table, "results/table_01_oracle_results.csv", row.names = FALSE)
make_figure_1(oracle_table, "figures/figure_01_oracle_efficiency.png")

message("Running Monte Carlo simulations (", n_rep, " replications per scenario)...")
all_res <- list()
all_summary <- list()
idx <- 1L

for (n in c(1000, 2000)) {
  for (pm in c("true", "correct")) {
    message("  N=", n, ", probability mode=", pm)
    rr <- run_reps(n, pm, n_rep = n_rep, cores = cores)
    if (length(rr) == 0) stop("No valid replications for N=", n, ", mode=", pm)
    all_res[[paste(n, pm, sep = "_")]] <- rr
    ss <- summarise_mc(rr)
    ss$N <- n
    ss$pi_mode <- pm
    all_summary[[idx]] <- ss
    idx <- idx + 1L
  }
}

summary_table <- do.call(rbind, all_summary)
write.csv(summary_table, "results/table_02_main_monte_carlo.csv", row.names = FALSE)
saveRDS(list(oracle = oracle, monte_carlo = all_res, summary = summary_table),
        "results/simulation_main_results.rds")
make_figure_2(summary_table, "figures/figure_02_main_monte_carlo.png")

message("Main simulation complete.")
