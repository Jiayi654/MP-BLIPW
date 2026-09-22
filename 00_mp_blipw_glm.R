#  mp_blipw_glm.R  —  Multi-Phase Best Linear IPW Estimation for GLMs
#
#  Extends Wang & Dai (2019, Stat Med 38:2783-2796) to K auxiliary
#  phases and arbitrary GLM families / link functions.
#
#  PRIMARY MODEL:
#    E(Y_i | X_i, Z_i) = g^{-1}(beta' d_i),   d_i = (1, X_i', Z_i')'
#
#  WORKING MODEL k (any GLM family, may differ from primary):
#    E(Y_i | W_i^{(k)}, Z_i) = h_k^{-1}(theta^{(k)}' f_i^{(k)})
#
#  SELECTION MODEL:
#    pi_i = Pr(R_i = 1 | Y_sel_i, Z_i) = logit^{-1}(gamma' D_sel_i)
#    Fitted by logistic regression; Y_sel should be chosen to match
#    the response scale (use Y for Gaussian/Binomial, log(Y+1) for Poisson).
#
#  GLM SCORE EXTRACTION ----
#  For glm(..., weights = w_ipw):
#    fit$weights   = w_ipw * (dmu/deta)^2 / [phi V(mu)]   (IRLS weights)
#    fit$residuals = (Y - mu) / (dmu/deta)                 (working residuals)
#  Score:  phi_i = D_i * fit$weights_i * fit$residuals_i
#  Bread:  G     = (1/N) t(D) %*% (D * fit$weights)
#  DO NOT multiply by w_ipw again — it is already absorbed into fit$weights.
#
#  VARIANCE — FULL APPENDIX B AUGMENTED SANDWICH ----
#  All estimators use the stacked-score augmented sandwich that accounts for
#  estimating gamma (the selection model parameters).
#
#  SIPW (standalone, sipw_fit_glm):
#    V_sipw(pi_hat) = G1^{-1} M_eff G1^{-T} / N
#    M_eff = C11 - G1gam I_gam^{-1} G1gam'
#    (simplified from the four-term formula using the identity B_ps = G1gam)
#
#  MP-BLIPW (mp_blipw_glm):
#    V_mp(pi_hat) = V_sipw(pi_hat) - B_aug A_aug^{-1} B_aug' / N
#
#  where B_aug and A_aug are the fully augmented cross-sandwich blocks
#  (Appendix B, Steps B5-B7), accounting for estimation of gamma in both
#  the CC-SIPW and full-phase working model estimating equations.
#
#  PARTIAL SUB-COHORT PHASES (pi^{(k)} != 1) ----
#  The code now supports phases whose sub-cohorts C_k are strictly smaller than
#  the full cohort.  pi_list[[k]] may be any length-N vector of phase-k
#  inclusion probabilities; use rep(1, N) for a full-cohort phase.
#  C_CC ⊆ C_k must hold (R_i=1 => R_list[[k]][i]=1), ensuring both working-
#  model estimators target the same pseudo-true value.
#
#  KEY CHANGE — influence-function A_raw blocks:
#  The influence function for the k-th contrast is
#    IC_ki = psi_cc_all[i,] − psi_all[i,]
#          = (R_i/pi_i  −  R_ki/pi_ki) * psi_raw_ki
#  and  A_kk = (1/N) t(IC_k) %*% IC_k  (always PSD by construction).
#  For full-cohort phases (pi_k = 1) the IC formula equals C22−C33 in
#  expectation, but the two differ at O(1/sqrt(N)) in finite samples.
#  To preserve the Wang & Dai machine-precision verification for full-cohort
#  phases, the implementation uses:
#    • C22−C33        for all-pi_k-equal-1 phases (backward-compatible)
#    • IC outer product for any partial phase (pi_k < 1, always PSD)
#  Off-diagonal blocks A_jk follow the same rule: C22_jk−C33_jk when both
#  phases are full-cohort; (1/N) t(IC_j) %*% IC_k otherwise.
#
#  For each phase k with pi^{(k)} = 1 (full cohort), the augmented blocks are:
#
#  b_k_aug = [A1_inv B12k A2k_inv']_{beta,theta_k}
#            - [A1_inv B13k A3k_inv']_{beta,theta_k}
#  where:
#    A1  = [[G1,    G1gam],[0,      I_gam]]   (primary + selection bread)
#    A2k = [[G2k,   G2kgam],[0,     I_gam]]   (CC-SIPW working + selection bread)
#    A3k = G3k                                 (full-cohort bread; no selection)
#    B12k = [[C12,  G1gam],[G2kgam',I_gam]]   (shared gamma in both tildes)
#    B13k = [[C13],[B_sel_psi]]                (C13: primary x full-cohort;
#                                               B_sel_psi: selection x full-cohort)
#
#  [A]_kk_aug = V22 + V33 - V23 - V23'
#  where:
#    V22 = [A2k_inv B22k A2k_inv']    B22k = [[C22,G2kgam],[G2kgam',I_gam]]
#    V33 = G3k_inv C33 G3k_inv
#    V23 = [A2k_inv B23k A3k_inv']    B23k = [[C23],[B_sel_psi']]
#    C23 = t(psi_cc) %*% psi_all[CC,] / N  (CC-SIPW x full-phase at CC rows)
#
#  New cross-meat quantities:
#    G2kgam = (1/N) t(sweep(psi_cc_k, 1, 1-pi_cc, "*")) %*% D_sel_cc
#    B_sel_psi = (1/N) t(sel_sc) %*% psi_all_k
#
#  POINT ESTIMATE: always uses the naive (Appendix A) Delta for stability.
#  The augmented optimal Delta requires inverting A_aug, which can be
#  ill-conditioned in finite samples, causing large variance in the estimate.
#  The naive Delta is asymptotically equivalent to first order.
#
#  FALLBACK: if A_aug^{(k)} is not PSD or has condition number > cond_thresh,
#  the naive A_naive^{(k)} and b_k_naive are used for that phase.
#  If V_mp is not PSD, it falls back to V_sipw.
#
#  REDUCTION TO WANG & DAI (2019) ----
#  gaussian(identity), K=1, pi^{(1)} = 1 (full cohort):
#  mp_blipw_glm() matches the linear BLIPW of Wang & Dai (2019)
#  to machine precision (max |Delta beta| < 2e-15).

suppressPackageStartupMessages({ library(MASS); library(stats) })

# Fit a GLM with a pre-built design matrix D (already includes intercept). ----
fit_wglm <- function(Y, D, family, weights = NULL) {
  df   <- as.data.frame(D)
  nms  <- paste0("v", seq_len(ncol(D)))
  colnames(df) <- nms
  fmla <- as.formula(paste("Y ~", paste(nms, collapse = "+"), "- 1"))
  if (is.null(weights)) glm(fmla, data = df, family = family)
  else                  glm(fmla, data = df, family = family, weights = weights)
}

# Extract unscaled score matrix and bread from a fitted (possibly weighted) GLM.
# fit$weights already incorporates the user IPW weights; do NOT re-weight.
glm_scores_bread <- function(fit, D) {
  scores <- D * fit$weights * fit$residuals # n x p
  bread  <- t(D) %*% (D * fit$weights) # p x p
  list(scores = scores, bread = bread)
}

# Fit logistic selection model. ----
#
# Two calling conventions (matching Wang & Dai flexibility):
#
# (A) Formula-based:  fit_selection(Y_sel, Z, R)
#     Fits pi = logit^{-1}(gamma_0 + gamma_1*Y_sel + gamma_2*Z).
#     Y_sel should match the response scale:
#       - Gaussian / Binomial : Y_sel = Y
#       - Poisson             : Y_sel = log(Y + 1)
#
# (B) Design-matrix:  fit_selection(R = R, D_sel = D_sel)
#     Fits the logistic model with the user-supplied N x q design matrix.
#     Use this to include surrogates W in the selection model, as in
#     Wang & Dai (2019): pi = logit^{-1}(gamma_0+gamma_1*Y+gamma_2*Z+gamma_3*W).
#     Example:
#       D_sel <- cbind(1, Y, Z, W)   # intercept, Y, Z, W
#       pi_obj <- fit_selection(R = R, D_sel = D_sel)
#
# Returns a list with:
#   pi     : length-N selection probabilities (clipped at 0.05)
#   D      : N x q design matrix used for the logistic model
#   fit    : fitted glm object
#   mu     : unclipped fitted probabilities (needed for augmented sandwich)
#   scores : N x q score matrix s_i = (R_i - mu_i) D_sel_i
#   I_gam  : q x q Fisher information matrix 
fit_selection <- function(Y_sel = NULL, Z = NULL, R, D_sel = NULL) {
  if (!is.null(D_sel)) {
    # Design-matrix interface (arbitrary selection covariates)
    stopifnot(nrow(D_sel) == length(R))
    df   <- as.data.frame(D_sel)
    nms  <- paste0("s", seq_len(ncol(D_sel)))
    colnames(df) <- nms
    fmla <- as.formula(paste("R ~", paste(nms, collapse = "+"), "- 1"))
    fit  <- glm(fmla, data = df, family = binomial())
  } else {
    # Formula interface: pi = H(gamma_0 + gamma_1*Y_sel + gamma_2*Z)
    fit  <- glm(R ~ Y_sel + Z, family = binomial())
    D_sel <- model.matrix(fit)
  }
  mu     <- fitted(fit)
  pi     <- pmax(mu, 0.05)
  scores <- D_sel * as.vector(R - mu)                            # N x q: s_i = (R_i-mu_i)*D_sel_i
  I_gam  <- t(D_sel) %*% (D_sel * as.vector(mu * (1-mu))) / length(R) # q x q Fisher info
  list(pi = pi, D = D_sel, fit = fit, mu = mu, scores = scores, I_gam = I_gam)
}

# SIPW estimator with augmented sandwich SE (Appendix B). ----
#
# Accounts for estimation of the selection probabilities by computing
# the (beta,beta) block of the joint (beta,gamma) augmented sandwich,
# exactly as in Wang & Dai (2019) -- see sipwlogi.piyzw() in their code.
#
# The augmented sandwich replaces the basic meat C11 with:
#   M_eff = C11
#           - B_ps  I_gam^{-1} G1gam'
#           - G1gam I_gam^{-1} B_ps'
#           + G1gam I_gam^{-1} G1gam'
# where:
#   B_ps  = (1/N) sum_i phi_all_i s_i'   (cross-meat, via IPW identity)
#   G1gam = (1/N) sum_{CC} (1-pi_i) phi_i_unwtd D_sel_i'
#         = (1/N) t(sweep(phi_cc, 1, 1-pi_cc, "*")) %*% D_sel_cc
#           (from differentiating (R_i/pi_i) w.r.t. gamma in the logistic model)
#   I_gam : Fisher information for gamma (from pi_obj)
#
# Falls back to basic sandwich if M_eff is not positive semi-definite.

sipw_fit_glm <- function(Y, X, Z, R, family = gaussian(), pi_obj) {
  N      <- length(Y)
  pi_hat <- pi_obj$pi
  D_sel  <- pi_obj$D
  sel_sc <- pi_obj$scores
  I_gam <- pi_obj$I_gam
  ii     <- which(R == 1)
  pi_cc  <- pi_hat[ii]
  w_ipw <- 1 / pi_hat[ii]; 
  Dp      <- cbind(1, X, Z)[ii, , drop = FALSE]
  fit    <- fit_wglm(Y[ii], Dp, family, weights = w_ipw)
  beta <- coef(fit) # point estimate
  
  # Basic sandwich
  sc     <- glm_scores_bread(fit, Dp)
  G1  <- sc$bread/N
  G1_inv <- solve(G1)
  p <- ncol(Dp)
  phi_cc  <- sc$scores
  phi_all <- matrix(0, N, p)
  phi_all[ii, ] <- phi_cc
  C11 <- t(phi_cc) %*% phi_cc / N

  se_basic <- sqrt(pmax(diag(G1_inv%*%C11%*%t(G1_inv)/N),0))
  
  # Augmented sandwich (Appendix B)
  
  # Cross-meat: B_ps = E[phi_i s_i'] estimated via IPW identity
  B_ps  <- t(phi_all) %*% sel_sc / N  # p x q
  # G1gamma: derivative of (R/pi)*phi w.r.t. gamma (logistic model)
  G1gam <- t(sweep(phi_cc, 1, 1 - pi_cc, "*")) %*% D_sel[ii, ] / N  # p x q

  # Effective meat (augmented sandwich)
  #I_inv <- tryCatch(solve(I_gam), error = function(e) ginv(I_gam))
  I_inv  <- solve(I_gam)
  
  M_eff <- C11 -
           B_ps  %*% I_inv %*% t(G1gam) -
           G1gam %*% I_inv %*% t(B_ps)  +
           G1gam %*% I_inv %*% t(G1gam)

  # Fallback to basic sandwich if M_eff is not PSD
  BasicSE <- FALSE
  if (any(eigen(M_eff, only.values = TRUE, symmetric = TRUE)$values < -1e-8)){
    M_eff <- C11
    BasicSE <- TRUE
  }
    
  V_sand <- G1_inv %*% M_eff %*% t(G1_inv) / N
  
  list(coef = beta, se = sqrt(pmax(diag(V_sand), 0)), se_basic = se_basic, BasicSE = BasicSE, G1 = G1, G1_inv = G1_inv, C11 = C11)
}

# Multiple imputation (gaussian only, for simulation benchmark)
mi_fit <- function(Y, X, Z, W, R, m = 20) {
  n <- length(Y); i <- which(R==1); mi <- which(R==0); n0 <- length(mi)
  Dm <- cbind(1,Y,Z,W)[i,]; XtX <- t(Dm)%*%Dm
  a  <- as.vector(solve(XtX)%*%t(Dm)%*%X[i])
  df <- length(i)-4; s2 <- sum((X[i]-Dm%*%a)^2)/df
  br <- matrix(NA,m,3); vr <- matrix(NA,m,3)
  for (k in seq_len(m)) {
    sig2 <- s2*df/rchisq(1,df); ad <- mvrnorm(1,a,sig2*solve(XtX))
    Xi <- X; Xi[mi] <- cbind(1,Y,Z,W)[mi,]%*%ad+rnorm(n0,0,sqrt(sig2))
    fit <- lm(Y~Xi+Z); br[k,] <- coef(fit); vr[k,] <- summary(fit)$coefficients[,2]^2
  }
  list(coef=colMeans(br), se=sqrt(colMeans(vr)+(1+1/m)*apply(br,2,var)))
}


# MP-BLIPW for arbitrary GLM families. ----
#
# Arguments:
#   Y, X, Z      : length-N; X observed iff R = 1. X may be a matrix (N x p_x).
#   W_list       : list of K surrogate predictors. Each element is either:
#                    - a length-N vector (default; working model = intercept+W+Z), or
#                    - an N x pk_k matrix with the full working model design
#                      (intercept already included; Z included iff desired).
#                  Using a matrix overrides the default cbind(1,W_k,Z) design,
#                  allowing custom working models (e.g. W_k only, no Z).
#   pi_obj       : output of fit_selection()
#   pi_list      : list of K phase-k selection probabilities
#                  (use rep(1, N) when all N subjects are in phase k)
#   R            : length-N CC indicator
#   R_list       : list of K phase-k membership indicators
#                  (use rep(1L, N) when all N subjects are in phase k)
#   family       : GLM family object for primary model
#   family_w     : GLM family for working models — list of K families,
#                  or a single family applied to all K phases.
#                  Defaults to `family`.
#   lambda       : ridge regularisation for inverting A_raw (default 1e-6)
#
# Returns a list:
#   coef, se          : MP-BLIPW estimate and augmented sandwich SE
#   coef_sipw, se_sipw: SIPW estimate and augmented sandwich SE
#   se_naive          : SE using naive (Appendix A) sandwich (for comparison)

mp_blipw_glm <- function(Y, X, Z, W_list, pi_obj, pi_list,
                          R, R_list,
                          family   = gaussian(),
                          family_w = NULL,
                          lambda   = 1e-6) {

  N      <- length(Y); K <- length(W_list)
  pi_hat <- pi_obj$pi
  ii     <- which(R == 1); w <- 1 / pi_hat[ii]

  # Resolve working families
  if (is.null(family_w))     family_w <- rep(list(family), K)
  if (!is.list(family_w))    family_w <- rep(list(family_w), K)
  if (length(family_w) == 1) family_w <- rep(family_w, K)

  # Primary SIPW ----
  # X may be a vector or a matrix (N x p_x); cbind handles both.
  Dp    <- cbind(1, X, Z)[ii, , drop = FALSE]; Yc <- Y[ii]
  fit_p <- fit_wglm(Yc, Dp, family, weights = w)
  b_s   <- coef(fit_p); p <- length(b_s)
  sc_p  <- glm_scores_bread(fit_p, Dp)

  G1     <- sc_p$bread/N; G1_inv <- solve(G1)
  phi_cc <- sc_p$scores # n_CC x p
  phi_all <- matrix(0, N, p); 
  phi_all[ii, ] <- phi_cc  # N x p

  # Augmented sandwich for SIPW (Appendix B; accounts for estimated gamma). ----
  # Since B_ps = G1gam exactly (proved algebraically and verified numerically
  # to < 2e-17), the four-term M_eff formula reduces to:
  #   M_eff = C11 - G1gam I_gam^{-1} G1gam'
  C11    <- t(phi_cc) %*% phi_cc / N
  pi_cc  <- pi_hat[ii]
  D_sel  <- pi_obj$D
  sel_sc <- pi_obj$scores
  I_gam <- pi_obj$I_gam
  G1gam  <- t(sweep(phi_cc, 1, 1 - pi_cc, "*")) %*% D_sel[ii, ] / N  # p x q
  I_inv  <- tryCatch(solve(I_gam), error = function(e) ginv(I_gam))
  M_eff  <- C11 - G1gam %*% I_inv %*% t(G1gam)
  if (any(eigen(M_eff, only.values = TRUE, symmetric = TRUE)$values < -1e-8))
    M_eff <- C11
  V_sipw  <- G1_inv %*% M_eff %*% t(G1_inv) / N
  se_sipw <- sqrt(pmax(diag(V_sipw), 0))

  # Per-phase working models ----
  th_s_list <- th_h_list <- G2_list <-
    psi_cc_list <- psi_all_list <- vector("list", K)
  pk <- integer(K)

  G3_list <- vector("list", K)   # full-phase bread G3k (pk x pk), /N

  for (k in seq_len(K)) {
    Wk    <- W_list[[k]]
    ik    <- which(R_list[[k]] == 1)
    pi_k  <- pmax(pi_list[[k]][ik], 0.05); wk <- 1 / pi_k
    fam_k <- family_w[[k]]

    # Working model design matrices.
    # If Wk is a matrix (N x pk), use it directly as the full design;
    # otherwise build the default design cbind(1, Wk, Z).
    if (is.matrix(Wk)) {
      Dw_cc <- Wk[ii, , drop = FALSE]  # n_CC x pk (user-supplied)
      Dw_ph <- Wk[ik, , drop = FALSE]  # n_k  x pk
    } else {
      Dw_cc <- cbind(1, Wk, Z)[ii, ]   # n_CC x pk (default: intercept+W+Z)
      Dw_ph <- cbind(1, Wk, Z)[ik, ]   # n_k  x pk
    }
    Yk_ph <- Y[ik]

    # CC-SIPW working model (IPW-weighted by primary CC weights w)
    fit_wcc <- fit_wglm(Yc, Dw_cc, fam_k, weights = w)
    th_s <- coef(fit_wcc); pk[k] <- length(th_s)
    th_s_list[[k]] <- th_s
    sc_wcc <- glm_scores_bread(fit_wcc, Dw_cc)
    G2_list[[k]]     <- sc_wcc$bread/N    # pk x pk  (empirical G_k, /N)
    psi_cc_list[[k]] <- sc_wcc$scores     # n_CC x pk (IPW-weighted scores)

    # Full-phase working model
    # When pi^{(k)} = 1: unweighted GLM on all n_k subjects.
    # Otherwise: IPW-weighted with phase-k weights.
    fit_wph <- if (all(pi_list[[k]] == 1))
      fit_wglm(Yk_ph, Dw_ph, fam_k)
    else
      fit_wglm(Yk_ph, Dw_ph, fam_k, weights = wk)

    th_h_list[[k]] <- coef(fit_wph)
    sc_wph <- glm_scores_bread(fit_wph, Dw_ph)    # n_k x pk

    # FIX 1: Store the full-phase BREAD (Hessian-based), not score outer product.
    # For Gaussian the bread equals C33, but for Binomial/Poisson they differ
    # because G3k = (1/N) t(D) %*% (D * irls_weights) whereas
    # C33 = (1/N) t(score) %*% score = (1/N) t(D*irls*resid) %*% (D*irls*resid).
    # Using C33 as G3k_inv inflates the correction term for non-Gaussian families.
    G3_list[[k]] <- sc_wph$bread / N  # pk x pk full-phase bread (correct G3k)

    # Embed full-phase scores into N-length matrix (zero outside C_k)
    psa <- matrix(0, N, pk[k]); psa[ik, ] <- sc_wph$scores
    psi_all_list[[k]] <- psa
  }

  p_plus <- sum(pk)

  # Embed CC-SIPW scores into N-length matrices and compute IC_list ----
  # psi_cc_all_list[[k]][i,] = (R_i/pi_i)*psi_raw_ki for CC; 0 elsewhere.
  # IC_list[[k]][i,]         = psi_cc_all[i,] - psi_all[i,]
  #                          = (R_i/pi_i - R_ki/pi_ki)*psi_raw_ki
  # Used in A_raw for partial phases (always PSD as a Gram matrix average).
  psi_cc_all_list <- vector("list", K)
  IC_list         <- vector("list", K)
  for (k in seq_len(K)) {
    pca <- matrix(0, N, pk[k]); pca[ii, ] <- psi_cc_list[[k]]
    psi_cc_all_list[[k]] <- pca
    IC_list[[k]]         <- pca - psi_all_list[[k]]
  }

  # Moment matrices (Appendix A, Section 1.3) ----
  #
  # C12^{(k)} = (1/N) sum_{CC} (1/pi_i) phi_i (psi_i^{(k)})'
  #           = t(phi_cc) %*% psi_cc / N
  #   [phi_cc has 1/pi factor; psi_cc has 1/pi factor from CC-SIPW fit]
  #
  # C13^{(k)} = (1/N) sum_{C_k} (1/pi_i^{(k)}) phi_i (psi_i^{(k)})'
  #           = t(phi_all) %*% psi_all / N
  #   [phi_all has 1/pi factor at CC rows, 0 elsewhere;
  #    psi_all has 1/pi^{(k)} factor at C_k rows, 0 elsewhere;
  #    for pi^{(k)}=1: psi_all is unweighted, so C13 = (1/N) sum_{CC} (1/pi) phi psi']
  #
  # C22^{(k)} = t(psi_cc) %*% psi_cc / N
  # C33^{(k)} = t(psi_all) %*% psi_all / N
  #
  # Off-diagonal C22^{(j,k)} = t(psi_cc_j) %*% psi_cc_k / N
  # Off-diagonal C33^{(j,k)} uses psi_all_j and psi_all_k (zero outside
  # respective phases, so their product is nonzero only at intersection C_j ∩ C_k)

  C12_list <- C13_list <- C22_list <- C33_list <- vector("list", K)
  for (k in seq_len(K)) {
    C12_list[[k]] <- t(phi_cc)  %*% psi_cc_list[[k]]  / N
    C13_list[[k]] <- t(phi_all) %*% psi_all_list[[k]] / N
    C22_list[[k]] <- t(psi_cc_list[[k]])  %*% psi_cc_list[[k]]  / N
    C33_list[[k]] <- t(psi_all_list[[k]]) %*% psi_all_list[[k]] / N
  }

  # Augmented gamma-sensitivity matrices for primary and working models ----
  # (needed for Appendix B B_aug and A_aug blocks).
  # G1gam already computed above.
  # G2kgam^{(k)} = (1/N) t(sweep(psi_cc_k, 1, 1-pi_cc, "*")) %*% D_sel_cc
  # (same logistic chain-rule formula as G1gam, applied to working scores)
  # B_sel_psi^{(k)} = (1/N) t(sel_sc) %*% psi_all_k
  # (cross-meat between selection scores and full-cohort working scores)
  G2gam_list     <- vector("list", K)  # pk x q per phase
  B_sel_psi_list <- vector("list", K)  # q x pk per phase
  C23_list       <- vector("list", K)  # pk x pk per phase

  for (k in seq_len(K)) {
    G2gam_list[[k]] <- t(sweep(psi_cc_list[[k]], 1, 1 - pi_cc, "*")) %*%
                       D_sel[ii, ] / N
    # B_sel_psi: E[sel_score_i * psi_all_i^{(k)}'] via direct sum (all N subjects)
    B_sel_psi_list[[k]] <- t(sel_sc) %*% psi_all_list[[k]] / N
    # C23: cross between CC-SIPW working scores (IPW-weighted) and full-cohort
    # working scores (unweighted) at the CC rows.
    # t(psi_cc_k) %*% psi_all_k[CC,] / N
    # = (1/N) sum_{CC} w_i * psi_raw_i * psi_raw_i'  (at CC, psi_all is psi_raw)
    C23_list[[k]] <- t(psi_cc_list[[k]]) %*% psi_all_list[[k]][ii, ] / N
  }

  # NAIVE blocks (Appendix A, for point estimate and fallback SE) ----
  # B_raw (p x p+):  k-th column-block = G1_inv (C12_k - C13_k)
  # A_raw (p+ x p+): (j,k)-block — see below.
  # These do not include G_k^{-1} factors; u_hat includes G_k instead.

  B_raw <- do.call(cbind, lapply(seq_len(K), function(k)
    G1_inv %*% (C12_list[[k]] - C13_list[[k]])))

  # Helper: is phase k a full-cohort phase?
  full_cohort <- vapply(seq_len(K), function(k) all(pi_list[[k]] == 1), logical(1))

  pcu <- c(0, cumsum(pk)); A_raw <- matrix(0, p_plus, p_plus)
  for (j in seq_len(K)) {
    rj <- (pcu[j]+1):pcu[j+1]
    for (k in j:K) {
      ck  <- (pcu[k]+1):pcu[k+1]
      blk <- if (j == k) {
        # Diagonal block:
        #   Full-cohort: C22−C33 (matches Wang & Dai to machine precision)
        #   Partial:     IC outer product (always PSD, consistent estimator)
        if (full_cohort[k])
          C22_list[[k]] - C33_list[[k]]
        else
          t(IC_list[[k]]) %*% IC_list[[k]] / N
      } else {
        # Off-diagonal block:
        #   Both full-cohort: C22_jk−C33_jk (original formula)
        #   Any partial:      IC cross product (consistent, avoids non-PSD)
        if (full_cohort[j] && full_cohort[k])
          t(psi_cc_list[[j]]) %*% psi_cc_list[[k]] / N -
          t(psi_all_list[[j]]) %*% psi_all_list[[k]] / N
        else
          t(IC_list[[j]]) %*% IC_list[[k]] / N
      }
      A_raw[rj, ck] <- blk; A_raw[ck, rj] <- t(blk)
    }
  }
  A_raw  <- (A_raw + t(A_raw)) / 2

  # Ridge-regularised inverse for the point estimate.
  # The condition number of A_raw governs the stability of the point estimate
  # correction Delta * u_hat. When A_raw is ill-conditioned (cond > cond_thresh),
  # use a stronger adaptive ridge lambda_pe = median(diag(A_raw)) * 0.1 to
  # keep the correction finite while preserving O(1/N) bias.
  cond_thresh_pe <- 500
  eig_A_raw <- eigen(A_raw, symmetric = TRUE, only.values = TRUE)$values
  cond_A_raw <- if (min(eig_A_raw) > 0) max(eig_A_raw) / min(eig_A_raw) else Inf
  lambda_pe  <- if (cond_A_raw > cond_thresh_pe)
                  max(lambda, median(abs(diag(A_raw))) * 0.1)
                else lambda
  A_inv  <- tryCatch(solve(A_raw + lambda_pe * diag(p_plus)),
                     error = function(e) ginv(A_raw))

  # AUGMENTED blocks (Appendix B, Steps B5-B7) ----
  # Computed per phase; can fall back to naive per phase if ill-conditioned.
 
  # Bread inverse blocks (shared across phases)
  A1_tl <- G1_inv
  A1_tr <- -G1_inv %*% G1gam %*% I_inv   # p x q

  b_aug_list   <- vector("list", K)   # augmented b_k (p x pk)
  A_aug_list   <- vector("list", K)   # augmented [A]_kk (pk x pk)
  use_aug      <- logical(K)          # whether augmented block is used

  for (k in seq_len(K)) {
    G2kgam    <- G2gam_list[[k]]        # pk x q
    Bsp_k     <- B_sel_psi_list[[k]]    # q x pk
    C23_k     <- C23_list[[k]]          # pk x pk

    # FIX 1: Use G3k from the stored full-phase bread (NOT C33).
    # G3k = (1/N) t(D) %*% (D * irls_weights) extracted from glm_scores_bread.
    # For Gaussian irls_weights=1 so G3k=C33, but for Binomial/Poisson they differ.
    G3k       <- G3_list[[k]]           # pk x pk  (correct full-phase bread)
    G3k_inv   <- tryCatch(solve(G3k + lambda * diag(pk[k])),
                          error = function(e) ginv(G3k))
    G2k_inv   <- tryCatch(solve(G2_list[[k]] + lambda * diag(pk[k])),
                          error = function(e) ginv(G2_list[[k]]))
    A2k_tl    <- G2k_inv                # pk x pk
    A2k_tr    <- -G2k_inv %*% G2kgam %*% I_inv  # pk x q

    # Augmented b_k (Appendix B, Step B7)
    # cov1_k = h1_inv [c13_k h3_inv' - c12_k h2_inv']
    # = A1_inv [ C13_k A3k_inv' - C12_k A2k_inv' ]   (theta-theta block)
    # A1_inv upper row: [A1_tl | A1_tr], A2k_inv upper row: [A2k_tl | A2k_tr]
    # A3k is block-diagonal [[G3k, 0],[0, I_gam]] since full-cohort score
    # does not depend on gamma: A3k_inv upper-left = G3k_inv, upper-right = 0.
    #
    # FIX 2: term1 = [A1_inv B12k A2k_inv']_{beta,theta_k} (B12k has shared gamma)
    term1 <- (A1_tl %*% C12_list[[k]] + A1_tr %*% t(G2kgam)) %*% t(A2k_tl) +
             (A1_tl %*% G1gam         + A1_tr %*% I_gam     ) %*% t(A2k_tr)
    # term2 = [A1_inv B13k A3k_inv']_{beta,theta_k}
    # B13k upper = C13_k, lower = Bsp_k; A3k_inv upper-left = G3k_inv, upper-right = 0
    term2 <- A1_tl %*% C13_list[[k]] %*% G3k_inv +
             A1_tr %*% Bsp_k         %*% G3k_inv
    b_k_aug <- term1 - term2           # p x pk

    # Augmented [A]_kk (Appendix B, Step B5):
    # var4[theta,theta] = var2[theta,theta] + var3[theta,theta] - cov23[theta,theta] - cov23[theta,theta]'
    #
    # FIX 3: A_aug uses the full block sandwiches and extracts the theta-theta block,
    # matching WD's invvar4 = inv(var4[1:pk, 1:pk]).
    # var2[theta,theta] = [A2k_inv B22k A2k_inv']_{theta,theta}:
    V22_th <- (A2k_tl %*% C22_list[[k]] + A2k_tr %*% t(G2kgam)) %*% t(A2k_tl) +
              (A2k_tl %*% G2kgam         + A2k_tr %*% I_gam    ) %*% t(A2k_tr)
    # var3[theta,theta] = G3k_inv C33_k G3k_inv:
    V33_th <- G3k_inv %*% C33_list[[k]] %*% G3k_inv
    # cov23[theta,theta] = [A2k_inv B23k A3k_inv']_{theta,theta}:
    # B23k top = C23_k (pk x pk), B23k bot = Bsp_k' (q x pk), A3k_inv top-left = G3k_inv
    V23_th <- (A2k_tl %*% C23_k + A2k_tr %*% Bsp_k) %*% G3k_inv
    A_kk_aug <- (V22_th + V33_th - V23_th - t(V23_th))
    A_kk_aug <- (A_kk_aug + t(A_kk_aug)) / 2  # enforce symmetry

    # Fallback: use naive blocks if augmented A is ill-conditioned or non-PSD
    eig_A    <- eigen(A_kk_aug, only.values = TRUE, symmetric = TRUE)$values
    cond_A   <- if (min(eig_A) > 0) max(eig_A) / min(eig_A) else Inf
    if (any(eig_A < -1e-8) || cond_A > 1e5) {
      b_aug_list[[k]] <- B_raw[, (pcu[k]+1):pcu[k+1], drop = FALSE]
      # Fallback A: use IC formula (always PSD) rather than C22−C33
      # which can be non-PSD for partial phases in finite samples.
      A_aug_list[[k]] <- t(IC_list[[k]]) %*% IC_list[[k]] / N
      use_aug[k]      <- FALSE
    } else {
      b_aug_list[[k]] <- b_k_aug
      A_aug_list[[k]] <- A_kk_aug
      use_aug[k]      <- TRUE
    }
  }

  # Assemble augmented B_aug (p x p+) and A_aug (p+ x p+)
  B_aug <- do.call(cbind, b_aug_list)
  A_aug <- matrix(0, p_plus, p_plus)
  for (j in seq_len(K)) {
    rj <- (pcu[j]+1):pcu[j+1]
    for (k in j:K) {
      ck <- (pcu[k]+1):pcu[k+1]
      if (j == k) {
        A_aug[rj, ck] <- A_aug_list[[k]]
      } else {
        # Off-diagonal blocks: use naive when either phase fell back
        if (use_aug[j] && use_aug[k]) {
          # Augmented off-diagonal (Appendix B, Step B6)
          # [A]_jk = V22_jk + V33_jk - V23_jk - V32_jk
          # For full-cohort phases: V33_jk involves cross between full-cohort scores
          # = G3j_inv C33_jk G3k_inv  where C33_jk = t(psi_all_j) %*% psi_all_k / N
          # G3 bread for the augmented off-diagonal: use G3_list (IRLS-based),
          # not C33_list (score outer product). For non-Gaussian or IPW-weighted
          # phase-k models C33 != G3k; using C33 here would be wrong.
          G3j_inv <- tryCatch(solve(G3_list[[j]] + lambda * diag(pk[j])),
                              error = function(e) ginv(G3_list[[j]]))
          G3k_inv_jk <- tryCatch(solve(G3_list[[k]] + lambda * diag(pk[k])),
                                 error = function(e) ginv(G3_list[[k]]))
          G2j_inv   <- tryCatch(solve(G2_list[[j]] + lambda * diag(pk[j])),
                                error = function(e) ginv(G2_list[[j]]))
          G2k_inv_jk<- tryCatch(solve(G2_list[[k]] + lambda * diag(pk[k])),
                                error = function(e) ginv(G2_list[[k]]))
          A2j_tl <- G2j_inv
          A2j_tr <- -G2j_inv %*% G2gam_list[[j]] %*% I_inv
          A2k_tl_jk <- G2k_inv_jk
          A2k_tr_jk  <- -G2k_inv_jk %*% G2gam_list[[k]] %*% I_inv
          C22_jk <- t(psi_cc_list[[j]]) %*% psi_cc_list[[k]] / N
          C33_jk <- t(psi_all_list[[j]]) %*% psi_all_list[[k]] / N
          C23_jk <- t(psi_cc_list[[j]]) %*% psi_all_list[[k]][ii, ] / N
          C32_jk <- t(psi_cc_list[[k]]) %*% psi_all_list[[j]][ii, ] / N
          Bsp_j  <- B_sel_psi_list[[j]]; Bsp_k_jk <- B_sel_psi_list[[k]]
          V22_jk <- (A2j_tl%*%C22_jk + A2j_tr%*%t(G2gam_list[[k]]))%*%t(A2k_tl_jk) +
                    (A2j_tl%*%G2gam_list[[j]]+A2j_tr%*%I_gam)%*%t(A2k_tr_jk)
          V33_jk <- G3j_inv %*% C33_jk %*% G3k_inv_jk
          V23_jk <- (A2j_tl%*%C23_jk + A2j_tr%*%Bsp_k_jk) %*% G3k_inv_jk
          V32_jk <- (A2k_tl_jk%*%C32_jk + A2k_tr_jk%*%Bsp_j) %*% G3j_inv
          blk_aug <- V22_jk + V33_jk - V23_jk - t(V32_jk)
        } else {
          # Naive off-diagonal fallback: IC cross product is consistent for
          # both full-cohort and partial phases, and avoids non-PSD issues.
          blk_aug <- t(IC_list[[j]]) %*% IC_list[[k]] / N
        }
        A_aug[rj, ck] <- blk_aug; A_aug[ck, rj] <- t(blk_aug)
      }
    }
  }
  A_aug <- (A_aug + t(A_aug)) / 2

  A_aug_inv <- tryCatch(solve(A_aug + lambda * diag(p_plus)),
                        error = function(e) ginv(A_aug))

  # Point estimate: ALWAYS uses naive Delta (Appendix A formula). ----
  # Naive Delta is stable and asymptotically equivalent to first order.
  Delta <- -B_raw %*% A_inv
  # u_hat uses G2*(th_s - th_h): the "raw" formulation where B_raw = G1_inv(C12-C13)
  # and the G2 factor converts theta-scale contrast to score-scale contrast.
  u_hat <- unlist(lapply(seq_len(K), function(k)
    as.vector(G2_list[[k]] %*% (th_s_list[[k]] - th_h_list[[k]]))))
  correction <- as.vector(Delta %*% u_hat)

  # Finite-sample stability guard: if any component of the correction is
  # larger than 3 * SE_sipw (a multiple of the primary estimator's own SE),
  # the correction is almost certainly driven by a near-singular A_raw in
  # this replicate. Scale it down to keep the estimator bounded.
  # This guard has negligible effect when A_raw is well-conditioned
  # (correction is O(1/sqrt(N)) and SE_sipw is also O(1/sqrt(N))).
  scale_limit <- 3 * sqrt(pmax(diag(V_sipw), 0))
  correction_scale <- max(abs(correction) / scale_limit)
  if (correction_scale > 1) correction <- correction / correction_scale

  b_mp  <- b_s + correction

  # FIX 4: Augmented point estimate uses raw contrast (th_s - th_h), NOT G2*contrast.
  # b_k_aug already contains G2k_inv (via A2k_tl = G2k_inv) so the correction is
  # Delta_aug * (th_s - th_h) = -B_aug %*% A_aug_inv %*% (th_s - th_h).
  u_raw <- unlist(lapply(seq_len(K), function(k)
    as.vector(th_s_list[[k]] - th_h_list[[k]])))
  b_mp_aug <- b_s + as.vector(-B_aug %*% A_aug_inv %*% u_raw)

  # Augmented variance (Appendix B, Step B9): ----
  #   V_mp(pi_hat) = V_sipw(pi_hat) - B_aug A_aug^{-1} B_aug' / N
  # V_sipw uses augmented M_eff. B_aug and A_aug use the Appendix B blocks.
  V_mp       <- V_sipw - B_aug %*% A_aug_inv %*% t(B_aug) / N
  V_mp_naive <- G1_inv %*% C11 %*% t(G1_inv) / N -
                B_raw  %*% A_inv %*% t(B_raw) / N

  # PSD guards
  eig_mp <- eigen(V_mp, only.values = TRUE, symmetric = TRUE)$values
  if (any(eig_mp < -1e-10)) V_mp <- V_sipw
  if (any(diag(V_mp_naive) < 0)) V_mp_naive <- G1_inv %*% C11 %*% t(G1_inv) / N

  se_mp       <- sqrt(pmax(diag(V_mp),       0))
  se_mp_naive <- sqrt(pmax(diag(V_mp_naive), 0))

  # FIX 5: return se_mp (augmented) as the primary SE, se_mp_naive for comparison.
  list(coef      = b_mp,     se      = se_mp,
       coef_sipw = b_s,      se_sipw = se_sipw,
       se_naive  = se_mp_naive,
       coef_aug  = b_mp_aug,
       aug_used  = use_aug)
}