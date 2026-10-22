# Layer-1 liability oracle (ADR 0001).
#
#   Rscript oracles/liability.R <input.json> <output.json>
#
# For each simACE pedigree, draws M independent realisations of A, C and E on
# the same pedigree, conditional on the recorded generation-0 values, and
# computes every statistic on simACE's values and on each draw with the same
# code. Writes, per pedigree and statistic, simACE's value and how many draws
# fall below and equal to it. Shares no code with simACE or pedigree-graph.
#
# The law, for generations g >= 1 (random mating):
#   A_i = (A_mother + A_father) / 2 + w_i,  w_i ~ N(0, d_i * A), corr rA across traits
#   C_i = bm * C_mother + bf * C_father + u_h for household h (the mother),
#         (u_h1, u_h2) ~ N(0, U), per trait bm, bf (simACE ADR 0025)
#   E_i ~ N(0, E), corr rE
# MZ co-twins share w and the household, not E. d_i = 1/2 ("constant") or
# 1/2 - (F_mother + F_father) / 4 ("infinitesimal"). Sigma is the configured C
# covariance (C1, C2, rC). U is, per generation:
#   "realised"  Sigma - the covariance of T = bm * C_mother + bf * C_father over
#               the generation's offspring rows (ADR 0025 section 3)
#   "analytic"  Sigma - diag(bm) Sigma diag(bm) - diag(bf) Sigma diag(bf)
# A negative eigenvalue of U is clipped to 0, as ADR 0025 does for sampling
# noise; the clips are counted. simACE instead stops a replicate whose excess
# is beyond its tolerance; the oracle has no such branch. Generation 1 is
# transmitted from recorded generation 0, stored as float32, where simACE used
# the unrounded values. With bm = bf = 0, C is one draw per household.

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
# nadiv is installed into .rlib/ by tools/build_nadiv.sh.
.libPaths(c(file.path(here, "..", ".rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(jsonlite)
  library(arrow)
  library(nadiv)
})
source(file.path(here, "pedigree.R"))

args <- commandArgs(trailingOnly = TRUE)
inp <- read_json(args[1], simplifyVector = TRUE)
p <- inp$params
M <- inp$M
mutation <- if (is.null(inp$mutation)) "none" else inp$mutation
stopifnot(mutation %in% c("none", "null", "seg_045", "rA_plus_01", "c_by_father", "bm_plus_01"))
stopifnot(inp$law %in% c("constant", "infinitesimal"))
c_law <- if (is.null(inp$c_law)) "realised" else inp$c_law
stopifnot(c_law %in% c("realised", "analytic"))
bm <- if (is.null(p$bm)) c(0, 0) else as.numeric(p$bm)
bf <- if (is.null(p$bf)) c(0, 0) else as.numeric(p$bf)
if (mutation == "bm_plus_01") bm <- bm + 0.1

# ---- statistics (columnwise over an n x M matrix) ---------------------------

col_cor <- function(u, v) {
  u <- sweep(u, 2, colMeans(u))
  v <- sweep(v, 2, colMeans(v))
  colSums(u * v) / sqrt(colSums(u * u) * colSums(v * v))
}

pair_r <- function(X, pr, symmetric) {
  a <- X[pr[, 1], , drop = FALSE]
  b <- X[pr[, 2], , drop = FALSE]
  if (symmetric) col_cor(rbind(a, b), rbind(b, a)) else col_cor(a, b)
}

col_cov <- function(X, Y) {
  X <- sweep(X, 2, colMeans(X))
  Y <- sweep(Y, 2, colMeans(Y))
  colSums(X * Y) / (nrow(X) - 1)
}

col_var <- function(X) col_cov(X, X)

# Pair correlations that are 1 by construction are not statistics. Maternal
# half-sibs share all of C only while the father transmits none of it.
degenerate <- function(cl, t) {
  switch(cl,
    MZ = c("A", "C"),
    FS = "C",
    MHS = if (bf[t] == 0) "C" else character(0),
    character(0)
  )
}

statistics <- function(V, ped, classes) {
  out <- list()
  for (t in 1:2) {
    comps <- list(A = V[[paste0("A", t)]], C = V[[paste0("C", t)]], E = V[[paste0("E", t)]])
    comps$L <- comps$A + comps$C + comps$E
    for (cl in names(classes)) {
      pr <- classes[[cl]]
      if (nrow(pr) < 10) next
      for (k in names(comps)) {
        if (k %in% degenerate(cl, t)) next
        out[[sprintf("r_%s%d_%s", k, t, cl)]] <- pair_r(comps[[k]], pr, symmetric = !(cl %in% c("MO", "FO")))
      }
    }
    for (g in sort(unique(ped$generation[ped$generation >= 1]))) {
      rows <- which(ped$generation == g)
      out[[sprintf("var_A%d_g%d", t, g)]] <- col_var(comps$A[rows, , drop = FALSE])
      out[[sprintf("var_C%d_g%d", t, g)]] <- col_var(comps$C[rows, , drop = FALSE])
      out[[sprintf("var_L%d_g%d", t, g)]] <- col_var(comps$L[rows, , drop = FALSE])
    }
  }
  L1 <- V$A1 + V$C1 + V$E1
  L2 <- V$A2 + V$C2 + V$E2
  for (g in sort(unique(ped$generation[ped$generation >= 1]))) {
    rows <- which(ped$generation == g)
    for (k in c("A", "C", "E")) {
      out[[sprintf("r_%s12_g%d", k, g)]] <- col_cor(V[[paste0(k, 1)]][rows, , drop = FALSE], V[[paste0(k, 2)]][rows, , drop = FALSE])
    }
    out[[sprintf("r_L12_g%d", g)]] <- col_cor(L1[rows, , drop = FALSE], L2[rows, , drop = FALSE])
  }
  out
}

# ---- segregation variance ---------------------------------------------------

# d_i, the Mendelian-sampling variance of A_i in units of A.
segregation <- function(ped) {
  if (inp$law == "constant") {
    return(rep(if (mutation == "seg_045") 0.45 else 0.5, nrow(ped)))
  }
  F <- inbreeding(ped)
  d <- rep(NA_real_, nrow(ped))
  nf <- !is.na(ped$mo)
  d[nf] <- 0.5 - (F[ped$mo[nf]] + F[ped$fa[nf]]) / 4
  d
}

# ---- the oracle draw --------------------------------------------------------

# Two standard-normal n x M matrices with correlation r.
corr_normals <- function(n, M, r) {
  z1 <- matrix(rnorm(n * M), n, M)
  z2 <- r * z1 + sqrt(1 - r^2) * matrix(rnorm(n * M), n, M)
  list(z1, z2)
}

# U for one generation, one value per draw column: the configured C covariance
# minus the covariance of the transmitted parts T1, T2 (rows x M).
residual_c <- function(T1, T2) {
  S <- if (c_law == "realised") {
    list(col_var(T1), col_var(T2), col_cov(T1, T2))
  } else {
    lapply(list(
      (bm[1]^2 + bf[1]^2) * p$C1, (bm[2]^2 + bf[2]^2) * p$C2,
      (bm[1] * bm[2] + bf[1] * bf[2]) * p$rC * sqrt(p$C1 * p$C2)
    ), rep, ncol(T1))
  }
  list(u11 = p$C1 - S[[1]], u22 = p$C2 - S[[2]], u12 = p$rC * sqrt(p$C1 * p$C2) - S[[3]])
}

# One N(0, U) draw per household and column, with U's negative eigenvalues
# clipped to 0. Returns the two h x M matrices and the number of columns clipped.
household_c <- function(U, h) {
  mid <- (U$u11 + U$u22) / 2
  rad <- sqrt(((U$u11 - U$u22) / 2)^2 + U$u12^2)
  theta <- atan2(2 * U$u12, U$u11 - U$u22) / 2
  s1 <- sqrt(pmax(mid + rad, 0))
  s2 <- sqrt(pmax(mid - rad, 0))
  M <- length(mid)
  z1 <- matrix(rnorm(h * M), h, M)
  z2 <- matrix(rnorm(h * M), h, M)
  col <- function(v) rep(v, each = h)
  list(
    z1 * col(s1 * cos(theta)) - z2 * col(s2 * sin(theta)),
    z1 * col(s1 * sin(theta)) + z2 * col(s2 * cos(theta)),
    clipped = sum(mid - rad < 0)
  )
}

oracle_draw <- function(ped, M) {
  n <- nrow(ped)
  sdA <- sqrt(c(p$A1, p$A2))
  sdC <- sqrt(c(p$C1, p$C2))
  sdE <- sqrt(c(p$E1, p$E2))
  d_seg <- segregation(ped)
  rA <- if (mutation == "rA_plus_01") p$rA + 0.1 else p$rA
  V <- list()
  g0 <- ped$generation == 0
  for (k in c("A1", "C1", "E1", "A2", "C2", "E2")) {
    X <- matrix(NA_real_, n, M)
    X[g0, ] <- ped[[k]][g0]
    V[[k]] <- X
  }
  hh_parent <- if (mutation == "c_by_father") ped$fa else ped$mo
  transmits <- any(bm > 0 | bf > 0)
  clipped <- 0
  for (g in sort(unique(ped$generation[!g0]))) {
    rows <- which(ped$generation == g)
    m <- length(rows)
    w <- corr_normals(m, M, rA)
    tw <- ped$tw[rows]
    copy <- which(!is.na(tw) & tw < rows)
    for (t in 1:2) {
      wt <- sqrt(d_seg[rows]) * sdA[t] * w[[t]]
      wt[copy, ] <- wt[match(tw[copy], rows), ]
      A <- V[[paste0("A", t)]]
      A[rows, ] <- (A[ped$mo[rows], , drop = FALSE] + A[ped$fa[rows], , drop = FALSE]) / 2 + wt
      V[[paste0("A", t)]] <- A
    }
    hh <- match(hh_parent[rows], unique(hh_parent[rows]))
    if (transmits) {
      tr <- lapply(1:2, function(t) {
        Ct <- V[[paste0("C", t)]]
        bm[t] * Ct[ped$mo[rows], , drop = FALSE] + bf[t] * Ct[ped$fa[rows], , drop = FALSE]
      })
      u <- household_c(residual_c(tr[[1]], tr[[2]]), max(hh))
      clipped <- clipped + u$clipped
      for (t in 1:2) V[[paste0("C", t)]][rows, ] <- tr[[t]] + u[[t]][hh, , drop = FALSE]
    } else {
      # The same law with U = Sigma. Kept as its own branch so cases without
      # transmission draw the same random numbers as before transmission existed.
      cz <- corr_normals(max(hh), M, p$rC)
      for (t in 1:2) V[[paste0("C", t)]][rows, ] <- sdC[t] * cz[[t]][hh, , drop = FALSE]
    }
    ez <- corr_normals(m, M, p$rE)
    for (t in 1:2) V[[paste0("E", t)]][rows, ] <- sdE[t] * ez[[t]]
  }
  attr(V, "c_clipped") <- clipped
  V
}

# The statistics rebuild liability as A + C + E, so first check that simACE's
# stored liability is that sum (float32 storage rounds each column).
observed <- function(ped) {
  V <- list()
  for (k in c("A1", "C1", "E1", "A2", "C2", "E2")) V[[k]] <- matrix(as.double(ped[[k]]), ncol = 1)
  for (t in 1:2) {
    sum_t <- V[[paste0("A", t)]] + V[[paste0("C", t)]] + V[[paste0("E", t)]]
    gap <- max(abs(as.double(ped[[paste0("liability", t)]]) - sum_t))
    if (gap > 1e-5) stop(sprintf("simACE liability%d differs from A + C + E by up to %g", t, gap))
  }
  V
}

# ---- main -------------------------------------------------------------------

set.seed(inp$oracle_seed)
result <- list()
for (s in seq_along(inp$pedigrees)) {
  ped <- read_ped(inp$pedigrees[[s]])
  classes <- classify(ped)
  # "null" replaces simACE by one more oracle draw: the power study's false-positive check.
  obs_values <- if (mutation == "null") oracle_draw(ped, 1) else observed(ped)
  obs_clipped <- if (is.null(attr(obs_values, "c_clipped"))) 0 else attr(obs_values, "c_clipped")
  obs <- statistics(obs_values, ped, classes)
  draws <- oracle_draw(ped, M)
  orc <- statistics(draws, ped, classes)
  rows <- lapply(names(obs), function(k) {
    x <- obs[[k]]
    y <- orc[[k]]
    list(
      stat = k, simace = x, n_less = sum(y < x), n_equal = sum(y == x),
      oracle_mean = mean(y), oracle_sd = sd(y)
    )
  })
  result[[s]] <- list(
    seed = inp$seeds[[s]],
    n_pairs = lapply(classes, nrow),
    c_clipped = attr(draws, "c_clipped") + obs_clipped,
    rows = rows
  )
}
write_json(result, args[2], auto_unbox = TRUE, digits = NA)
