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
# The law, for generations g >= 1 (no C transmission, no assortment):
#   A_i = (A_mother + A_father) / 2 + w_i,  w_i ~ N(0, d_i * A), corr rA across traits
#   C_i = c_h for household h (the mother), c_h ~ N(0, C), corr rC
#   E_i ~ N(0, E), corr rE
# MZ co-twins share w and the household, not E. d_i = 1/2 ("constant") or
# 1/2 - (F_mother + F_father) / 4 ("infinitesimal").

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
stopifnot(mutation %in% c("none", "null", "seg_045", "rA_plus_01", "c_by_father"))
stopifnot(inp$law %in% c("constant", "infinitesimal"))

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

col_var <- function(X) {
  X <- sweep(X, 2, colMeans(X))
  colSums(X * X) / (nrow(X) - 1)
}

# Pair correlations that are 1 by construction are not statistics.
degenerate <- list(MZ = c("A", "C"), FS = "C", MHS = "C")

statistics <- function(V, ped, classes) {
  out <- list()
  for (t in 1:2) {
    comps <- list(A = V[[paste0("A", t)]], C = V[[paste0("C", t)]], E = V[[paste0("E", t)]])
    comps$L <- comps$A + comps$C + comps$E
    for (cl in names(classes)) {
      pr <- classes[[cl]]
      if (nrow(pr) < 10) next
      for (k in names(comps)) {
        if (k %in% degenerate[[cl]]) next
        out[[sprintf("r_%s%d_%s", k, t, cl)]] <- pair_r(comps[[k]], pr, symmetric = cl != "PO")
      }
    }
    for (g in sort(unique(ped$generation[ped$generation >= 1]))) {
      rows <- which(ped$generation == g)
      out[[sprintf("var_A%d_g%d", t, g)]] <- col_var(comps$A[rows, , drop = FALSE])
      out[[sprintf("var_L%d_g%d", t, g)]] <- col_var(comps$L[rows, , drop = FALSE])
    }
  }
  for (g in sort(unique(ped$generation[ped$generation >= 1]))) {
    rows <- which(ped$generation == g)
    for (k in c("A", "C", "E")) {
      out[[sprintf("r_%s12_g%d", k, g)]] <- col_cor(V[[paste0(k, 1)]][rows, , drop = FALSE], V[[paste0(k, 2)]][rows, , drop = FALSE])
    }
    L1 <- V$A1 + V$C1 + V$E1
    L2 <- V$A2 + V$C2 + V$E2
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
    cz <- corr_normals(max(hh), M, p$rC)
    ez <- corr_normals(m, M, p$rE)
    for (t in 1:2) {
      V[[paste0("C", t)]][rows, ] <- sdC[t] * cz[[t]][hh, , drop = FALSE]
      V[[paste0("E", t)]][rows, ] <- sdE[t] * ez[[t]]
    }
  }
  V
}

observed <- function(ped) {
  V <- list()
  for (k in c("A1", "C1", "E1", "A2", "C2", "E2")) V[[k]] <- matrix(as.double(ped[[k]]), ncol = 1)
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
  obs <- statistics(obs_values, ped, classes)
  orc <- statistics(oracle_draw(ped, M), ped, classes)
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
    rows = rows
  )
}
write_json(result, args[2], auto_unbox = TRUE, digits = NA)
