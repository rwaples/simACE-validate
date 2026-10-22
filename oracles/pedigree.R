# Pedigree reading, relationship classes and kinship, shared by the oracles.
# Shares no code with simACE or pedigree-graph.

read_ped <- function(path) {
  ped <- as.data.frame(read_parquet(path))
  n <- nrow(ped)
  stopifnot(all(ped$id == seq_len(n) - 1L))
  ped$mo <- ifelse(ped$mother < 0, NA_integer_, ped$mother + 1L)
  ped$fa <- ifelse(ped$father < 0, NA_integer_, ped$father + 1L)
  ped$tw <- ifelse(ped$twin < 0, NA_integer_, ped$twin + 1L)
  ped
}

# All unordered pairs within each group of row indices.
within_pairs <- function(groups) {
  groups <- groups[lengths(groups) >= 2]
  if (!length(groups)) {
    return(matrix(integer(0), ncol = 2))
  }
  do.call(rbind, lapply(groups, function(g) t(combn(sort(g), 2))))
}

pair_key <- function(pr, n) (pr[, 1] - 1) * n + pr[, 2]

# Relationship classes from the parent links: MZ, FS, MHS, PHS (both members
# in generation >= 1), MO (mother, child), FO (father, child), 1C (children of
# two full sibs, one parent each). Rows are 1-based; sib pairs are sorted (lo, hi).
classify <- function(ped) {
  n <- nrow(ped)
  nf <- which(ped$generation >= 1)
  by_mother <- within_pairs(split(nf, ped$mo[nf]))
  by_father <- within_pairs(split(nf, ped$fa[nf]))
  km <- pair_key(by_mother, n)
  kf <- pair_key(by_father, n)
  full <- by_mother[km %in% kf, , drop = FALSE]
  is_mz <- !is.na(ped$tw[full[, 1]]) & ped$tw[full[, 1]] == full[, 2]
  mz <- full[is_mz, , drop = FALSE]
  fs <- full[!is_mz, , drop = FALSE]
  mhs <- by_mother[!(km %in% kf), , drop = FALSE]
  phs <- by_father[!(kf %in% km), , drop = FALSE]
  mo <- cbind(ped$mo[nf], nf)
  fo <- cbind(ped$fa[nf], nf)

  # 1C: for each FS pair (a, b), children(a) x children(b).
  child_list <- vector("list", n)
  for (i in nf) {
    child_list[[ped$mo[i]]] <- c(child_list[[ped$mo[i]]], i)
    child_list[[ped$fa[i]]] <- c(child_list[[ped$fa[i]]], i)
  }
  fs_par <- fs[!vapply(child_list[fs[, 1]], is.null, TRUE) & !vapply(child_list[fs[, 2]], is.null, TRUE), , drop = FALSE]
  c1 <- if (nrow(fs_par)) {
    do.call(rbind, lapply(seq_len(nrow(fs_par)), function(r) {
      as.matrix(expand.grid(child_list[[fs_par[r, 1]]], child_list[[fs_par[r, 2]]]))
    }))
  } else {
    matrix(integer(0), ncol = 2)
  }
  if (nrow(c1)) {
    c1 <- t(apply(c1, 1, sort))
    sib_keys <- c(km, kf)
    c1 <- c1[!(pair_key(c1, n) %in% sib_keys) & c1[, 1] != c1[, 2], , drop = FALSE]
    c1 <- unique(c1)
  }
  list(MZ = mz, FS = fs, MHS = mhs, PHS = phs, MO = mo, FO = fo, `1C` = c1)
}

# MZ co-twins are one genetic individual: each row's keeper is itself, or the
# lower-id twin for a co-twin.
twin_keeper <- function(ped) {
  keeper <- seq_len(nrow(ped))
  co <- which(!is.na(ped$tw) & ped$tw < keeper)
  keeper[co] <- ped$tw[co]
  keeper
}

# Inbreeding coefficients from nadiv (makeAinv's f), on the pedigree with each
# co-twin folded into its twin and its children re-pointed to that twin; a
# co-twin then gets its twin's F. Generation 0 is unrelated and non-inbred:
# its hidden parents are unknown (ADR 0001 §3).
inbreeding <- function(ped) {
  keeper <- twin_keeper(ped)
  kept <- which(keeper == seq_len(nrow(ped)))
  folded <- data.frame(id = kept, dam = keeper[ped$mo[kept]], sire = keeper[ped$fa[kept]])
  nadiv::makeAinv(folded, det = FALSE)$f[match(keeper, kept)]
}
