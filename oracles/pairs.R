# Relationship-class pairs and their 2 * phi, for the gene-drop layer.
#
#   Rscript oracles/pairs.R <pedigree.parquet> <pairs.parquet>
#
# Writes one row per pair (0-based rows i, j), its class, and 2 * phi from
# kinship2 (which returns phi), with recorded generation 0 unrelated. MZ
# co-twins are one genetic individual for gene-drop, so MZ pairs are left out.

suppressPackageStartupMessages({
  library(arrow)
  library(kinship2)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "pedigree.R"))

args <- commandArgs(trailingOnly = TRUE)
ped <- read_ped(args[1])
classes <- classify(ped)
classes$MZ <- NULL
# Gene-drop scores parent-offspring IBD as one class.
classes <- c(classes[c("FS", "MHS", "PHS")], list(PO = rbind(classes$MO, classes$FO)), classes["1C"])

# MZ co-twins are one genetic individual: kinship2 sees only the lower-id
# twin, with the co-twin's children re-pointed to it (as in genedrop/drop.py),
# and the co-twin then copies its twin's kinship.
keeper <- twin_keeper(ped)
kept <- keeper == seq_len(nrow(ped))
mo <- keeper[ped$mo]
fa <- keeper[ped$fa]
K0 <- kinship(id = which(kept), dadid = fa[kept], momid = mo[kept], sex = ifelse(ped$sex[kept] == 1, 1, 2))
pos <- match(keeper, which(kept))
K <- as.matrix(K0)[pos, pos]

out <- do.call(rbind, lapply(names(classes), function(cl) {
  pr <- classes[[cl]]
  if (!nrow(pr)) {
    return(NULL)
  }
  data.frame(i = pr[, 1] - 1L, j = pr[, 2] - 1L, class = cl, phi2 = 2 * K[pr])
}))
write_parquet(out, args[2])
