
library(GOSemSim)
library(GO.db)

hsGO_MF <- godata("org.Hs.eg.db", ont = "MF", computeIC = TRUE)
hsGO_BP <- godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE)
hsGO_CC <- godata("org.Hs.eg.db", ont = "CC", computeIC = TRUE)

compute_random_ns <- function(semdata, ancestor_db, n_terms=6, n_reps=50) {
  all_terms <- keys(ancestor_db)
  sims <- c()
  for (i in 1:n_reps) {
    set.seed(i)
    random_terms <- sample(all_terms, min(n_terms, length(all_terms)))
    gt_terms <- sample(all_terms, 3)
    sim <- tryCatch({
      mgoSim(random_terms, gt_terms, semData = semdata, measure = "Wang", combine = "BMA")
    }, error = function(e) NA)
    sims <- c(sims, sim)
  }
  return(mean(sims, na.rm = TRUE))
}

cat("Computing MF random baseline...
")
mf_random <- compute_random_ns(hsGO_MF, GOMFANCESTOR)
cat(sprintf("MF random baseline: %.3f
", mf_random))

cat("Computing BP random baseline...
")
bp_random <- compute_random_ns(hsGO_BP, GOBPANCESTOR)
cat(sprintf("BP random baseline: %.3f
", bp_random))

cat("Computing CC random baseline...
")
cc_random <- compute_random_ns(hsGO_CC, GOCCANCESTOR)
cat(sprintf("CC random baseline: %.3f
", cc_random))
