library(jsonlite)
library(GOSemSim)
library(GO.db)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(org.Mm.eg.db)
library(org.Dr.eg.db)

cat("Loading GO semantic data...\n")
hsGO_MF <- godata('org.Hs.eg.db', ont = "MF", computeIC = TRUE)
hsGO_BP <- godata('org.Hs.eg.db', ont = "BP", computeIC = TRUE)
hsGO_CC <- godata('org.Hs.eg.db', ont = "CC", computeIC = TRUE)
mmGO_MF <- godata('org.Mm.eg.db', ont = "MF", computeIC = TRUE)
mmGO_BP <- godata('org.Mm.eg.db', ont = "BP", computeIC = TRUE)
mmGO_CC <- godata('org.Mm.eg.db', ont = "CC", computeIC = TRUE)
scGO_MF <- godata('org.Sc.sgd.db', ont = "MF", computeIC = TRUE)
scGO_BP <- godata('org.Sc.sgd.db', ont = "BP", computeIC = TRUE)
scGO_CC <- godata('org.Sc.sgd.db', ont = "CC", computeIC = TRUE)
drGO_MF <- godata('org.Dr.eg.db', ont = "MF", computeIC = TRUE)
drGO_BP <- godata('org.Dr.eg.db', ont = "BP", computeIC = TRUE)
drGO_CC <- godata('org.Dr.eg.db', ont = "CC", computeIC = TRUE)
ecGO_MF <- hsGO_MF
ecGO_BP <- hsGO_BP
ecGO_CC <- hsGO_CC
cat("GO data loaded!\n")

get_semdata <- function(organism, namespace) {
  if (organism == "human") {
    if (namespace == "molecular_function") return(hsGO_MF)
    if (namespace == "biological_process") return(hsGO_BP)
    if (namespace == "cellular_component") return(hsGO_CC)
  } else if (organism == "mouse") {
    if (namespace == "molecular_function") return(mmGO_MF)
    if (namespace == "biological_process") return(mmGO_BP)
    if (namespace == "cellular_component") return(mmGO_CC)
  } else if (organism == "yeast") {
    if (namespace == "molecular_function") return(scGO_MF)
    if (namespace == "biological_process") return(scGO_BP)
    if (namespace == "cellular_component") return(scGO_CC)
  } else if (organism == "zebrafish") {
    if (namespace == "molecular_function") return(drGO_MF)
    if (namespace == "biological_process") return(drGO_BP)
    if (namespace == "cellular_component") return(drGO_CC)
  } else {
    if (namespace == "molecular_function") return(ecGO_MF)
    if (namespace == "biological_process") return(ecGO_BP)
    if (namespace == "cellular_component") return(ecGO_CC)
  }
}

extract_go_terms <- function(response) {
  matches <- regmatches(response, gregexpr("GO:[0-9]{7}", response))[[1]]
  return(matches)
}

extract_gt_terms <- function(ground_truth) {
  if (is.na(ground_truth) || ground_truth == "") return(character(0))
  matches <- regmatches(ground_truth, gregexpr("GO:[0-9]{7}", ground_truth))[[1]]
  return(matches)
}

compute_sim <- function(predicted, groundtruth, semdata) {
  if (length(predicted) == 0 || length(groundtruth) == 0) return(NA)
  tryCatch({
    mgoSim(predicted, groundtruth, semData = semdata, measure = "Wang", combine = "BMA")
  }, error = function(e) return(NA))
}

compute_random_baseline <- function(groundtruth, semdata, n_terms, n_reps = 10) {
  if (length(groundtruth) == 0) return(NA)
  all_terms <- keys(GOMFANCESTOR)
  sims <- c()
  for (i in 1:n_reps) {
    set.seed(i)
    random_terms <- sample(all_terms, min(n_terms, length(all_terms)))
    sim <- tryCatch({
      mgoSim(random_terms, groundtruth, semData = semdata, measure = "Wang", combine = "BMA")
    }, error = function(e) NA)
    sims <- c(sims, sim)
  }
  return(mean(sims, na.rm = TRUE))
}

process_model <- function(filepath, model_name) {
  cat(sprintf("Processing %s...\n", model_name))
  data <- stream_in(file(filepath))
  results <- data.frame()
  for (i in 1:nrow(data)) {
    row <- data[i, ]
    gt_terms <- extract_gt_terms(row$ground_truth)
    if (length(gt_terms) == 0) next
    pred_terms <- extract_go_terms(row$response)
    if (length(pred_terms) == 0) next
    semdata <- get_semdata(row$organism, row$namespace)
    llm_sim <- compute_sim(pred_terms, gt_terms, semdata)
    random_sim <- compute_random_baseline(gt_terms, semdata, length(pred_terms))
    results <- rbind(results, data.frame(
      model = model_name,
      accession = row$accession,
      gene_name = row$gene_name,
      organism = row$organism,
      namespace = row$namespace,
      prompt_id = row$prompt_id,
      llm_sim = llm_sim,
      random_sim = random_sim,
      stringsAsFactors = FALSE
    ))
    if (i %% 100 == 0) cat(sprintf("  Processed %d/%d rows\n", i, nrow(data)))
  }
  return(results)
}

base_path <- "/Users/kallolnaha/Documents/mix_projects/mindrouter/PROBE/probe_results/"

all_results <- rbind(
  process_model(paste0(base_path, "mistral-large-123b.jsonl"), "Mistral Large 123B"),
  process_model(paste0(base_path, "llama3.3-70b.jsonl"), "Llama 3.3 70B"),
  process_model(paste0(base_path, "qwen2.5-72b.jsonl"), "Qwen2.5 72B")
)

write.csv(all_results, "/Users/kallolnaha/Documents/mix_projects/mindrouter/PROBE/semantic_similarity_results.csv", row.names = FALSE)
cat("Results saved!\n")

summary_table <- aggregate(cbind(llm_sim, random_sim) ~ model + namespace,
                          data = all_results,
                          FUN = function(x) round(mean(x, na.rm=TRUE), 3))
print(summary_table)
