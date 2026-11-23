#############################################################
# run_all.R  -  runs the full analysis in order
# Usage (from this folder):  Rscript run_all.R
#############################################################
steps <- c("R/01_prepare_data.R",
           "R/02_diversity.R",
           "R/03_stability.R",
           "R/04_differential_abundance.R",
           "R/05_clinical_associations.R",
           "R/06_networks.R",
           "R/07_sensitivity.R")
for (s in steps) {
  cat("\n==============================\nRunning", s, "\n==============================\n")
  source(s, local = new.env())
}
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("\nAll steps finished. Package versions written to sessionInfo.txt\n")