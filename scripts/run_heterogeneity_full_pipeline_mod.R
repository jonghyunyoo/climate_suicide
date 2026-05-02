# run_heterogeneity_full_pipeline_mod.R
# -----------------------------------------------------------------------------
# One-command runner for the revised heterogeneity analysis and table generation.
# Run from the project root, where the inputs/ directory is located.
# -----------------------------------------------------------------------------

script_candidates <- list(
  functions = c("scripts/functions_mod.R", "functions_mod.R"),
  region    = c("scripts/heterogeneity_region_mod.R", "heterogeneity_region_mod.R"),
  sex       = c("scripts/sex_mod.R", "sex_mod.R"),
  age       = c("scripts/age_mod.R", "age_mod.R"),
  tables    = c("scripts/make_heterogeneity_tables_mod.R", "make_heterogeneity_tables_mod.R")
)

find_script <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Could not find script. Tried: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

message("Running revised district-level heterogeneity models...")
source(find_script(script_candidates$region))

message("Running revised sex-specific heterogeneity models...")
source(find_script(script_candidates$sex))

message("Running revised age-specific heterogeneity models...")
source(find_script(script_candidates$age))

message("Building publication-ready heterogeneity tables...")
source(find_script(script_candidates$tables))

message("Pipeline complete.")
