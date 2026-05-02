# heterogeneity_master_mod.R
# Convenience runner for all revised heterogeneity scripts.
# Place the *_mod.R files in the same directory, or in the scripts/ directory.

script_candidates <- function(fname) {
  candidates <- c(file.path("scripts", fname), fname)
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Cannot find ", fname)
  hit
}

source(script_candidates("heterogeneity_region_mod.R"))
source(script_candidates("sex_mod.R"))
source(script_candidates("age_mod.R"))
