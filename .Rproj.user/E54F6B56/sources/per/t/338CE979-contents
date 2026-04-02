library(rstan)
library(tidyverse)


source("1.generating_data.R")
source("2.summary_result.R")
source("4.linegraph_CI.R")

## calculate coverage of the true parameter values

compute_coverage <- function(sim_results, true_vals) {
  
  # Build one long data frame from all simulation summaries
  summary_long_df <- map_dfr(seq_along(sim_results), function(i) {
    res <- sim_results[[i]]
    
    # Filter unwanted parameters
    keep_params <- !grepl("log_lik|lp_", names(res$estimates))
    param_names <- names(res$estimates)[keep_params]
    
    tibble(
      sim_id = i,
      parameter = param_names,
      mean = res$estimates[keep_params],
      lower = res$lower_ci[keep_params],
      upper = res$upper_ci[keep_params]
    )
  })
  
  # Filter only parameters with known true values
  coverage_df <- summary_long_df %>%
    filter(parameter %in% names(true_vals)) %>%
    mutate(
      true_value = true_vals[parameter],
      covered = (lower <= true_value) & (upper >= true_value)
    )
  
  # Summarize coverage rate
  coverage_summary <- coverage_df %>%
    group_by(parameter) %>%
    summarise(
      coverage = mean(covered),
      .groups = "drop"
    )
  
  return(coverage_summary)
}


# set names for vector of true values
true_wb.s <- true_wb.rp <- true_wb.srp
true_ln.s <- true_ln.rp <- true_ln.srp
true_llg.s <- true_llg.rp <- true_llg.srp

names(true_wb.l) <- names(true_ln.l) <- names(true_llg.l) <-c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "sigma")
names(true_wb.s) <- names(true_ln.s) <- names(true_llg.s) <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "alpha_s", "sigma") # for slogit
names(true_wb.rp) <- names(true_ln.rp) <- names(true_llg.rp) <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "lambda_rp", "sigma") # for rplogit
names(true_wb.fg) <- names(true_ln.fg) <- names(true_llg.fg) <- c("betaU[1]", "betaU[2]", "betaU[3]", "betaC[1]", "betaC[2]", "alpha_fg", "lambda_fg", "sigma") # for fglogit




# example use ###################################################################################################################

# wb1: WBL, wb2: WBS, wb3: WBRP, wb4: WBFG
# ln1: LNL, ln2: LNS, ln3: LNRP, ln4: LNFG
# llg1: LLGL, llg2: LLGS, llg3: LLGRP, llg4: LLGFG

coverage_wbfg <- lapply(list(simulation_summaries.wb4.500), compute_coverage, true_vals = true_wb.fg)

print(coverage_wbfg)

