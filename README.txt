README: Code and Data for Simulation and Real Data Analysis

This repository contains the R code, Stan models, and data required to replicate the simulation studies and real data analysis presented in the manuscript. The project includes Bayesian mixture cure rate models with different AFT latency distributions and link functions for the cure proportion, as well as associated model assessment and diagnostics. 

Open the `Supplementary Material.Rproj` before running any scripts. This sets the working directory and ensures relative paths function correctly. 

If the results/ subfolders are not present,
create the following folders in the same directory as `Supplementary Material.Rproj` before running any scripts:

  results/data_result/
  results/result500/
  results/result1000/


========================================
Required R Packages:
========================================
install.packages(c("rstan", "tidyverse", "loo", "survminer",
                   "survival", "ggpubr", "ggplot2"))


========================================
R Scripts Overview (Files 1–7):
========================================
1. generating_data.R
   └── Functions for generating datasets for simulation studies

2. summary_result.R
   └── Loads Stan files, extracts and summarizes key performance
       metrics from simulation output files

3. Model assessment scripts
   ├── 3_1.Weibull_model_assessment.R
   ├── 3_2.Lognormal_model_assessment.R
   └── 3_3.Loglogistic_model_assessment.R
       └── Run simulations for diagnosing model assessment criteria

4. linegraph_CI.R
   └── Defines functions to generate line plots of 95% credible
       intervals across simulation replications

5. single_cox_snell.R
   └── Generates Cox-Snell residual KM plots for a single
       simulation replicate

6. coverage.R
   └── Computes coverage probabilities for all model parameters
       based on simulation results

7. augmented_cox_snell.R
   └── Constructs Cox-Snell residual KM plots using the full
       set of simulation replications



========================================
 FOLDER STRUCTURE
========================================

SIMULATION CODE
├── sim500/
│   ├── R scripts for simulations with sample size n = 500
│   └── Scripts follow naming: sim_[model].500.R
│       where [model] is wb1–wb4, ln1–ln4, llg1–llg4
└── sim1000/
    ├── R scripts for simulations with sample size n = 1000
    └── Scripts follow naming: sim_[model].1000.R
        where [model] is wb1–wb4, ln1–ln4, llg1–llg4

SIMULATION AND REAL DATA APPLICATION RESULTS
├── data_rcode/
│   ├── iowa18.csv (real dataset)
│   ├── *_stan.R: R scripts to fit twelve model specifications
│   └── iowa_coxsnell_ppp.R: Cox-Snell residual KM plots for fitted models
├── stan_final/
│   └── Stan model files for Bayesian mixture cure rate models
└── results/
    ├── data_result/ 
    │   └── Populated after running *_stan.R scripts in data_rcode/
    ├── result500/
    │   └── Populated after running all scripts in sim500/
    └── result1000/
        └── Populated after running all scripts in sim1000/



========================================
Reproducing the Simulation Results:
========================================
NOTE: Pre-computed simulation results (.rds files) are not 
included in this repository due to file size constraints.
To reproduce all results from scratch, follow the steps below.

1. Run generating_data.R and summary_result.R
2. Run ALL scripts in sim500/ and sim1000/ to generate 
   the required .rds result files
   (24 scripts total: sim_[model].500.R and sim_[model].1000.R
   for wb1–wb4, ln1–ln4, llg1–llg4)
3. Once all .rds files from Step 2 are generated, the 
   following scripts can be run:
   - 4. linegraph_CI.R
   - 6. coverage.R
   - 7. augmented_cox_snell.R
4. The following scripts can be run independently 
   without prior steps:
   - 3_1.Weibull_model_assessment.R
   - 3_2.Lognormal_model_assessment.R
   - 3_3.Loglogistic_model_assessment.R
   - 5. single_cox_snell.R

Note: Step 2 is computationally intensive. Each script runs
1000 simulation replications and may take several hours.
Running on an HPC cluster is recommended, or adjust
num_reps in the script to reduce the number of replications.


========================================
Reproducing the Real data analysis:
========================================
NOTE: Pre-computed model fit results (.rds files) are not
included in this repository due to file size constraints.
To reproduce all results from scratch, follow the steps below.

1. Navigate to the data_rcode/ directory
2. Run ALL twelve model fitting scripts (*_stan.R) to generate
   the required .rds result files saved in results/data_result/
   (12 scripts total: one for each combination of
   latency distribution and link function)
3. Once all .rds files are generated, run iowa_coxsnell_ppp.R
   to generate Cox-Snell residual KM plots and posterior
   predictive p-values for the fitted models


========================================
Supplementary Material:
========================================
- Supplementary_JDS.pdf  
  - Supplementary figures and tables referenced in the manuscript.
