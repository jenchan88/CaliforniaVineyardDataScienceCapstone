Files needed for running module input document:

Function files:
- abundance_helper_functions.r: contains functions to generate correlation matrices and differential correlation matrix
- functionality_helper_functions.r: contains functions for generating functions for the modules
- correlation_contribution_helper_functions.r: contains functions for correlation contribution and edaphic data

Data files:
- "rarified_counts_by_(taxa level).csv": data file containing rarified counts for different classes/order/... under all treatments and all replicates.

- "filtered_modules_rds": folder with rds files for each combination of treatments for each experiment. Files are of form "T1-T2.rds", ex "H-NH.rds"

- sample module files, for samples of size 2-8 in rds form, in the a folder called "sample_modules" containing pre-sampled functions for significance check on functions in the module.

- "rainfall_data.csv": raw rainfall data for the SLO Airport

- "cndata.xlsx" containing 3 separate sheets for each year, with each sheet being labeled as the year the data was collected. So the sheet "2016" would contain the 2016 CN data.


We also have a fully automated version that doesn't require the pre-sampled module functionalities that can be used on taxonomy levels other than the class level. This document takes a long time to run, usually 3-4 minutes.


