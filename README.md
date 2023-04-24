# btbspread
Code for fitting dynamic occupancy model to bovine tuberculosis surveillance data in white-tailed deer

The code is organized into scripts and functions.

## Scripts (main)

### Feuka et al model fitting
Runs data cleaning/organization functions to load data into environment, specifies dynamic occupancy model using model-building function, fits model to data over multiuple chains in parallel, and saves outputs. Also calculates within-sample model fit evaluation post-hoc.

### Feuka et al figures
Uses model outputs to create non-pdf figures from manuscript

### Feuka et al cross validation
Prepares data and performs k-fold cross validation over all years in the study (leaving one year out for each fold), saves outputs, and calculates cross-validation metrics post-hoc.

## Functions

### dynOccSim_fun
Simulates data from a dynamic occupancy model with the same structure as the model. Outputs simulated data, covariates, simulation paramters, and simulated raster image off occupancy probabilities.

### fitDynOccMod_fun
Fits single chain of dynamic occupancy model using Markov-chian Monte Carlo algorithm. Moodel details described in Feuka et al. Can be fun in parallel (see script). outputs MCMC samples for all parameters.

### fitDynOccMod_loo_fun
Same MCMC algorithm as model-fitting function, but removes one fold of data. Intended to be run in a for loop or in parallel for all K folds. Outputs MCMC samples for all parameters used in cross-validation metrics.

### modDataPrep_fun
Prepares data for model fitting. Outputs list of data objects. Requires inputs from spatDataPrep_fun

### spatDataPrep_fun
Spatially references tabular data. Outputs spatial one dataframe by occupancy site level and year and one by individual sample.
