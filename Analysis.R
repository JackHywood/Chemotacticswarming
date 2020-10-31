# `ChemotaxisCharacteriser` can be used to analyse any agent-based process consisting of a agents 
# performing a stochastic random walk-like motility process within a constrained 2D circular region,
# and for which the evolution of the agent population is well approximated by an advection-diffusion 
# partial differential equation (PDE).
# 
# `Analysis.R` analyses the provided example data sets found in the `Data` folder.
#
# Please see README.txt file for detailed explanation of the components of this script.
#
# To use:
# - Download repository, move to Desktop, and ensure the folder is titled `Chemotacticswarming`.
# - Uncomment the `dataset` you wish to analyse and ensure remaining `datasets` remain commented.
# - If desired, change values for `rR`, `rB`, and input variables called by the various functions.
# - Source `Analysis.R`
#
# Jack Hywood and Mark N. Read, 2019.

rm(list=ls())  # Ensure a clean environment.

require(fda)
require(rainbow)
require(fields)
require(refund)

setwd("~/Desktop/Chemotacticswarming/")
data_location = "~/Desktop/Chemotacticswarming/Data/"

# Uncomment one dataset for analaysis:
# Simulations:

dataset = "Unbiased"
# dataset = "Attraction_const"
seed = 1 # 1, 2, 3

# dataset = "Attraction_var"

# Experiments:

# dataset = "Preembedcognate"  # See if statements below for text to put here, to run a specific dataset.
# dataset = "Noncognate"  # See if statements below for text to put here, to run a specific dataset.

if (dataset == "Unbiased")
{
  source('FunEstim_const.R')
  rR = 100 # Radius of R.
  rB = 20 # Radius of B.
  if(seed == 1){
    working_dir = paste(data_location, "/RandWalk/Unbiased/seed1/", sep="") 
  }else if(seed == 2){
    working_dir = paste(data_location, "/RandWalk/Unbiased/seed2/", sep="")
  }else if(seed == 3){
    working_dir = paste(data_location, "/RandWalk/Unbiased/seed3/", sep="")
  }
  result = FunEstim_const(
    working_dir = working_dir,
    tmin_iter = NULL,
    tmax_iter = NULL,
    rR = rR,
    rB = rB,
    r_resolution = 0.1,
    normalise = TRUE,
    splines_K = 40,
    lambda_K = NULL,
    min_loglam = -5,
    max_loglam = 10,
    Norder = 6,
    Lfdval = 4,
    rmin_lambda_plot = 5,
    rmax_lambda_plot = rR-5,
    kymograph_zlim = NULL,
    rmin_est = 5,
    rmax_est = rR-2,
    rmin_plot = rB,
    rmax_plot = rR,
    splines_f = 20,
    splines_D = 20,
    lambda_f = 10^2,
    lambda_D = 10^2
  )
  r_trunc_values_plot = result$r_trunc_values_plot
  f_df = result$f_df
  D_df = result$D_df
  rmin_plot = result$rmin_plot
  rmax_plot = result$rmax_plot
  source(paste(data_location, '/RandWalk/Unbiased/Unbiased_plots.R', sep=''))
}

if (dataset == "Attraction_const")
{
  source('FunEstim_const.R')
  rR = 100
  rB = 20
  if(seed == 1){
    working_dir = paste(data_location, "/RandWalk/Attraction_const/seed1/", sep="") 
  }else if(seed == 2){
    working_dir = paste(data_location, "/RandWalk/Attraction_const/seed2/", sep="")
  }else if(seed == 3){
    working_dir = paste(data_location, "/RandWalk/Attraction_const/seed3/", sep="")
  }
  result = FunEstim_const(
    working_dir = working_dir,
    tmin_iter = NULL,
    tmax_iter = NULL,
    rR = rR,
    rB = rB,
    r_resolution = 0.1,
    normalise = TRUE,
    splines_K = 40,
    lambda_K = NULL,
    min_loglam = -10,
    max_loglam = 10,
    Norder = 6,
    Lfdval = 4,
    rmin_lambda_plot = 5,
    rmax_lambda_plot = rR-5,
    kymograph_zlim = NULL,
    rmin_est = 5,
    rmax_est = rR-2,
    rmin_plot = rB,
    rmax_plot = rR,
    splines_f = 20,
    splines_D = 20,
    lambda_f = 10^2,
    lambda_D = 10^2
  )
  r_trunc_values_plot = result$r_trunc_values_plot
  f_df = result$f_df
  D_df = result$D_df
  rmin_plot = result$rmin_plot
  rmax_plot = result$rmax_plot
  source(paste(data_location, '/RandWalk/Attraction_const/Attraction_const_plots.R', sep=''))
}

if (dataset == "Attraction_var")
{
  source('FunEstim_var.R')
  rR = 100
  rB = 20
  working_dir = paste(data_location, "/RandWalk/Attraction_var/", sep="")
  result = FunEstim_var(
    working_dir = working_dir,
    tmin_iter = NULL,
    tmax_iter = NULL,
    half_time_window_width = NULL,
    rR = rR,
    rB = rB,
    r_resolution = 0.1,
    normalise = TRUE,
    splines_K = 40,
    lambda_K = NULL,
    min_loglam = -5,
    max_loglam = 10,
    Norder = 6,
    Lfdval = 4,
    rmin_lambda_plot = 10,
    rmax_lambda_plot = rR-5,
    kymograph_zlim = NULL,
    rmin_est = 5,
    rmax_est = rR-2,
    rmin_plot = rB,
    rmax_plot = rR,
    splines_f = 20,
    splines_D = 20,
    lambda_f = 10^1,
    lambda_D = 10^1
  )
  r_trunc_values_plot = result$r_trunc_values_plot
  f_fds = result$f_fds
  f_var = result$f_var
  D_fds = result$D_fds
  D_var = result$D_var
  rmin_plot = result$rmin_plot
  rmax_plot = result$rmax_plot
  source(paste(data_location, '/RandWalk/Attraction_var/Attraction_var_plots.R', sep=''))
}

if (dataset == "Preembedcognate")
{
  source('FunEstim_const.R')
  source('FunEstim_var.R')
  rR = 2985
  rB = (1082.5+1115.5)/2
  working_dir = paste(data_location, "/Preembedcognate/", sep="")
  fD =  FunEstim_const(
    working_dir = working_dir,
    tmin_iter = 11,
    tmax_iter = 100,
    rR = rR,
    rB = rB,
    r_resolution = 1,
    normalise = TRUE,
    splines_K = 250,
    lambda_K = NULL,
    min_loglam = -10,
    max_loglam = 10,
    Norder = 6,
    Lfdval = 4,
    rmin_lambda_plot = 100,
    rmax_lambda_plot = rR-50,
    kymograph_zlim = c(-0.5*10^-4,9*10^-4),
    rmin_est = 5,
    rmax_est = rR-5,
    rmin_plot = rB,
    rmax_plot = rR,
    splines_f = 20,
    splines_D = 20,
    lambda_f = 10^1,
    lambda_D = 10^1
  )
  
  FunEstim_var(
    working_dir = working_dir,
    tmin_iter = 11,
    tmax_iter = 100,
    half_time_window_width = NULL,
    rR = rR,
    rB = rB,
    r_resolution = 5,
    normalise = TRUE,
    splines_K = 250,
    lambda_K = NULL,
    min_loglam = -5,
    max_loglam = 10,
    Norder = 6,
    Lfdval = 4,
    rmin_lambda_plot = 100,
    rmax_lambda_plot = rR-50,
    kymograph_zlim = c(-0.5*10^-4,9*10^-4),
    rmin_est = 5,
    rmax_est = rR-5,
    rmin_plot = rB,
    rmax_plot = rR,
    splines_f = 20,
    splines_D = 20,
    lambda_f = 10^1,
    lambda_D = 10^1
  )
}

if (dataset == "Noncognate")
{
  source('FunEstim_const.R')
  rR = 2855
  rB = (1043+1110.5)/2
  working_dir = paste(data_location, "/Noncognate/", sep="")
  FunEstim_const(
    working_dir = working_dir,
    tmin_iter = NULL,
    tmax_iter = 97,
    rR = rR,
    rB = rB,
    r_resolution = 1,
    normalise = TRUE,
    splines_K = 250,
    lambda_K = NULL,
    min_loglam = -5,
    max_loglam = 10,
    Norder = 6,
    Lfdval = 4,
    rmin_lambda_plot = 100,
    rmax_lambda_plot = rR-50,
    kymograph_zlim = NULL,
    rmin_est = 5,
    rmax_est = rR-5,
    rmin_plot = rB,
    rmax_plot = rR,
    splines_f = 20,
    splines_D = 20,
    lambda_f = 10^2,
    lambda_D = 10^2
  )
}

