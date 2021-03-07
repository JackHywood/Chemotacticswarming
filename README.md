# Chemotacticswarming

## Table of Contents

   * [Intro and context](#intro-and-context)
   * [What the analysis does](#what-the-analysis-does)
   * [Usage](#usage)
      * [Input data format](#input-data-format)
      * [Entry point to analysis](#entry-point-to-analysis)
      * [Usage tips](#usage-tips)
   * [Module FunEstim_const.R](#module-funestim_constr)
      * [Inputs for module FunEstim_const.R](#inputs-for-module-funestim_constr)
      * [Outputs for module FunEstim_const.R](#outputs-for-module-funestim_constr)
      * [Notes on choice of inputs for module FunEstim_const.R](#notes-on-choice-of-inputs-for-module-funestim_constr)
   * [Module FunEstim_var.R](#module-funestim_varr)
      * [Inputs for module FunEstim_var.R](#inputs-for-module-funestim_varr)
      * [Outputs for module FunEstim_var.R](#outputs-for-module-funestim_varr)
      * [Guidance on setting input parameters](#guidance-on-setting-input-parameters)
   * [Module Prelim.R; Preprocessing and preliminary analysis](#module-prelimr-preprocessing-and-preliminary-analysis)
      * [Function read_distance](#function-read_distance)
      * [Function Chemotax_preprocessing](#function-chemotax_preprocessing)
      * [Function initialise_Kstep](#function-initialise_kstep)
      * [Function focal_L_M](#function-focal_l_m)
      * [Function plot_surface](#function-plot_surface)
   * [Module Analysis.R](#module-analysisr)
   * [Example data](#example-data)
      * [Simulations](#simulations)
      * [In vitro experiments](#in-vitro-experiments)
         * [Non-Cognate experiment](#non-cognate-experiment)
         * [Pre-Embedded-Cognate experiment](#pre-embedded-cognate-experiment)

## Intro and context

This analysis is described in detail in the academic article (which is currently in press):

```
Detection and characterisation of chemotactic swarming without cell tracking. Jack D. Hywood, Gregory Rice, Sophie V. Pageon, Mark N. Read, Maté Biro.
```

`Chemotacticswarming` can be used to analyse a system of motile "agents" performing a random walk, and it decomposes the population's collective movement into directed motion (*chemotaxis*) and diffusion components.
We've provided an implementation of the method here, and sample data so you can see how to get going with your own (more details below).

Here are two key benefits to our approach:

1. You don't need to visualise the chemokine field/factors/molecules. This is difficult to do, and any effort will likely disrupt agent movements.

2. You don't need to track individual agents between time frames. Our analysis only requires agent positions, we don't need to know which agent was which. 

There are a few keys assumptions/restrictions to our method.

1. The analysis assumes a circular environment (which we call `R`), and chemotaxis is analysed with respect to the centre.
   
   We expect that a slight deviation from a circular geometry will not greatly impact analysis.
   
   Our approach deals with both a closed and open boundary for `R`. If your environment is not circular, you can filter out any agent positions beyond a circular region that you define. We expect our method can be extended to deal with non-circular environments by adding an edge correction term as is normally done for spatial statistics.

2. The analysis assumes that agent population movement is well-approximated by an advection-diffusion partial differential equation (PDE), the *Fokker-Planck* equation - this stems from an assumption in our method. It's worked well for everything we've tried, but just be warned that if your population motility departs substantially from such an equation, the results could be non-representative.

## What the analysis does

The analysis focusses on the behaviour of agents within `R` with respect to a smaller concentric circular region `B`: are agents heading into or away from B?
As an example, we might be interested in the behaviour of immune cells outside of a tumour.
If the behaviour of agents within the entirety of `R` is of interest, then `B` can be set be a point of zero area.
The real experiment doesn't have to be exact in these terms, but the analysis will be more accurate as the experiment better approximates them.

Our analysis is focused on measuring the bias of agent movements towards the origin (the centre of `R`) as a function of the distance from the origin, `r`.

The analysis decomposes the movement of a population of agents into a directed-motion component, `f(r,t)`, and an undirected-diffusion component, `D(r,t)`.
In terms of agents performing a stochastic movement process, `f(r,t)` indicates the average velocity of agents towards the origin, while `D(r,t)` indicates the random diffusivity of agent movement.

Data for analysis is to consist of distances from the origin for agents within `R` for a set of time points. We assume time points occur at approximately regular intervals.
There is no requirement that the identity of the agent associated with each distance be known. I.e. agents do not need to be identified or tracked for analysis.
There is no requirement to track agents between time points.

We direct interested parties to the accompanying paper, but in essence observed agent positions are used to fit the advection and diffusion terms of a PDE using spatial statistics and functional regression.

The package is written in R.
This package is based on code/approaches outlined in the following text:

```
Ramsay JO, Hooker G, Graves S. Functional Data Analysis with R and MATLAB. 1st ed. Springer Publishing Company, Incorporated; 2009.
```

Other packages are required by `Chemotacticswarming`:

- For analysis:
`fda`,
`rainbow`,
`fields`,
`refund`.

- For simulations we used to create test data:
`plotrix`,
`Rfast`,
`useful`.

## Usage

Download a copy of this repository. This can be done directly with git, or by downloading a ZIP file (in which case, unzip it).

At present it assumes the repository is downloaded as a ZIP file, and that the folder title is changed to `Chemotacticswarming` and moved to `Desktop`. The directory paths within the code have been written with this assumption. 

However, you can put this wherever you want and change the main folder name, but you'll need to adjust the directory paths in the code we supply.

The entry points to our analyses are the functions captured in the modules `FunEstim_const.R` and `FunEstim_var.R`; these represent subtly different approaches to calculating diffusion and directed motion components of a population's motility under different assumptions.

We have included example data in this repository (those data from the accompanying paper), both from real world experiments and from simulations.

Experimental data sets include those found in `Data/Preembedcognate` and `Data/Noncognate`.
These are real experimental biological data generated from live cell imaging experiments of cytotoxic T cells interacting with tumour cells.
For details see the accompanying paper.

The `Data/RandWalk` datasets are simulated data used to verify the performance of the analysis for agent populations with known motility behaviour.
Code used to generate the simulated data and to generate graphs contrasting analysis results with known theoretical results are provided in the relevant folders.

The file `Analysis.r` captures how we run this analyses on each of these datasets.

### Input data format

A single comma-separated variable (csv) file is needed, named `distances.csv`. 
This captures distances of agents from the origin (centre of `R`) at corresponding time points. 
The `distances.csv` csv file must contain two columns: `distances` and `time`.

Each row represents an agent observed at a given point in time:

- `distance` is distance from the origin.
- `time` gives the corresponding time stamp.

Units for both columns are of your own choosing - the outputs and graphs generated by `Chemotacticswarming` are implicitly in these units (and units are not given in the graphs generated).

The analysis assumes an unchanging time interval between the times of population observation.

Other columns can be included (e.g. metadata), and are ignored by the analysis.

### Entry point to analysis

Two types of analyses of `f(r,t)` and `D(r,t)` are supplied in the following modules:

- `FunEstim_const.R`
- `FunEstim_var.R` 

Each analysis produces estimates for `f(r,t)` and `D(r,t)` under different assumptions detailed in the following sections.

Please see accompanying paper for specific details regarding the methodology used by these modules. 
In brief, these modules use smoothed focal K functions to calculate empirical estimates for `gamma(r,t)` (see accompanying paper fo definition), `lambda(r,t)` (agent density), and `lambda'(r,t)` (the derivative of lambda(r,t) with respect to `r`), and uses these as the response variable and covariates in a functional linear regression model to estimate `f(r,t)` and `D(r,t)`.

Each module calls `Prelim.R`, which performs preliminary preprocessing steps including producing the empirical focal K step functions, smoothing these functions using B-splines, and using these smoothed functions to produce empirical `lambda(r,t)` and `lambda'(r,t)` functional time series. Details regarding `Prelim.R` are provided in a following section.

### Usage tips

We note the following regarding the accuracy of estimates from these modules.

We find from investigating simulations studies, for which approximations of the true `f(r,t)` and `D(r,t)` functions are known, that estimates for `f(r,t)` and `D(r,t)` using the above modules are in general robust and accurate.

This is especially true for estimates for `f(r,t)`, which do not appear to be significantly influenced by changes in simulation parameters or any of the module parameters listed below related to fitting the functional linear models.

We note one limitation regarding these modules is that estimates for `D(r,t)` can be dependent on the smoothness of the focal K functions. 
This appears to be due to the fact that whilst the estimates for `lambda(r,t)`, the covariate for `f(r,t)`, are stable with respect to the smoothness of focal K functions, the magnitude of `lambda'(r,t)`, the covariate for `D(r,t)`, is dependent on the smoothness of focal K function.
This issue can lead to over or under estimation for `D(r,t)` depending on the smoothness of the focal K functions. 

We find that the using a large number of splines, without much smoothing, to capture the variability of the empirical focal K functions is in general preferable for producing accurate estimates for `D(r,t)`. A number of different splines should be compared with respect to the effect on fitting the empirical focal K functions and the estimates for `f(r,t)` and `D(r,t)`.

## Module `FunEstim_const.R`

This module estimates f and D assuming that both functions are time-invariant, i.e. f(r,t)=f(r) and D(r,t)=D(r).

### Inputs for module `FunEstim_const.R`

- `working_dir`; Location of data and where files/plots are to be saved.
- `tmin_iter`; Minimum time iteration to be used for analysis. If NULL, `tmin_iter` set to 1. Default is NULL. 
- `tmax_iter`; Maximum time iteration to be used for analysis. If NULL, `tmax_iter` set to maximum possible iteration. Default is NULL.
- `rR`; radius of the observed region, `R`.
- `rB`; radius of the region of interest, `B`.
- `r_resolution`; the resolution at which `r` is discretised for use in analysis. Default is 1.
- `normalise`; set to TRUE if using normalised focal K functions, FALSE if otherwise. Default is TRUE.
- `splines_K`; the number of B-splines used.
- `lambda_K`; the order of the smoothing parameter used to smooth K functions. If NULL, lambda_K is automatically selected by minimising a generalised cross-validation (GCV) measure. Default is NULL. 
- `min_loglam`; minimum value log(lambda_K) value used when searching for minimum GCV value. Default is -5.
- `max_loglam`; maximum value log(lambda_K) value used when searching for minimum GCV value.  Default is 10.
  - A warning message is produced if either `min_loglam` or `max_loglam` is selected in minimising GCV, which suggests that this range of values be extended.
- `Norder`; order of spline functions used to smooth empirical K step functions. Default is 6.
- `Lfdval`; derivative order to be penalised in smoothing. Default is 4.
  - Norder=6 and Lfdval=4 are set as defaults since we require the second derivative of the smoothed focal K functions in order to estimate both `lambda(r,t)` and `lambda'(r,t)`. Hence we smooth the empirical focal K functions using an order 6 B-spline basis such that the second derivative is order 4 (i.e. composed of a cubic B-spline basis) and penalise the fourth derivative such that the second derivative is smooth.
- `rmin_lambda_plot`; minimum `r` value used when plotting `lambda(r,t)` and `lambda'(r,t)`. If NULL, `rmin_lambda_plot` set to 0. Default is NULL.
- `rmax_lambda_plot`; maximum `r` value used when plotting `lambda(r,t)` and `lambda'(r,t)`. If NULL, `rmin_lambda_plot` set to `rR`. Default is NULL.
- `kymograph_zlim`; set range of z values for the heat map/kymograph of `lambda(r,t)`. Set to NULL if no specific limit set. Default is NULL.
- `rmin_est`; minimum r value used in estimating `f(r)` and `D(r)`. If NULL, `rmin_est` set to 0. Default is NULL.
- `rmax_est`; maximum r value used in estimating `f(r)` and `D(r)`. If NULL, `rmin_est` set to `rR`. Default is NULL.
- `rmin_plot`; minimum r value used in plotting `f(r)` and `D(r)`. If NULL, `rmin_plot` set to `rmin_est`. Default is NULL.
- `rmax_plot`; maximum r value used in plotting `f(r)` and `D(r)`. If NULL, `rmin_plot` set to `rmax_est`. Default is NULL.
- `splines_f`; set the number of splines used for estimating `f(r)`.
- `splines_D`; set the number of splines used for estimating `D(r)`.
- `lambda_f`; set the smoothing parameter used for estimating `f(r)`.
- `lambda_D`; set the smoothing parameter used for estimating `D(r)`.

### Outputs for module `FunEstim_const.R`

Note: We take `rmin_plot` and `rmax_plot` to define the `r` values of interest, and all csv files and plots are produced using only `r` within [`rmin_plot`, `rmax_plot`].

The following files are saved to the working directory, `working_dir`:

- `f.csv`; csv file for the estimate for f(r) and associated standard error.
- `D.csv`; csv file for the estimate for D(r) and associated standard error.
- `gamma.pdf`; a rainbow plot of empirical gamma(r,t) functions.
- `f.pdf`; a plot of the estimate for f(r).
- `D.pdf`; a plot of the estimate for D(r).

Returned values include:

* `gamma_fRegress`; fitted functional regression model.
* `rmin_plot`
* `rmax_plot`
* `r_trunc_values_plot`; `r` values within [`rmin_plot`,`rmax_plot`].
* `f_df`; data frame containing estimate for f(r), and associated standard error.
* `D_df`; data frame containing estimate for D(r), and associated standard error.

### Notes on choice of inputs for module `FunEstim_const.R`

`r_resolution` changes the resolution of `r` used in analysis and for plotting. This can be changed in order to account for the size of `rR`. The number of `r` points used should be (significantly) larger than the number of splines used to smooth empirical K functions. We find that estimates are robust to changes in `r_resolution`.

For this module, `FunEstim_const.R`, `normalise` is typically set to TRUE. 

In the accompanying paper focal K functions are in general normalised by (total number of agents)/(Volume of `R`) to ensure the theoretical focal K function associated with complete spatial randomness is pi*r^2, conforming with existing norms in spatial statistics regarding Ripley's K function. 
Normalisation also is required to produce meaningful focal L functions and the swarming metric `M`.

We allow for `normalise` to be FALSE to facilitate analysis of the case of `R` having an open boundary, with agents able to move across this boundary. See accompanying paper for details.

The parameters `splines_K` and `lambda_K` determine the smoothness of the smoothed focal K functions, and the resulting empirical gamma(r,t), lambda(r,t), and lambda'(r,t) functions that are then used to estimate f(r) and D(r). 

`kymograph_zlim` can be set by the user in order to produce heat maps/kymographs of lambda(r,t), which are useful if the user wishes to compare such plots for different data sets:

We find from analysis of our simulated data (for which approximate theoretical values of f(r) and D(r) are known against which we can compare our estimates) that taking `rmin_est` to be greater than 0 and `rmax_est` to be less than `rR` leads to more accurate estimates for f(r) and D(r) for the truncated `r` values. This is due to the fact that smoothed functions tend to have unstable fits at the start and end of the interval over which they are defined, and this instability is still greater for their derivatives. The functional linear regression uses covariates that are derived from the derivatives of the smoothed K functions. As such the estimates for f(r) and D(r) at or around r=0 and r=rR may be inaccurate.

In general, we find that making `rmin_est` slightly larger than 0, and `rmax_est` slightly smaller than `rR`, is enough to ensure robust accuracy of estimates.

Choices for `splines_f`, `lambda_f`, `splines_D`, and `lambda_D` should be made based on expectations regarding the smoothness of f(r) and D(r). This is dependent on the data set being analysed. 

## Module `FunEstim_var.R`

This module estimates f(r,t) and D(r,t) assuming dependance on both time and `r`.

As detailed in the accompanying paper this module estimates f(r,t) and D(r,t) by sequentially estimating them as if they were constant on a moving time window over the time series data. This involves setting the width of the moving window, specifically we set `h` to be half the window width, and then estimating f(r,t) and D(r,t) for a given time point, `t_{i}`, using time points within this moving window, i.e. all times within [t_{i-h},t_{i+h}].
Note that this approach can not estimate the first and last `h` time points for a given data set.

### Inputs for module FunEstim_var.R

Inputs for this module are identical to those for `FunEstim_const.R` with the following additional input:

- half_time_window_width; set `h`, half the width of the moving window. If NULL, `h` is set to floor(sqrt(number of time points)). Default is NULL.

### Outputs for module FunEstim_var.R

Note: We again take `rmin_plot` and `rmax_plot` to define the `r` values of interest, and all csv files and plots are produced using only for `r` within [`rmin_plot`,`rmax_plot`].
Furthermore, csv files and plots of f(r,t) and D(r,t) only contain data for the time indices for which estimates exist (see `window_centre_indices` below).

The following files are saved to the working directory, `working_dir`:

- `f_var.csv`; csv file for the estimate for f(r,t).
- `f_sde_var.csv`; csv file for the standard error for estimates f(r,t).
- `D_var.csv`; csv file for the estimate for D(r,t).
- `D_sde_var.csv`; csv file for the standard error for estimates D(r,t).
- `gamma.pdf`; a rainbow plot of empirical gamma(r,t) functions.
- `f_var.pdf`; a rainbow plot of the estimate for f(r,t).
- `D_var.pdf`; a rainbow plot of the estimate for D(r,t).

Returned values include:

* `gamma_fRegress`; fitted functional regression model.
* `rmin_plot`
* `rmax_plot`
* `r_trunc_values_plot`; `r` values within [`rmin_plot`,`rmax_plot`].
* `window_centre_indices`; time indices for which f(r,t) and D(r,t) are estimated.
* `f_fds`; a functional data object for estimates f(r,t).
* `f_var`; a matrix for estimates f(r,t).
* `f_sde_var`; a matrix for standard error for estimates f(r,t).
* `D_fds`; a functional data object for estimates D(r,t).
* `D_var`; a matrix for estimates D(r,t).
* `D_sde_var`; a matrix for standard error for estimates D(r,t).

### Guidance on setting input parameters

As briefly discussed in the accompanying paper, `half_time_window_width` should be chosen based on the expected underlying local variability of f(r,t) and D(r,t).
A rule of thumb when this is not known is to select half_time_window_width=floor(sqrt(N)), where N gives the number of time points for the data.
It is recommended that users try several values for `half_time_window_width` and that results are stable over these choices.

## Module `Prelim.R`; Preprocessing and preliminary analysis

`Prelim.R` performs essential preprocessing of the input data.
It calculates empirical focal K functions and other time series and functional time series used in the analysis.

It smooths the empirical focal K functions using B-splines. Parameters used for smoothing these functions are outlined in above sections.

Note that we find B-splines provide an easy to use approach to this analysis despite the fact that focal K functions are by definition monotonically increasing.

Of note, `Prelim.R` determines the 'swarming metric', `M`, a time series that gives a dynamic measure of agent aggregation, which we have found to be highly informative. For details regarding the swarming metric see accompanying paper.

Unless otherwise stated inputs listed below with are the same those described in the previous sections.

Functions included in `Prelim.R` include:

### Function `read_distance`

Converts csv file to data frame, truncated appropriately using user defined minimum and maximum time iterations to be included in analysis.

**Inputs**:

- `working_dir`; Location of data and where files/plots are to be saved.
- `tmin_iter`; Minimum time iteration to be used for analysis. If NULL, `tmin_iter` set to 1. Default is NULL. 
- `tmax_iter`; Maximum time iteration to be used for analysis. If NULL, `tmax_iter` set to maximum possible iteration. Default is NULL.

**Outputs**:

- `dist_data`; data frame of distances from the origin for each time point.
- `times_unique`; vector or individual time points.

### Function `Chemotax_preprocessing`

Uses `distances.csv` file to generating functional time series required for further analysis.

**Inputs**:

- `working_dir`; Location of data and where files/plots are to be saved.
- `dist_data`; data frame of distances from the origin for each time point.
- `times_unique`; vector or individual time points.
- `rR`; radius of the observed region, `R`.
- `rB`; radius of the region of interest, `B`.
- `r_resolution`; the resolution at which `r` is discretised for use in analysis. Default is 1.
- `normalise`; set to TRUE if using normalised focal K functions, FALSE if otherwise. Default is TRUE.
- `splines_K`; the number of B-splines used.
- `lambda_K`; the order of the smoothing parameter used to smooth K functions. If NULL, lambda_K is automatically selected by minimising the generalised cross-validation (GCV) measure. Default is NULL.
- `min_loglam`; minimum value log(lambda_K) value used when searching for minimum GCV value. Default is -5.
- `max_loglam`; maximum value log(lambda_K) value used when searching for minimum GCV value.  Default is 10.
- `Norder`; order of spline functions used to smooth empirical K step functions. Default is 6.
- `Lfdval`; derivative order to be penalised in smoothing. Default is 4.
- `rmin_lambda_plot`; minimum `r` value used when plotting lambda(r,t) and lambda'(r,t). Default is 0.
- `rmax_lambda_plot`; maximum `r` value used when plotting lambda(r,t) and lambda'(r,t). Default is `rR`.
- `kymograph_zlim`; set range of z values for the heat map/kymograph of lambda(r,t). Set to NULL if no specific limit set.
 
**Outputs**:

The following files are saved to the working directory, `working_dir`:
- `Kstep.csv`; csv file for empirical focal K functions.
- `M.csv`; csv file for swarming metric `M`.
- `M.pdf`; plot of swarming metric `M`.
- `L.pdf`; rainbow plot of empirical focal L functions.
- `L-r.pdf`; rainbow plot of empirical focal L functions - r + rB (can be easier to visualise).
- `K_step.pdf`; rainbow plot of empirical focal K functions.
- `n.pdf`; plot of number of agents within `R`.
- `nB.pdf`; plot of the number of agents within `B`.
- `K_smooth.pdf`; a rainbow plot of smoothed focal K functions.
- `L_smooth.pdf`; a rainbow plot of smoothed focal L functions.
- `L_smooth-r.pdf`; a rainbow plot of smoothed focal L functions - r + rB (can be easier to visualise).
- `K_res.pdf`; a rainbow plot of residuals between unsmoothed and smoothed focal K functions.
- `L_res.pdf`;  a rainbow plot of residuals between unsmoothed and smoothed focal L functions.
- `K_dr.pdf`; a rainbow plot of the first derivative of the smoothed focal K function.
- `lambda.pdf`; a rainbow plot of empirical estimates for lambda(r,t) for a range of `r` values determined by the user defined inputs `rmin_lambda_plot` and `rmax_lambda_plot`.
- `lambda_dr.pdf`; a rainbow plot of empirical estimates for lambda'(r,t) for a range of `r` values determined by the user defined inputs `rmin_lambda_plot` and `rmax_lambda_plot`.
- `lambda_surface.pdf`; a heat map/kymograph of lambda(r,t), an alternative means of plotting lambda(r,t) that we find easier to interpret than the rainbow plots. Again uses `rmin_lambda_plot` and `rmax_lambda_plot`.

Returned values include:

* `r_full`; vector of `r` values from 0 to `rR`, with step size determined by `r_resolution`.
* `n_time_points`; number of time points.
* `n_agents_R`; vector of number of agents in `R` at each time point.
* `n_agents_R_mat`; matrix of number of agents at each time point, repeats n_agents_R x (number of time points).
* `n_agents_B`; vector of number of agents within `B` at each time point.
* `R_vol`; volume of `R`.
* `Kstep_mat`; matrix of empirical focal K functions (step functions).
* `K_splinebasis`; B-spline basis used to smooth focal K functions.
* `Ksmooth`; smoothed focal K functions.
* `lambda`; a matrix of empirical estimates for lambda(r,t) functions.
* `lambda_dr`; a matrix of empirical estimates for lambda'(r,t) functions.
* `delta_times_vec`; vector of differences between time points.
* `delta_times_mat`; matrix of differences between time points.
* `Ksmooth_dt`; matrix of estimate of dK/dt for `r_full` values.

### Function `initialise_Kstep`

Called by Chemotax_preprocessing. 
Produces unsmoothed empirical focal K(r,t) step functions and associated output.

**Inputs**:

* `dist_data`; data frame of distances from the origin for each time point.
* `times_unique`; vector or individual time points.
* `rR`; radius of the observed region, `R`.
* `rB`; radius of the region of interest, `B`.
* `r_resolution`; the resolution at which `r` is discretised for use in analysis.
* `normalise`; set to TRUE if using normalised focal K functions, FALSE if otherwise.

**Outputs**:

* `r_full`; vector of `r` values from 0 to `rR`, with step size determined by `r_resolution`.
* `n_time_points`; number of time points.
* `Kstep_mat`; matrix of empirical focal K functions (step functions).
* `n_agents_R`; vector of number of agents in `R` at each time point.
* `n_agents_B`; vector of number of agents within `B` at each time point.
* `R_vol`; volume of `R`.

### Function `focal_L_M`
Converts unsmoothed, normalised, focal K step functions into focal L functions and swarming metric for each time point.

**Inputs**:

* `r_full`; vector of `r` values from 0 to `rR`, with step size determined by `r_resolution`.
* `Kstep_mat`; matrix of empirical focal K functions (step functions).
* `rR`; radius of the observed region, `R`.
* `rB`; radius of the region of interest, `B`.
* `n_agents_R`; vector of number of agents in `R` at each time point.
* `n_time_points`; number of time points.
* `R_vol`; volume of `R`.

**Outputs**:

* `L_mat`; matrix of empirical focal L functions (step functions).
* `M`; vector for swarming metric.

5. smooth_Kstep
Smoothes unsmoothed focal K functions using user supplied number of B-splines. 

**Inputs**:

* `r_full`; vector of `r` values from 0 to `rR`, with step size determined by `r_resolution`.
* `Kstep_mat`; matrix of empirical focal K functions (step functions).
* `splines_K`; the number of B-splines used.
* `lambda_K`; the order of the smoothing parameter used to smooth K functions.
* min_loglam; minimum value log(lambda_K) value used when searching for minimum GCV value. 
* max_loglam; maximum value log(lambda_K) value used when searching for minimum GCV value. 
* Norder; order of spline functions used to smooth empirical K step functions. 
* Lfdval; derivative order to be penalised in smoothing. 

**Outputs**:

- `Ksmooth`; smoothed focal K functions.
- `K_splinebasis`; B-spline basis for associated smoothed focal K functions.

### Function `plot_surface`

Function used to produce heatmap/kymograph of agent density, lambda(r,t).

**Inputs**:

- `x`; x values for heat map.
- `y`; y values for heat map.
- `z`; z values for heat map.
- `xlab`; label for x axis. Default is 'r'.
- `ylab`; label for y axis. Default is 'time'.
- `centred_val`; Value of `z` that is assigned white (e.g. use to indicate no biased movement w.r.t bolus). Default is NULL.
- `path`; Filename (including directories) to save graph to.  
- `xlim`; range of x values plotted. Default is range(x).
- `ylim`; range of y values plotted. Default is range(y).
- `zlim`; range of z values plotted. Default is range(z).
- `plot_contours`; TRUE/FALSE variable determining whether contour plots are added. Default is TRUE.

**Outputs**:

- ``path`.pdf; surface/heat map of lambda(r,t).

## Module `Analysis.R`

This file is used to call the previously outlined functions for the analysis of the provided example data sets.

To use:

- Move renamed `Chemotacticswarming` folder to Desktop. Alternatively, place folder elsewhere and change the following lines at the start of `Analysis.R` to ensure correct working directory:

  ```{r}
  setwd("~/Desktop/Chemotacticswarming/")
  data_location = "~/Desktop/Chemotacticswarming/Data/"
  ```

- Uncomment the `dataset` you wish to analyse. Ensure other `datasets` remain commented.

- If required, change input variables called by the various functions.

- Source `Analysis.R`

NOTE:
For simulations, for which we have known approximations for true f(r,t), and D(r,t) functions, `Analysis.R` calls additional R files to plot these known approximations against the estimates produced through functional linear regression.

## Example data

Example data sets are provided in the `Data` folder.
Simulations and experimental data are provided.

### Simulations

R files that produce simulations found in the accompanying paper are found in the `Data/RandWalk`.

Simulations consist of agent-based models with agents performing random walks under different chemotactic conditions.

Rules governing agent motility simulations available are those outlined in the accompanying paper.

Different types of simulations are made available and are found in the respective folders:

- `Unbiased`:
Simulates unbiased agent movement.

- `Attraction_const`:
Simulates time-invariant chemotaxis.

- `Attraction_var`:
Simulates time-varying chemotaxis.

Seeds can be set as indicated to obtain the simulations included in the accompanying publication.

Each of the above mentioned folders contains an R file (with an identical name) that produces output from a single simulation, a `distances.csv` csv file associated with a single simulation (for the given seed in the R script). If desired figures demonstrating agent positions for each time point can be generated by adding a folder '/X', to the given folder and uncommenting the relevant code in the simulation code.

### In vitro experiments

Experimental data from two in vitro live cell imaging experiments are provided. These are found in the `Data/Noncognate` and `Data/Preembededcognate` folders.

These experiments are a central component of a paper cited within the main manuscript that uncovers swarming behaviour in T cells.

In brief each experiment consists of cytotoxic T cells migrating through a matrix in vitro around a central collection of tumour cells, i.e. a "tumouroid". Tumour cells are either cognate or non-cognate. Upon contact with cognate tumour cells T cells become activated, whilst contact with non-cognate tumour cells does NOT cause T cells to become activated. Activated T cells produce a chemoattractant that induces chemotaxis in neighbouring T cells.

#### Non-Cognate experiment

The `Data/Noncognate` folder contains the `distances.csv` csv file associated with the experimental condition in which the central tumouroid contains non-cognate tumour cells. 

For both experimental data sets the `distances.csv` file is produced from cell position data such that the distances from the "origin", `r`, are the distances from the centre of geometric mass of the tumouroid.

#### Pre-Embedded-Cognate experiment

The `Data/Preembedcognate` folder contains the `distances.csv` csv file  assoiated with the experimental condition in which the central tumouroid contains a mixture of cognate tumour cells and pre-embedded T cells. 
