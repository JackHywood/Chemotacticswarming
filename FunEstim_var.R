# File produces time-varying estimates for f and D for given data set.
# 
# Please see README.txt file for detailed explanation of the components of this script.
#
# Jack Hywood, Mark N. Read, 2019

source('./Prelim.R')

FunEstim_var = function(
  working_dir,
  tmax_iter = NULL,
  tmin_iter = NULL,
  half_time_window_width = NULL,
  rR,
  rB,
  r_resolution = 1,
  normalise = TRUE,
  splines_K,
  lambda_K = NULL,
  min_loglam = -5,
  max_loglam = 10,
  Norder = 6,
  Lfdval = 4,
  rmin_lambda_plot = NULL,
  rmax_lambda_plot = NULL,
  kymograph_zlim = NULL,
  rmin_est = NULL,
  rmax_est = NULL,
  rmin_plot = NULL,
  rmax_plot = NULL,
  splines_f,
  splines_D,
  lambda_f,
  lambda_D
)
{
  if(is.null(rmin_lambda_plot)){
    rmin_lambda_plot = 0
  }
  if(is.null(rmax_lambda_plot)){
    rmax_lambda_plot = rR
  }
  if(is.null(rmin_est)){
    rmin_est = 0
  }
  if(is.null(rmax_est)){
    rmax_est = rR
  }
  if(is.null(rmin_plot)){
    rmin_plot = rmin_est
  }
  if(is.null(rmax_plot)){
    rmax_plot = rmax_est
  }
  
  # Read csv file and truncate time as defined by user:
  result = read_distance(working_dir=working_dir,
                         tmax_iter = tmax_iter,
                         tmin_iter = tmin_iter)
  dist_data         = result$dist_data
  times_unique      = result$times_unique
  
  # If h NULL, set h (half the moving window width used for analysis) to sqrt of number of time steps:
  if(is.null(half_time_window_width)){
    half_time_window_width = floor(sqrt(length(times_unique))) 
  }
  
  # Convert agent positions into estimates for K(r,t), lambda(r,t), lambda'(r,t), swarming metric, etc.:
  result = Chemotax_preprocessing(working_dir=working_dir,
                                  dist_data=dist_data, 
                                  times_unique=times_unique, 
                                  rB=rB, 
                                  rR=rR,
                                  r_resolution=r_resolution,
                                  normalise=normalise,
                                  splines_K=splines_K,
                                  lambda_K=lambda_K,
                                  min_loglam = min_loglam,
                                  max_loglam = max_loglam,
                                  Norder=Norder,
                                  Lfdval=Lfdval,
                                  rmin_lambda_plot=rmin_lambda_plot,
                                  rmax_lambda_plot=rmax_lambda_plot,
                                  kymograph_zlim=kymograph_zlim)
  r_full                 = result$r_full
  n_time_points          = result$n_time_points
  n_agents_R             = result$n_agents_R
  n_agents_R_mat         = result$n_agents_R_mat
  n_agents_B             = result$n_agents_B
  R_vol                  = result$R_vol
  Kstep_mat               = result$Kstep_mat
  K_splinebasis          = result$K_splinebasis
  Ksmooth                = result$Ksmooth
  lambda               = result$lambda
  lambda_dr              = result$lambda_dr
  delta_times_vec        = result$delta_times_vec
  delta_times_mat        = result$delta_times_mat
  Ksmooth_dt             = result$Ksmooth_dt
  
  # Define time indicies (i.e. time points) to be included for analysis: 
  window_centre_indices = (1 + half_time_window_width):(n_time_points - 1 - half_time_window_width)
  
  # Define gamma functions associated with unnormalised and normalised focal K functions respectively:
  gamma_full = - Ksmooth_dt / (2 * pi * (r_full))
  if (normalise) {
    gamma_full = gamma_full * (n_agents_R_mat[,-n_time_points] / R_vol)
  }
  
  # Set truncated r values to be used for estimation of f and D:
  r_usr_abs_indices = which(max(r_full[2],rmin_est) <= r_full & r_full <= rmax_est) 
  rmin_index = min(r_usr_abs_indices)
  rmax_index = max(r_usr_abs_indices)
  rmin = r_full[rmin_index]
  rmax = r_full[rmax_index]
  r_trunc_indices = rmin_index:rmax_index
  r_trunc_values = r_full[r_trunc_indices]
  
  # Define basis and functional parameters to be used to define response and covariate functional data to be used in regression:
  gamma_splinebasis_trc = create.bspline.basis(c(rmin, rmax), splines_K, norder=4)
  gamma_fdPar = fdPar(fdobj=gamma_splinebasis_trc, Lfdobj=NULL)
  
  # Convert gamma_full, lambda, lambda_dr to functional data objects on truncated r range:
  gamma_fd = smooth.basis(argvals=r_trunc_values,
                          y=gamma_full[r_trunc_indices,],
                          fdParobj=gamma_fdPar)
  lambda_fd = smooth.basis(argvals=r_trunc_values,
                             y=lambda[r_trunc_indices, -n_time_points],
                             fdParobj=gamma_fdPar)
  lambda_dr_fd = smooth.basis(argvals=r_trunc_values,
                              y=lambda_dr[r_trunc_indices, -n_time_points],
                              fdParobj=gamma_fdPar) 
  
  # Set functional parameters for estimated f and D functions to be used by fRegress function:
  f_splinebasis_trc = create.bspline.basis(range(r_trunc_values), splines_f)
  betafdPar_f = fdPar(f_splinebasis_trc, Lfdobj=NULL, lambda=lambda_f)
  D_splinebasis_trc = create.bspline.basis(range(r_trunc_values), splines_D)
  betafdPar_D = fdPar(D_splinebasis_trc, Lfdobj=NULL, lambda=lambda_D)
  
  # Define matrices for estimates f and D:
  emptymat = matrix(NA, length(r_full), length(window_centre_indices))
  f = emptymat
  f_sde = emptymat
  D = emptymat
  D_sde = emptymat
  index = 1
  for(window_centre in window_centre_indices) {
    # Time indices the window is currently operating over:
    window_coverage = (window_centre - half_time_window_width):(window_centre + half_time_window_width)
    print(paste("Window centred at time step ", window_centre, '; r values used = [', rmin, ', ', rmax, '].', sep=''))
    
    unused = tryCatch(
      {
        # Fit functional linear model for given window:
        gamma_fRegress = fRegress(gamma_fd$fd[window_coverage],
                                  xfdlist=list(lambda_fd$fd[window_coverage],
                                               lambda_dr_fd$fd[window_coverage]),
                                  betalist=list(betafdPar_f,
                                                betafdPar_D))
        gamma_betaestlist = gamma_fRegress$betaestlist 
        
        # Produce 95% confidence intervals:
        gamma_yhat_mat = eval.fd(r_trunc_values, gamma_fRegress$yhatfdobj$fd)
        gamma_rmatb   = gamma_full[r_trunc_indices, window_coverage] - gamma_yhat_mat
        gamma_sigma_eb = var(t(gamma_rmatb))
        gamma_y2cMap = gamma_fd$y2cMap
        gamma_betastderrlist = fRegress.stderr(gamma_fRegress, gamma_y2cMap, gamma_sigma_eb)
        
        f[r_trunc_indices, index] = eval.fd(r_trunc_values,gamma_betaestlist[[1]]$fd)
        f_sde[r_trunc_indices, index] = eval.fd(r_trunc_values, gamma_betastderrlist$betastderrlist[[1]])
        
        D[r_trunc_indices, index] = eval.fd(r_trunc_values, -gamma_betaestlist[[2]]$fd)
        D_sde[r_trunc_indices, index] = eval.fd(r_trunc_values, gamma_betastderrlist$betastderrlist[[2]])
        
      }, error = function(err) {
        print("Error generated whilst performing functional regression with window coverage and r values described above")
        print("This time series will not be included in the final analysis.")
      }
    )  # END tryCatch
    index = index + 1
  }
  
  # Set truncated r values to be used for plotting f and D:
  if(rmax_plot>rmax_est){
    r_usr_abs_indices_plot = which(rmin_plot <= r_full & r_full <= rmax_est)
    rmin_index_plot = min(r_usr_abs_indices_plot)
    rmax_index_plot = max(r_usr_abs_indices_plot)
    r_trunc_indices_plot = rmin_index_plot:rmax_index_plot
    r_trunc_values_plot = r_full[r_trunc_indices_plot]
  }else{
    r_usr_abs_indices_plot = which(rmin_plot <= r_full & r_full <= rmax_plot)
    rmin_index_plot = min(r_usr_abs_indices_plot)
    rmax_index_plot = max(r_usr_abs_indices_plot)
    r_trunc_indices_plot = rmin_index_plot:rmax_index_plot
    r_trunc_values_plot = r_full[r_trunc_indices_plot]
  }
  
  # Define data frame of f estimate, and save as csv file:
  f_var = f[r_trunc_indices_plot,]
  write.csv(f_var, file=paste(working_dir, "/f_var", ".csv", sep=''), row.names=FALSE)
  f_sde_var = f_sde[r_trunc_indices_plot,]
  write.csv(f_sde_var, file=paste(working_dir, "/f_sde_var", ".csv", sep=''), row.names=FALSE)
  
  # Define data frame of D estimate, and save as csv file:
  D_var = D[r_trunc_indices_plot,]
  write.csv(D_var, file=paste(working_dir, "/D_var", ".csv", sep=''), row.names=FALSE)
  D_sde_var = D_sde[r_trunc_indices_plot,]
  write.csv(D_sde_var, file=paste(working_dir, "/D_sde_var", ".csv", sep=''), row.names=FALSE)
  
  # Plot gamma, and estimates for f and D:
  
  fds_gamma_fd = fds(x=r_trunc_values_plot, y=gamma_full[r_trunc_indices_plot,])
  f_fds = fds(x=r_trunc_values_plot, y=f_var)
  D_fds = fds(x=r_trunc_values_plot, y=D_var)
  
  pdf(paste(working_dir, '/gamma.pdf', sep=''))
  plot(fds_gamma_fd, xlab="r", ylab="gamma",xlim=c(rmin_plot,rmax_plot),
       las=1)
  dev.off()
  
  pdf(paste(working_dir, '/f_var.pdf', sep=''))
  plot(f_fds, xlab="r", ylab="f", 
       xlim=c(rmin_plot,rmax_plot),
       ylim=c(min(f_var), 
              max(f_var,0)),
       las=1)
  abline(h=0,lwd=2,lty=1,col="blue")
  dev.off()
  
  pdf(paste(working_dir, '/D_var.pdf', sep=''))
  plot(D_fds, xlab="r", ylab="D",
       xlim=c(rmin_plot,rmax_plot),
       ylim=c(min(D_var,0), 
              max(D_var)),
       las=1)
  abline(h=0,lwd=2,lty=1,col="blue")
  dev.off()
  
  result = list(gamma_fRegress = gamma_fRegress,
                rmin_plot = rmin_plot,
                rmax_plot = rmax_plot,
                r_trunc_values_plot = r_trunc_values_plot,
                window_centre_indices = window_centre_indices,
                f_fds = f_fds,
                f_var = f_var,
                f_sde_var = f_sde_var,
                D_fds = D_fds,
                D_var = D_var,
                D_sde_var = D_sde_var
  )
  return(result)
  
}
