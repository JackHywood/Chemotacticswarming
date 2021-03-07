# `Prelim.R` contains functions that convert csv file with agent distances from origin (i.e. r values) into functional time series for preliminary analysis and for functional linear regression.
#
# Please see README.txt file for detailed explanation of the components of this script.
#
# Jack Hywood, Greg Rice, Mark N. Read, 2020

read_distance = function(
  working_dir,
  tmin_iter = NULL,
  tmax_iter = NULL
)
{
  # Read and define the "distance" files:
  dist_file = paste(working_dir, "/distances.csv", sep='')
  if (!file.exists(dist_file)) {
    print(dist_file)
    stop('specified distance to bolus file does not exist')
  }
  
  # Table must contain two columns, distance and then times at which those distances were observed:
  dist_data = read.csv(dist_file,
                       skip=0,
                       colClasses=c(NA,NA))
  dist_data = dist_data[order(dist_data$time), ]  # Sort by monotonically increasing time.
  dist_data$time = dist_data$time
  
  # Truncate time as required by user:
  times_unique = unique(dist_data$time)
  if (is.null(tmin_iter))
  {tmin_iter = 1}
  if (is.null(tmax_iter))
  {tmax_iter = length(times_unique)}
  
  dist_data = dist_data[which((dist_data$time>=times_unique[tmin_iter])&(dist_data$time<=times_unique[tmax_iter])),]
  times_unique = unique(dist_data$time)
  
  result = list(dist_data=dist_data,
                times_unique=times_unique)
  return(result)
}

Chemotax_preprocessing = function(
  working_dir,
  dist_data,
  times_unique,
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
  kymograph_zlim = NULL
) 
{
  if(is.null(rmin_lambda_plot)){
    rmin_lambda_plot = 0
  }
  if(is.null(rmax_lambda_plot)){
    rmax_lambda_plot = rR
  }
  
  result = initialise_Kstep(dist_data=dist_data,
                            times_unique=times_unique, 
                            rR=rR, 
                            rB=rB, 
                            r_resolution=r_resolution,
                            normalise=normalise)
  r_full          = result$r_full
  n_time_points   = result$n_time_points
  Kstep_mat       = result$Kstep_mat
  n_agents_R      = result$n_agents_R
  n_agents_B      = result$n_agents_B
  R_vol           = result$R_vol
  
  # Write csv file for K functions:
  write.csv(Kstep_mat, file=paste(working_dir, "/Kstep",".csv", sep=''))
  
  # Write csv file for M metric and plot L functions and M metric:
  if (normalise) {
    print("Calculating M and L")
    result = focal_L_M(r_full=r_full, 
                       Kstep_mat=Kstep_mat, 
                       rR=rR, 
                       rB=rB,
                       n_agents_R=n_agents_R, 
                       n_time_points=n_time_points,
                       R_vol=R_vol)
    L_mat = result$L_mat
    M     = result$M
    
    M_df = data.frame(times_unique, M)
    names(M_df)[1] = "time_step"
    write.csv(M_df, file=paste(working_dir, "/M", ".csv", sep=''), row.names=FALSE)
    
    pdf(paste(working_dir, '/M.pdf', sep=''))
    plot(times_unique, M, 
         type="l", 
         lwd=2, 
         ylim=c(-1,1), 
         xlab="Time", ylab="M",
         las=1)
    dev.off()
    
    fds_L_mat = fds(x=r_full[which(r_full>=rB)], y=t(L_mat))  # (step) L function values
    pdf(paste(working_dir, '/L.pdf', sep=''))
    plot(fds_L_mat,
         xlab="r", 
         ylab="L", 
         xlim=c(rB,rR), 
         ylim=c(0, max(L_mat)),
         las=1)
    dev.off()
    
    fds_L_mat_r = fds(x=r_full[which(r_full>=rB)], y=t(L_mat)-r_full[which(r_full>=rB)]+rB)  # (step) L function values
    pdf(paste(working_dir, '/L-r.pdf', sep=''))
    plot(fds_L_mat_r,
         xlab="r", 
         ylab="L-r", 
         xlim=c(rB,rR),
         las=1)
    dev.off()
    
  } else {
    print("Not calculating M and L")
  }
  
  # Generate smoothed K functions
  result = smooth_Kstep(r_full=r_full,
                        Kstep_mat=Kstep_mat, 
                        splines_K=splines_K, 
                        lambda_K=lambda_K,
                        min_loglam=min_loglam, 
                        max_loglam=max_loglam, 
                        Norder=Norder,
                        Lfdval=Lfdval)
  

    Ksmooth        = result$Ksmooth
    K_splinebasis  = result$K_splinebasis
    
    # First and second derivatives of K (smoothed) with respect to r.
    Ksm_dr = deriv.fd(Ksmooth$fd, 1)
    Ksm_dr2 = deriv.fd(Ksmooth$fd, 2)
    
    # Find differences in smoothed K between r and r+1, for all time points (and r). Has shape [r-1, n_time_points-1].
    delta_times_vec = times_unique[-1] - times_unique[-n_time_points]
    delta_times_mat = t(matrix(rep(delta_times_vec, length(r_full)), nrow=n_time_points-1))
    Ksmooth_dt = (eval.fd(r_full, Ksmooth$fd)[,-1]
                  - eval.fd(r_full, Ksmooth$fd)[, -n_time_points]
    ) / delta_times_mat
    
    # Repeats n_agents_R x (number of time points)
    n_agents_R_mat = t(matrix(rep(n_agents_R, length(r_full)), nrow=n_time_points))
    
    # Estimates for lambda(r,t) and its first derivative w.r.t r
    if (normalise) {
      lambda = eval.fd(r_full, Ksm_dr) * n_agents_R_mat / ( (r_full) * 2 * pi * R_vol )
      lambda_dr = ( (eval.fd(r_full, Ksm_dr2) * (r_full)) - eval.fd(r_full, Ksm_dr) )  *
        ( n_agents_R_mat / ( (r_full)^2 * 2 * pi * R_vol ) )
    } else {
      lambda = eval.fd(r_full, Ksm_dr) / ( (r_full) * 2 * pi )
      lambda_dr = ( (eval.fd(r_full, Ksm_dr2) * (r_full)) - eval.fd(r_full, Ksm_dr) )  *
        ( 1 / ( (r_full)^2 * 2 * pi ) )
    }
  
  # "fds" variables arrange data for rainbow plots. 
  fds_Kstep = fds(x=r_full, y=t(Kstep_mat))
  pdf(paste(working_dir, '/K_step.pdf', sep=''))
  plot(fds_Kstep, xlab="r", ylab="K (step)", ylim=c(0, max(Kstep_mat)), xlim=range(r_full),las=1)
  dev.off()
  
  pdf(paste(working_dir, '/n.pdf', sep=''))
  plot(times_unique,n_agents_R,ylim=c(0,max(n_agents_R)),type="l", xlab="Time", ylab="Agents",las=1)
  dev.off()
  
  pdf(paste(working_dir, '/nB.pdf', sep=''))
  plot(times_unique,n_agents_B,ylim=c(0,max(n_agents_B)),type="l", xlab="Time", ylab="Agents",las=1)
  dev.off()
  
  fds_Ksm = fds(x=r_full, y=eval.fd(r_full, Ksmooth$fd))  # Smoothed K function values
  pdf(paste(working_dir, '/K_smooth.pdf', sep=''))
  plot(fds_Ksm, xlab="r", ylab="K (smoothed)", ylim=c(min(Kstep_mat),max(Kstep_mat)), xlim=range(r_full),las=1)
  dev.off()
  
  smooth_result = focal_L_M(r_full=r_full, 
                     Kstep_mat= t(eval.fd(r_full, Ksmooth$fd)), 
                     rR=rR, 
                     rB=rB,
                     n_agents_R=n_agents_R, 
                     n_time_points=n_time_points,
                     R_vol=R_vol)
  L_smooth_mat = smooth_result$L_mat
  
  fds_L_smooth = fds(x=r_full[which(r_full>=rB)], y=t(L_smooth_mat))
  pdf(paste(working_dir, '/L_smooth.pdf', sep=''))
  plot(fds_L_smooth, xlab="r", ylab="L (smooth)", ylim=c(0, max(L_smooth_mat)), xlim=range(r_full[which(r_full>=rB)]),las=1)
  dev.off()
  
  fds_L_smooth_r = fds(x=r_full[which(r_full>=rB)], y=t(L_smooth_mat)-r_full[which(r_full>=rB)]+rB)
  pdf(paste(working_dir, '/L_smooth-r.pdf', sep=''))
  plot(fds_L_smooth_r, xlab="r", ylab="L (smooth) - r",xlim=range(r_full[which(r_full>=rB)]),las=1)
  dev.off()
  
  Kres = t(Kstep_mat) - eval.fd(r_full, Ksmooth$fd)
  fds_Kres = fds(x=r_full, y=Kres)
  pdf(paste(working_dir, '/Kres.pdf', sep=''))
  plot(fds_Kres, xlab="r", ylab="K Difference (smoothed)", xlim=range(r_full))
  dev.off()
  
  Lres = t(L_mat) - t(L_smooth_mat)
  fds_Lres = fds(x=r_full[which(r_full>=rB)], y=Lres)
  pdf(paste(working_dir, '/Lres.pdf', sep=''))
  plot(fds_Lres, xlab="r", ylab="L Difference (smoothed)", xlim=range(r_full[which(r_full>=rB)]))
  dev.off()
  
  fds_Kdr = fds(x=r_full, y=eval.fd(r_full, Ksm_dr))
  pdf(paste(working_dir, '/K_dr.pdf', sep=''))
  plot(fds_Kdr, xlab="r", ylab="K'", xlim=range(r_full),las=1)
  dev.off()
  
  fds_lambda_plot = fds(x=r_full[which(r_full>=rmin_lambda_plot & r_full<= rmax_lambda_plot)], 
                        y=lambda[which(r_full>=rmin_lambda_plot & r_full<= rmax_lambda_plot),])
  pdf(paste(working_dir, '/lambda.pdf', sep=''))
  plot(fds_lambda_plot, xlab="r", ylab="lambda",las=1)
  dev.off()
  
  fds_lambda_dr_plot = fds(x=r_full[which(r_full>=rmin_lambda_plot & r_full<= rmax_lambda_plot)], 
                           y=lambda_dr[which(r_full>=rmin_lambda_plot & r_full<= rmax_lambda_plot),])
  pdf(paste(working_dir, '/lambda_dr.pdf', sep=''))
  plot(fds_lambda_dr_plot, xlab="r", ylab="lambda'",las=1)
  dev.off()
  
  # Plot lambda(r,t) as surface:
  plot_surface(x=r_full[which(r_full>=rmin_lambda_plot & r_full<= rmax_lambda_plot)], 
               y=times_unique, 
               z=lambda[which(r_full>=rmin_lambda_plot & r_full<= rmax_lambda_plot),],
               path=paste(working_dir, '/lambda_surface.pdf', sep=''),
               zlim=kymograph_zlim, 
               plot_contours=FALSE)
  
  result = list(r_full=r_full,
                n_time_points=n_time_points,
                n_agents_R=n_agents_R,
                n_agents_R_mat=n_agents_R_mat,
                n_agents_B=n_agents_B,
                R_vol=R_vol,
                Kstep_mat=Kstep_mat,
                K_splinebasis=K_splinebasis,
                Ksmooth=Ksmooth,
                lambda=lambda,
                lambda_dr=lambda_dr,
                delta_times_vec=delta_times_vec,
                delta_times_mat=delta_times_mat,
                Ksmooth_dt=Ksmooth_dt
  )
  print("Completed common pre-processing stage.")
  return(result)
}


initialise_Kstep = function(dist_data, 
                            times_unique, 
                            rR, 
                            rB, 
                            r_resolution, 
                            normalise)
{
  n_time_points = length(times_unique)  # Number of time steps in data set
  r_full = seq(0, rR, length.out=floor(rR / r_resolution) + 1) # r values used to produce empirical K functions
  R_vol = pi * rR^2  # Volume of region R; 2D.
  Kstep_mat = matrix(0, n_time_points, length(r_full)) # Unsmoothed K functions
  n_agents_B = rep(0, n_time_points)  # Number of agents in B, i.e. inside bolus, at time step j
  n_agents_R = rep(0, n_time_points)  # Number of agents in R at time step j
  
  # Populate Kstep_mat, iterate over time steps. 
  for (j in 1:n_time_points) 
  {
    time_sample = times_unique[j]
    # Collect cell distances, for cells observed at time j
    distances_j = dist_data$distance[which(dist_data$time == time_sample)]
    
    n_agents_B[j] = sum(distances_j <= rB)  
    n_agents_R[j] = sum(distances_j <= rR) 
    
    # Produce a vector counting the number of agents within "r" at current time step.
    k = sapply(r_full, function(r) sum(distances_j <= r))  # Focal K function for this point in time. 
    if(normalise) {
      # normalise K functions
      k = k * R_vol / n_agents_R[j]
    }
    Kstep_mat[j, ] = k
    print(paste('Building focal K functions, processing time index ', j, sep=''))
  }
  
  result = list(r_full=r_full, 
                n_time_points=n_time_points,
                Kstep_mat=Kstep_mat,
                n_agents_R=n_agents_R,  
                n_agents_B=n_agents_B,
                R_vol=R_vol)
  return(result)
}

focal_L_M = function(r_full,
                     Kstep_mat, 
                     rR,
                     rB, 
                     n_agents_R, 
                     n_time_points,
                     R_vol)
{
  # Produce L functions:
  L_mat = sqrt((Kstep_mat[,which(rB <= r_full)] / R_vol*(pi*(rR^2-rB^2)) + pi * rB^2) /pi) - rB
  
  # Produce metric M:
  dr = r_full[2] - r_full[1]  # Difference in r, for integration below.
  M = (2 / (rR - rB)^2) *  
    (rowSums(t(  
      t(L_mat) - dr*(0:(length(which(rB <= r_full))-1)))  
    ) * dr  
    )
  
  result = list(L_mat=L_mat, 
                M=M)
  return(result)
}

smooth_Kstep = function(r_full,
                        Kstep_mat, 
                        splines_K,
                        lambda_K, 
                        min_loglam, 
                        max_loglam, 
                        Norder,
                        Lfdval)
{
  K_splinebasis = create.bspline.basis(range(r_full), splines_K, norder=Norder)
    # Choose the the optimal smoothing parameter lambda_K by minimising GCV.
    if (is.null(lambda_K)) {
      print("User has not selected the smoothing penalty to use in smoothing (lambda_K).")
      # Gives ranges of log(lambda_K) values searched. 
      # May need to be adjusted for individual data sets.
      # Warning messages are generated if the max/min value specified by user is selected; the optimal value may lie beyond 
      # the range explored. 
      loglam = seq(min_loglam, max_loglam, 0.25)
      nlam = length(loglam)
      gcvsave = rep(NA, nlam)
      for (ilam in 1:nlam) {  # cycle through each (log) lambda_K value
        lambda_K_i = 10^loglam[ilam]
        fdParobj = fdPar(K_splinebasis, Lfdobj=Lfdval, lambda=lambda_K_i)
        smoothlist = smooth.basis(argvals=r_full, y=t(Kstep_mat), fdParobj)
        gcvsave[ilam] = sum(smoothlist$gcv)
      }
      
      pdf(paste(working_dir, '/smoothing_parameter.pdf', sep=''))
      plot(loglam, gcvsave, xlab="log(lambda)", ylab="GCV(lambda)")
      dev.off()
      
      optimal_loglam = loglam[which.min(gcvsave)]
      print(paste('Optimal log10 lambda =', optimal_loglam))
      lambda_K = 10^(optimal_loglam)  # The optimal lambda_K value
      if (min_loglam == optimal_loglam) { 
        print(paste("WARNING! GCV loglam procedure selected the minimum specified possibility, 10^", min_loglam,
                    "\n  You may wish to select a LOWER min_loglam value to ensure there isn't a better possibility", sep=''))
      }
      if (max_loglam == optimal_loglam) { 
        print(paste("WARNING! GCV loglam procedure selected the maximum specified possibility, 10^", min_loglam,
                    "\n  You may wish to select a HIGHER min_loglam value to ensure there isn't a better possibility", sep=''))
      }
    }
    print(paste('Smoothing parameter =', round(lambda_K,5)))
    # Using optimal lambda_K as established above.
    fdParobj = fdPar(K_splinebasis, Lfdobj=Lfdval, lambda=lambda_K)
    Ksmooth = smooth.basis(argvals=r_full, y=t(Kstep_mat), fdParobj)# Smoothed K function values
    
  result = list(Ksmooth=Ksmooth, 
                K_splinebasis=K_splinebasis)
  return(result)
}

plot_surface = function(
  x, y, z, 
  xlab='r', 
  ylab='time', 
  centred_val=NULL, # Value of z that is assigned white (e.g. use to indicate no biased movement w.r.t bolus)
  path, # Filename (including directories) to save graph to. 
  xlim=range(x), ylim=range(y), zlim=range(z),  # Set to c(lower, upper) to specify limits. 
  plot_contours=TRUE
)
{
  min_val = min(z, na.rm=TRUE)
  max_val = max(z, na.rm=TRUE)
  colmap = NULL
  if (! is.null(centred_val)) {
    colmap = designer.colors(250, c("blue","white", "red"))
  }
  pdf(path)
  # Automatically includes a colour bar
  image.plot(x=x, y=y, z=z, xlab=xlab,  ylab=ylab, col=colmap,
             xlim=xlim, ylim=ylim, zlim=zlim)
  
  if (plot_contours) {
    contour(x=x, y=y, z=z, add=TRUE, nlevels=1,
            levels=pretty(range(min_val, max_val), 8), col='dark grey')
    
    # Heavier weight for contour z=centred_val, if selected.
    if (! is.null(centred_val)) {
      contour(x=x, y=y, z=z, add=TRUE, nlevels=1, levels=c(centred_val), lwd=1.0)
    }
  }
  dev.off()
}

