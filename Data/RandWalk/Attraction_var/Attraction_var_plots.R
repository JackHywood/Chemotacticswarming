# Plots estimated f and D functions against approximates for true f and D functions.
#
# Please see README.txt file for detailed explanation of the components of this script.
#
# Jack Hywood, Greg Rice, Mark N. Read, 2020

dir = "Attraction_var"

Del = 1 # Distance of each movement.
tau = 1 # Length of time step. Must not be an infinate decimal.
beta = 40 # Parameter for `strength of bias`. 
sig = 30 # Parameter for `range of bias`.

T_Full = 200 # Time of simulation.
times_unique = seq(0,T_Full,tau) # Vector of time points.
time_m = length(times_unique)
htwin= floor(sqrt(T_Full)) 

# Derivative dV/dr as defined in associated publicaion, which determines bias of agent movements:
dVdr = function(r){
  return(-beta*r/sig^2*exp(-(r^2)/(2*sig^2)))
}

# T(t) function:
Tfun = function(t){
    return(1+t*0.02)
  }

# Functions for c(r,t) and c_k(r,t) from associated publication:
cunb = function(r){
  integrand = function(theta){
    1/(2*pi)*(sqrt((r+Del*cos(theta))^2+(Del*sin(theta))^2)-r)^1
  }
  return(integrate(integrand,lower=0,upper=2*pi)$value)
}

c1rt = function(r,k){
  integrand = function(theta){
    1/(2*pi*besselI(abs(dVdr(r)*Tfun(k)),0))*exp(abs(dVdr(r)*Tfun(k))*cos(theta-pi))*((sqrt((r+Del*cos(theta))^2+(Del*sin(theta))^2)-r)^1-cunb(r))
  }
  return(integrate(integrand,lower=0,upper=2*pi)$value)
}

c2rt = function(r,k){
  integrand = function(theta){
    1/(2*pi*besselI(abs(dVdr(r)*Tfun(k)),0))*exp(abs(dVdr(r)*Tfun(k))*cos(theta-pi))*((sqrt((r+Del*cos(theta))^2+(Del*sin(theta))^2)-r)^1-cunb(r))^2
  }
  return(integrate(integrand,lower=0,upper=2*pi)$value)
}

c1 = matrix(0,length(r_trunc_values_plot),length((1+htwin):(time_m-1-htwin)))
c2 = matrix(0,length(r_trunc_values_plot),length((1+htwin):(time_m-1-htwin)))
for(i in 1:length(r_trunc_values_plot)){
  for(t2 in (1+htwin):(time_m-1-htwin)){
    c1[i,(t2-htwin)] = c1rt(r_trunc_values_plot[i],t2)
    c2[i,(t2-htwin)] = c2rt(r_trunc_values_plot[i],t2)
  }
}

# Theoretical f and D functions:
f_theor = c1/tau
D_theor = c2/(2*tau)

fds_ff_mat = fds(x=r_trunc_values_plot, y=f_theor)
fds_fdif_mat = fds(x=r_trunc_values_plot, y=f_var-f_theor)

fds_DD_mat = fds(x=r_trunc_values_plot, y=D_theor)
fds_Ddif_mat = fds(x=r_trunc_values_plot, y=D_var-D_theor)

pdf(paste(working_dir, '/f',dir,'.pdf', sep=''))
plot(f_fds,xlab="r", ylab="f",ylim=c(min(f_var,f_theor,0),max(f_var,f_theor,0)),
     xlim=c(rmin_plot,rmax_plot),
     las=1)
abline(h=0,lwd=2,lty=1,col="blue")
dev.off()

pdf(paste(working_dir, '/fexact',dir,'.pdf', sep=''))
plot(fds_ff_mat,xlab="r", ylab="f",ylim=c(min(f_var,f_theor,0),max(f_var,f_theor,0)),
     xlim=c(rmin_plot,rmax_plot),
     las=1)
abline(h=0,lwd=2,lty=1,col="blue")
dev.off()

pdf(paste(working_dir, '/ferror',dir,'.pdf', sep=''))
plot(fds_fdif_mat,xlab="r", ylab="Error",
     xlim=c(rmin_plot,rmax_plot),
     las=1)
dev.off()

pdf(paste(working_dir, '/D',dir,'.pdf', sep=''))
plot(D_fds,xlab="r", ylab="D",ylim=c(min(D_var,D_theor,0),max(D_var,D_theor,0)),
     xlim=c(rmin_plot,rmax_plot),
     las=1)
abline(h=0,lwd=2,lty=1,col="blue")
dev.off()

pdf(paste(working_dir, '/Dexact',dir,'.pdf', sep=''))
plot(fds_DD_mat,xlab="r", ylab="D",ylim=c(min(D_var,D_theor,0),max(D_var,D_theor,0)),
     xlim=c(rmin_plot,rmax_plot),
     las=1)
abline(h=0,lwd=2,lty=1,col="blue")
dev.off()

pdf(paste(working_dir, '/Derror',dir,'.pdf', sep=''))
plot(fds_Ddif_mat,xlab="r", ylab="Error",
     xlim=c(rmin_plot,rmax_plot),
     las=1)
dev.off()
