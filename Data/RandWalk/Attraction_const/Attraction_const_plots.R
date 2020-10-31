# Plots estimated f and D functions against approximates for true f and D functions.
#
# Please see README.txt file for detailed explanation of the components of this script.
#
# Jack Hywood, Mark N. Read, 2019

dir = "Attraction_const"

Del = 1 # Distance of each movement.
tau = 1 # Length of time step.
beta = 40 # Parameter for strength of bias. 
sig = 30 # Parameter for range of bias.

# Derivative dV/dr as defined in associated publicaion, which determines bias of agent movements:
dVdr = function(r){
  return(-beta*r/sig^2*exp(-(r^2)/(2*sig^2)))
}

# Functions for c(r,t) and c_k(r,t) from associated publication:
cunb = function(r){
  integrand = function(theta){
    1/(2*pi)*(sqrt((r+Del*cos(theta))^2+(Del*sin(theta))^2)-r)^1
  }
  return(integrate(integrand,lower=0,upper=2*pi)$value)
}

c1r = function(r){
  integrand = function(theta){
    1/(2*pi*besselI(abs(dVdr(r)),0))*exp(abs(dVdr(r))*cos(theta-pi))*((sqrt((r+Del*cos(theta))^2+(Del*sin(theta))^2)-r)^1-cunb(r))
  }
  return(integrate(integrand,lower=0,upper=2*pi)$value)
}

c2r = function(r){
  integrand = function(theta){
    1/(2*pi*besselI(abs(dVdr(r)),0))*exp(abs(dVdr(r))*cos(theta-pi))*((sqrt((r+Del*cos(theta))^2+(Del*sin(theta))^2)-r)^1-cunb(r))^2
  }
  return(integrate(integrand,lower=0,upper=2*pi)$value)
}

c1 = rep(0,length(r_trunc_values_plot))
c2 = rep(0,length(r_trunc_values_plot))
for(i in 1:length(r_trunc_values_plot)){
  c1[i] = c1r(r_trunc_values_plot[i])
  c2[i] = c2r(r_trunc_values_plot[i])
}

# Theoretical f and D functions:
f_theor = c1/tau
D_theor = c2/(2*tau)

pdf(paste(working_dir, '/f',dir,'.pdf', sep=''))
plot(r_trunc_values_plot, f_df$f, type="l", lwd=3,lty=1,
     ylim=c(-0.5,0),
     xlim=c(rmin_plot,rmax_plot),
     xlab="r", ylab="f",
     las=1)
lines(r_trunc_values_plot, f_df$f-2*f_df$f_sde, lwd=2, lty=3)
lines(r_trunc_values_plot, f_df$f+2*f_df$f_sde, lwd=2, lty=3)
abline(h=0, lwd=2, lty=1, col="blue")
lines(r_trunc_values_plot,f_theor,lwd=2,col="red")
dev.off()

pdf(paste(working_dir, '/D',dir,'.pdf', sep=''))
plot(r_trunc_values_plot, D_df$D, type="l", lwd=3,lty=1,
     ylim=c(0, 0.5),
     xlim=c(rmin_plot,rmax_plot),
     xlab="r", ylab="D",
     las=1)
lines(r_trunc_values_plot, D_df$D-2*D_df$D_sde, lwd=2, lty=3)
lines(r_trunc_values_plot, D_df$D+2*D_df$D_sde, lwd=2, lty=3)
abline(h=0,lwd=2,lty=1,col="blue")
lines(r_trunc_values_plot,D_theor,lwd=2,col="red")
dev.off()
