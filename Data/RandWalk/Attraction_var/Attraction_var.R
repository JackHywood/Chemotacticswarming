# Simuation of individual based model with increasing attraction.
# Simulates a lattice-free random walk model without space occupying properties.
# Rules governing agent movement as described in associated publication.
#
# Please see README.txt file for detailed explanation of the components of this script.
#
# Jack Hywood, Mark N. Read, 2019

require(plotrix)
require(Rfast)
require(useful)

setwd("~/Desktop/Chemotacticswarming-master/Data/RandWalk/Attraction_var/")
working_dir = "~/Desktop/Chemotacticswarming-master/Data/RandWalk/Attraction_var/"
dir = "Attraction_var"

SEED = 1 # Set to NULL if no seed used. Set to 1 to generate the simulation from the associated publication.

N = 1000 # Number of agents uniformly distributed in region R.
B.empty = FALSE # Set B.empty = FALSE to keep all agents uniformly distributed in R. Set B.empty = TRUE to remove all agents within B.

Del = 1 # Distance of each movement.
tau = 1 # Duration of each time step.
T_Full = 200 # Duration of simulation.
times_unique = seq(0,T_Full,tau) # Vector of time points.
time_m = length(times_unique) # Number of time points

rR = 100 # Radius of R.
rB = 20 # Radius of B.

beta = 40 # Parameter for `strength of bias`. 
sig = 30 # Parameter for `range of bias`.

# Derivative dV/dr as defined in associated publicaion, which determines bias of agent movements:
dVdr = function(z){
  return(-beta*z/sig^2*exp(-(z^2)/(2*sig^2)))
}

# T(t) function:
Tfun = function(t){
  return(1+t*0.02)
}

if(is.null(SEED)==FALSE){
  if(SEED==1){
    set.seed(101)
  }
}

# Populate R uniformly with agents:
rhoX = rR*sqrt(runif(N, min = 0, max = 1))
thetaX = runif(N, min = 0, max = 2*pi)
X.x = rhoX*cos(thetaX)
X.y = rhoX*sin(thetaX)
X = cbind(X.x,X.y)

# If B.empty==TRUE then agents within B are removed:
if(B.empty==TRUE && sum(((X[,1])^2+(X[,2])^2)<=rB^2)>0){
  X = X[-which(((X[,1])^2+(X[,2])^2)<=rB^2),]
  N = dim(X)[1]
}

# Plot initial agent positions:
# pdf(paste(working_dir, '/X', '/X',0,'.pdf', sep=''))
# plot(X,xlim=c(-rR,rR), ylim=c(-rR,rR),asp=1,pch=20,
#      ylab="",
#      xlab="",
#      xaxt="n", yaxt="n",axes=FALSE,
#      cex.axis=1.5,
#      cex.lab=2,
#      las=1)
# draw.circle(0,0,rB,nv=100,border=NULL,col="black",lty=1,density=NULL,
#             angle=45,lwd=1)
# draw.circle(0,0,rR,nv=100,border=NULL,lty=1,density=NULL,
#             angle=45,lwd=1)
# points(X[which(((X[,1])^2+(X[,2])^2)<=rB^2),],
#        bg="red",col="red",pch=20)
# dev.off()


X.mat = array(0,dim=c(time_m,N,2))
X.mat[1,,] = X
for(k in 1:(time_m-1)){
  mv.agents = sample(1:N, N)
  for(l in mv.agents){
    r = sqrt(X[l,1]^2+X[l,2]^2) # r value
    dv = dVdr(r)*Tfun(k) # Gradient of v(r,t) 
    mu = cart2pol(-X[l,1],-X[l,2])[,2][[1]] # Mean direction in radians of direction moved for each agent
    
    # Generate direction moved via von Mises distribution:
    mv.theta = rvonmises(n=1,m = mu, k = abs(dv))
    mv = Del*c(cos(mv.theta),sin(mv.theta)) # Vector of movement (prior to reflection)
    X.mv = X[l,] + mv  # New agent positions prior to being corrected for reflections off R boundary
    
    # The following reflects agents that have moved out of R back into R:
    
    if((X.mv[1]^2+X.mv[2]^2)>rR^2){
      dx = X.mv[1]-X[l,1]
      dy = X.mv[2]-X[l,2]
      dr = sqrt(dx^2+dy^2)
      D = X[l,1]*X.mv[2] - X.mv[1]*X[l,2]
      xint1 = (D*dy+sign(dy)*dx*sqrt(rR^2*dr^2-D^2))/dr^2
      xint2 = (D*dy-sign(dy)*dx*sqrt(rR^2*dr^2-D^2))/dr^2
      yint1 = (-D*dx+abs(dy)*sqrt(rR^2*dr^2-D^2))/dr^2
      yint2 = (-D*dx-abs(dy)*sqrt(rR^2*dr^2-D^2))/dr^2
      d1=sqrt((xint1-X.mv[1])^2+(yint1-X.mv[2])^2)
      d2=sqrt((xint2-X.mv[1])^2+(yint2-X.mv[2])^2)
      if(d1<d2){
        P=c(xint1,yint1)
      }else{
        P=c(xint2,yint2)
      }
      a=c(X.mv[1]-P[1],X.mv[2]-P[2])
      X.mv = X.mv-2*(P%*%a/rR^2)%*%P 
    }
    X[l,]=X.mv
  }
  X.mat[k+1,,] = X
  
  # Plot agents positions:
  # pdf(paste(working_dir, '/X', '/X',k,'.pdf', sep=''))
  # plot(X,xlim=c(-rR,rR), ylim=c(-rR,rR),asp=1,pch=20,
  #      xaxt="n", yaxt="n",axes=FALSE,
  #      ylab="",
  #      xlab="",
  #      cex.axis=1.5,
  #      cex.lab=2,
  #      las=1)
  # draw.circle(0,0,rB,nv=100,border=NULL,col="black",lty=1,density=NULL,
  #             angle=45,lwd=1)
  # draw.circle(0,0,rR,nv=100,border=NULL,lty=1,density=NULL,
  #             angle=45,lwd=1)
  # points(X[which(((X[,1])^2+(X[,2])^2)<=rB^2),],
  #        bg="red",col="red",pch=20)
  # dev.off()
}

# Save csv of agent distances from origin, i.e. associated r values:
dist_file = data.frame(distance_bolus = as.vector(t(sqrt((X.mat[,,1])^2+(X.mat[,,2])^2))),time = rep(0:(time_m-1),each=N))
write.csv(dist_file,file="distances.csv")
