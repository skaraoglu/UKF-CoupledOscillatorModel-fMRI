# Load necessary libraries
library(devtools)
install_github("insilico/UKF")  
library(UKF)
install.packages(c('KernSmooth'))
library(KernSmooth)
library(deSolve)

# parameters to optimize, p1 p2, in coupled_osc_model
param_guess <- c(10,1,4)

# first row is time, then voxel A and voxel B rows
# UKF_blend assumes column data, so transpose
vox.A.B.data <- read.delim("data/voxel_A&B_data.txt",sep="",header=F)
# make column-wise
vox.A.B.data <- t(vox.A.B.data)
# plot_and_smooth.R
smoothed_data <- plot_voxels_and_smooth(vox.A.B.data)

t_vec <- vox.A.B.data[,1] # time is first column
dT <- t_vec[2]-t_vec[1]   # assume uniform time steps dT=1
# smaller steps size for propagating model between dT steps
dt <- 0.1*dT
#nn <- round(dT/dt)
# num observed ind vars, first col is time, so -1
N_y <- ncol(vox.A.B.data)-1  # N_y=2
# number of unknown model parameters to be estimated
N_p <- 3
# size of augmented state vector
N_x <- N_p + N_y


# New Coupled Oscillator Model
new_coupled_osc_model <- function(t, x, p){
  # Parameters
  g <- p[1,]  # acceleration due to gravity
  L <- p[2,]  # length of pendulums
  k <- p[3,]  # spring constant

  # State variables
  theta1 <- x[1,]
  theta2 <- x[2,]

  # Equations of motion
  theta1_dot_dot <- -g/L*sin(theta1) - k/L*(sin(theta1) - sin(theta2))
  theta2_dot_dot <- -g/L*sin(theta2) + k/L*(sin(theta1) - sin(theta2))

  # Return rate of change
  return(rbind(theta1_dot_dot, theta2_dot_dot))
}
param_guess <- c(1,1,1)

ukf_out <- UKF_blend(t_dummy,vox.A.B.data,new_coupled_osc_model,
                     N_p,N_y,param_guess,dt,dT)
# ukf_out <- ukf(coupled_osc_model(t_vec, N_y, N_p))
ukf_out$param_est
ukf_out$chisq
# function in plot_and_smooth.R
plot_ukf_and_smoothed(ukf_out, smoothed_data,
                      top_title='One UKF Step, Raw Data')


iter_opt <- iterative_param_optim(param_guess,
                      t_dummy, smoothed_data,
                      #coupled_osc_model,
                      new_coupled_osc_model,
                      N_p,N_y,dt,dT,
                      param_tol=.01,MAXSTEPS=3000)
iter_opt$par   # params
iter_opt$value # chi-square
iter_opt$steps
iter_opt$param_norm
# plug iterative optim parameters back into UKF and plot
ukf_iter_opt <- UKF_blend(t_dummy,smoothed_data,
                           #coupled_osc_model,
                            new_coupled_osc_model,
                           N_p,N_y,iter_opt$par,dt,dT)
ukf_iter_opt$param_est
ukf_iter_opt$chisq
plot_ukf_and_smoothed(ukf_iter_opt, smoothed_data,
                      top_title = 'Iterative Optimization')
