# Simulating runoff errors with CoSMoS

library(raster)
library(rgdal)
library(gstat)
library(CoSMoS)

setwd("/hdd/ISR/02_Basins/HH2")

# r <- raster("./Data/basinmask.tif")

# Simulate a time series of error
n <- 365 # sample size/number of time steps
marginal_dist <- "norm"
param <- list(mean = 0, sd = 1)
normal_sim <- generateTS(n, marginal_dist, param, acsvalue = c(1,0.8))
quickTSPlot(normal_sim[[1]])

# Simulate an spatial field of error
rf_sim <- generateRFFast(n = 2001, spacepoints = 7, p0 = 0, 
                         margdist ="lnorm",
                         margarg = list(meanlog = 0, sdlog = 1),
                         stcsarg = list(scfid = "weibull", tcfid = "weibull",
                         scfarg = list(scale = 1, shape = 1), # scale = L
                         tcfarg = list(scale = 24, shape = 1)) # scale = T
                         )

checkRF(rf_sim, nfields = 4, method = "field")

# Does anisotropy and/or advection occur? Do we need to include it?

# Is it possible to simulate errors whose marginal distribution changes in time/space?

# Simulate lognormal errors for the Hetch Hetchy basin -----------------------------

# Determining parameters
mean1 <- 2; # mean of perturbations
# mean1 <- 0.7536 # (using the total mean runoff for the time period a la PW13)
stddev <- 2
# stddev <- 100*0.7536; # standard deviation of perturbations
mu1 <- log((mean1^2)/sqrt(stddev^2+mean1^2));
sigma1 = sqrt(log(stddev^2/(mean1^2)+1));
v <- rlnorm(1e6, meanlog = mu1, sdlog  = sigma1)
mean(v)
sd(v)
hist(v)

# theoretical standard deviation
sqrt(exp(2*mu1 + sigma1^2)*(exp(sigma1^2)-1))

# Simulate an spatial field of error (M times)

# We may consider temporal correlations on the order of 0-168 hours
# We may consider spatial correlations on the order of 0-5 grid cells

# dimension of ohio: 72x104
# need spacepoints = 104

cell_length <- 6.17 # km
A <- 761.48 # basin area (km^2)
Lmin <- 0.01
Lmax <- sqrt(A)/cell_length
numTests <- 10
Linc <- seq(Lmin, Lmax, length.out = numTests)

# Tmin <- 0.001
# Tmax <- 24*7
# numTests <- 10
# Tinc <- seq(Tmin, Tmax, length.out = numTests)

for (mm in 1:500)
{
  
  rf_sim <- generateRFFast(n = 722, spacepoints = 7, p0 = 0, 
                           margdist ="lnorm",
                           margarg = list(meanlog = mu1, sdlog = sigma1),
                           stcsarg = list(scfid = "weibull", tcfid = "weibull",
                                          scfarg = list(scale = 10000, shape = 1), # scale = L
                                          tcfarg = list(scale = 0.01, shape = 1)) # scale = T
  )
  
  dirname <-paste0("/hdd/ISR/02_Basins/HH2/sensitivity_to_prior/varyingL/L", 10000)
  # dirname <-paste0("/hdd/ISR/02_Basins/HH2/sensitivity_to_prior/varyingL/L", Linc[testnum])
  if (mm==1)
  {
    dir.create(dirname)
  }
  write.table(as.matrix(rf_sim), 
              file = paste0(dirname, "/sim_", mm, ".txt"), 
              sep = ",", 
              row.names = FALSE, col.names = FALSE
  )
  
}  


# basin is 6x7, 1/16 degree grid cells
checkRF(rf_sim, nfields = 16, method = "field")
   


