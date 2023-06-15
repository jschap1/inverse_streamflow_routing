# Simulating runoff errors with CoSMoS

library(raster)
library(rgdal)
library(gstat)
library(CoSMoS)

setwd("/hdd/ISR/02_Basins/HH2")

# Simulate lognormal errors for the Hetch Hetchy basin -----------------------------

# cell_length <- 6.17 # km
# A <- 761.48 # basin area (km^2)
# Lmin <- 0.01
# Lmax <- sqrt(A)/cell_length
# numTests <- 10
# Linc <- seq(Lmin, Lmax, length.out = numTests)

# Parameter combos to test
L_vals <- c(0.01, 1.57, 7, 14, 21, 28, 56, 83, 538) # i
T_vals <- c(0.01, 1.5, 1.875, 3, 4.5, 6, 12, 18, 117) # j
COV <- 5
mu_vals <- c(0.5, 1, 2) # k
sigma_vals <- COV*mu_vals
# It takes about 8 hours to run this on my home desktop
# It produces about 75 GB of data 

M <- 500
ncombos <- length(L_vals)*length(T_vals)*length(mu_vals)

# cell_length <- 6.17 # km
# A <- 761.48 # basin area (km^2)
# sigma_min <- 0.1
# sigma_max <- 3
# # mu_min <- 0.01
# # mu_max <- 5
# numTests <- 10
# sigma_inc <- seq(sigma_min, sigma_max, length.out = numTests)
# mu_inc <- seq(mu_min, mu_max, length.out = numTests)

for (i in 1:length(L_vals))
{
  for (j in 1:length(T_vals))
  {
    for (k in 1:length(mu_vals))
    {
      
      # Determining parameters
      mean1 <- mu_vals[k]; # mean of perturbations
      stddev <- sigma_vals[k] # standard deviation of perturbations
      mu1 <- log((mean1^2)/sqrt(stddev^2+mean1^2));
      sigma1 = sqrt(log(stddev^2/(mean1^2)+1));
      
      for (mm in 1:M)
      {
        
        rf_sim <- generateRFFast(n = 722, spacepoints = 7, p0 = 0, 
                                 margdist ="lnorm",
                                 margarg = list(meanlog = mu1, sdlog = sigma1),
                                 stcsarg = list(scfid = "weibull", tcfid = "weibull",
                                                scfarg = list(scale = L_vals[i], shape = 1), # scale = L
                                                tcfarg = list(scale = T_vals[j], shape = 1)) # scale = T
        )
        
        fname <- paste0("L", L_vals[i], "_T", T_vals[j], "_mu", mu_vals[k], "_sigma", sigma_vals[k])
        dirname <-paste0("/hdd/ISR/02_Basins/HH2/Data/runoff_priors_cov5/", fname)
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
      
    }
  }
}

# basin is 6x7, 1/16 degree grid cells
# checkRF(rf_sim, nfields = 16, method = "field")



