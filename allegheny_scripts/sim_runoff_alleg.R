# Simulating runoff with CoSMoS
#
# Allegheny River Basin (7/27/2023)
# Goal is to simulate runoff time series for each cell in Allegheny, to serve as
# synthetic true runoff for correlation experiments.
#
# We will generate runoff fields with different L, T values and see how well ISR
# is able to estimate from synthetic truth discharge
#
# We can also check whether the optimal L, T values are those that were used to
# generate the synthetic true runoff fields in the first place

library(gstat)
library(CoSMoS)
# setwd("/home/jschap/Documents/ISR")
setwd("/Volumes/HD3/ISR/inverse_streamflow_routing")
# dir.create("./allegheny_data/runoff_sim")

n <- 365 # sample size/number of time steps
#-----------------------------

Cv <- 1; mean1<-1;
mu1 <- log((mean1^2)/sqrt(Cv^2+mean1^2));
sigma1 = sqrt(log(Cv^2/(mean1^2)+1));
# v <- rlnorm(1e6, meanlog = mu1, sdlog  = sigma1)
# hist(v, breaks = "FD")

# dimension of allegheny: 20x21
# One grid cell is approximately 12.31 km

# To run in background:
# use nohup Rscript ./simulate*.R in command line

Lstar <- sqrt(33750)/12.31 # grid cells
Tstar <- 4 # approx time of concentration

ns <- 21

#Lvals <- c(0.01, 3, 6, 9, 12, 16) # grid cells
#Tvals <- c(0.01, 1, 2, 3, 4) # days

Lvals <- seq(0.01, 50, length = 50)
Tvals <- seq(0.01, 25, length = 50)

M <- 1 # sample size

for (i in 1:length(Lvals))
{
  for (j in 1:length(Tvals))
  {
    for (mm in 1:M)
    {

      # normal
      rf_sim <- generateRFFast(n = 365, spacepoints = ns, p0 = 0,
                               margdist ="norm",
                               margarg = list(mean = 6*Cv, sd = Cv),
                               stcsarg = list(scfid = "weibull", tcfid = "weibull",
                                              scfarg = list(scale = Lvals[i], shape = 1), # scale = L
                                              tcfarg = list(scale = Tvals[j], shape = 1)) # scale = T
      # lognormal
      # rf_sim <- generateRFFast(n = 365, spacepoints = ns, p0 = 0,
      #                          margdist ="lnorm",
      #                          margarg = list(meanlog = mu1, sdlog = sigma1),
      #                          stcsarg = list(scfid = "weibull", tcfid = "weibull",
      #                                         scfarg = list(scale = Lvals[i], shape = 1), # scale = L
      #                                         tcfarg = list(scale = Tvals[j], shape = 1)) # scale = T
      )

      LL <- Lvals[i]
      TT <- Tvals[j]
      write.table(as.matrix(rf_sim),
                  file = paste0("./allegheny_data/runoff_sim_norm/L", sprintf("%.3f", LL), "_T", sprintf("%.3f", TT), "_", mm, ".txt"),
                  sep = ",",
                  row.names = FALSE, col.names = FALSE
      )

    }

  }
}

# Takes about 1 second per replicate

# Plot time series at grid cells:
# plot(rf_sim[,1], type="l")
