# Simulating runoff errors with CoSMoS

setwd("/home/jschap/Documents/ISR/inverse_streamflow_routing/")

library(gstat)
library(CoSMoS)

# r <- raster("./basin_setup/basinmask.tif")

# Simulate a time series of error
n <- 120 # sample size/number of time steps
# -----------------------------

# Determining parameters
mean1=1
Cv <- 1
mu1 <- log((mean1^2)/sqrt(Cv^2+mean1^2));
sigma1 = sqrt(log(Cv^2/(mean1^2)+1));

# One grid cell is approximately 12.31 km
# To run in background:
# use nohup Rscript ./simulate*.R in command line

Lval <- 3.249391 # grid cells (40 km)
Tval <- 5 # days

M <-2000 # number of replicates

ns <- 21
# takes about 30 minutes per replicate on my computer
# requires about 12 GB of RAM
for (mm in 1:M)
{

  rf_sim <- generateRFFast(n=n, spacepoints = ns, p0 = 0,
                           margdist ="lnorm",
                           margarg = list(meanlog = mu1, sdlog = sigma1),
                           stcsarg = list(scfid = "weibull", tcfid = "weibull",
                                          scfarg = list(scale = Lval, shape = 1), # scale = L
                                          tcfarg = list(scale = Tval, shape = 1)) # scale = T
  )

  write.table(as.matrix(rf_sim),
              file = paste0("./allegheny_data/errors/m1a1L40T5_LM/alleg_err_sim_", mm, ".txt"),
              # file = paste0("/media/jschap/HD_ExFAT/Ohio/initial_runoff/cosmos_errors/Case4/sim_", mm, ".txt"),
              sep = ",",
              row.names = FALSE, col.names = FALSE
              )

}
