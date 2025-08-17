# Simulating runoff errors with CoSMoS

library(gstat)
library(CoSMoS)

# r <- raster("./basin_setup/basinmask.tif")

# Simulate a time series of error
n <- 365 # sample size/number of time steps
# -----------------------------

# Determining parameters
mean1=1
Cv <- 10
mu1 <- log((mean1^2)/sqrt(Cv^2+mean1^2));
sigma1 = sqrt(log(Cv^2/(mean1^2)+1));
# v <- rlnorm(1e6, meanlog = mu1, sdlog  = sigma1)
# hist(v, breaks = "FD")

# dimension of ohio: 72x104
# need spacepoints = 104
# One grid cell is approximately 12.31 km

# To run in background:
# use nohup Rscript ./simulate*.R in command line

Lval <- 10 # grid cells
Tval <- 2 # days

ns <- 21
# takes about 30 minutes per replicate on my computer
# requires about 12 GB of RAM
for (mm in 1:200)
{

  rf_sim <- generateRFFast(n = 365, spacepoints = ns, p0 = 0,
                           margdist ="lnorm",
                           margarg = list(meanlog = mu1, sdlog = sigma1),
                           stcsarg = list(scfid = "weibull", tcfid = "weibull",
                                          scfarg = list(scale = Lval, shape = 1), # scale = L
                                          tcfarg = list(scale = Tval, shape = 1)) # scale = T
  )

  write.table(as.matrix(rf_sim),
              file = paste0("./ohio_err_sim_", mm, ".txt"),
              # file = paste0("/media/jschap/HD_ExFAT/Ohio/initial_runoff/cosmos_errors/Case4/sim_", mm, ".txt"),
              sep = ",",
              row.names = FALSE, col.names = FALSE
              )

}
