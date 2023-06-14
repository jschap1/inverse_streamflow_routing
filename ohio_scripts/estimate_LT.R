# Estimating L and T for NLDAS vs. TMPA runoff
#
# 6/12/2023 JRS

library(sp)
library(gstat)
library(rgdal)

# Load in NLDAS - TMPA errors ----------------------------------------

setwd("/home/jschap/Documents/ISR/inverse_streamflow_routing/")
err <- as.matrix(read.table("./ohio_data/nldas_tmpa_error.txt", 
                  header = FALSE
                  ))

# Determine the dimensions of the 2D matrix
nr <- 72
nc <- 104
nt <- 365

# Reconstruct the 3D array from the 2D matrix
err.3D <- array(err, dim = c(nr, nc, nt))

# Do time series analysis -----------------------------------------

gi <- which(!is.na(err.3D[,,1]), arr.ind = TRUE) # index of cells in basin
n <- length(gi[,1]) # number of grid cells in basin

errmat <- matrix(0, nrow = n, ncol = nt)
for (i in 1:n)
{
  errmat[i,] <- err.3D[gi[i,1], gi[i,2],]  
}
errdf <- as.data.frame(t(errmat))

tv <- seq(from = as.Date("2009-01-01"), by = "day", length.out = nt)
errdf <- cbind(tv, errdf)
colnames(errdf) <- c("Time", paste0("Cell_", 1:3681))

maxlag <- vector(length = n)
for (c in 1:n)
{
  val <- acf(errdf[,c+1], plot = FALSE)
  consecutive_numbers <- which(val$acf>2/sqrt(nt))
  # For white noise, 95% of ACF values should be in ±2/√T, where T = 365
  maxlag[c] <- which(diff(consecutive_numbers)>1)[1]
  if (is.na(maxlag[c]))
  {
    maxlag[c] <- consecutive_numbers[length(consecutive_numbers)]
  }
}

# Plot ACF for a particular cell
val <- acf(errdf[,4], plot = TRUE)

# ACF for basin average (and lag)
err.basin.avg <- colMeans(errmat)
acf.basin.avg <- acf(as.data.frame(err.basin.avg), plot=TRUE)

# Fit a correlation function of the form rho(tau) = exp(-tau/T)
acf.fit.data <- data.frame(lag = acf.basin.avg$lag, 
                           acf = acf.basin.avg$acf,
                           logacf = log(acf.basin.avg$acf)
                           )

fit1 <- lm(log(acf)~lag+0, acf.fit.data)
summary(fit1)

# Obtain the predicted values from the linear regression model
predicted <- exp(predict(fit1))

# Create a scatter plot of the original data
plot(acf.basin.avg$lag, acf.basin.avg$acf, 
     xlab = "lag", 
     ylab = "Temporal correlation", 
     main = "Fit vs Data"
     )

# Add the fitted values to the plot
lines(acf.basin.avg$lag, predicted, col = "red")
T.basin.avg <- -1/fit1$coefficients[1]

# Do this for each grid cell and report distribution of T values

T <- vector(length = n)
for (c in 1:n)
{
  T[c] <- calc_T(errmat[c,])  
}

mean(T,na.rm=TRUE)
summary(T)
hist(T, "fd")

# Functions -------------------------------------------------------

calc_T <- function(ts)
{
  
  nt <- 365

  acf.cell <- acf(as.data.frame(ts), plot=FALSE)
  
  # # Keep only significant lag correlations (if desired)
  # consecutive_numbers <- which(acf.cell$acf>2/sqrt(nt))
  # # For white noise, 95% of ACF values should be in ±2/√T, where T = 365
  # maxlag <- which(diff(consecutive_numbers)>1)[1]
  # if (is.na(maxlag))
  # {
  #   maxlag <- consecutive_numbers[length(consecutive_numbers)]
  # }
  # acf.cell <- acf.cell[0:maxlag,]
  
  # Fit a correlation function of the form rho(tau) = exp(-tau/T)
  acf.cell.data <- data.frame(lag = acf.cell$lag, 
                             acf = acf.cell$acf
  )
  
  fit1 <- lm(log(acf)~lag+0, acf.cell.data)
  # summary(fit1)
  
  # Obtain the predicted values from the linear regression model
  predicted <- exp(predict(fit1))
  
  # # Create a scatter plot of the original data
  # plot(acf.cell$lag, acf.cell$acf, 
  #      xlab = "lag", 
  #      ylab = "Temporal correlation", 
  #      main = "Fit vs Data"
  # )
  # 
  # # Add the fitted values to the plot
  # lines(acf.cell$lag, predicted, col = "red")
  T <- -1/fit1$coefficients[1]
  
  return(T)
  
}

# Spatial correlation analysis ----------------------------------------

distmat <- read.table("./ohio_data/distmat.txt", 
                      header = FALSE)

L <- vector(length = nt)
for (tt in 1:nt)
{
  runoff_error_snapshot <- paste0("./ohio_data/runoff_error_xyz/xyz", tt, ".txt")
  L[tt] <- calc_L(runoff_error_snapshot)
}

calc_L <- function(runoff_error_snapshot)
{
  
  # runoff_error_snapshot <- "./ohio_data/runoff_error_xyz/xyz1.txt"
  
  xyz <- read.table(runoff_error_snapshot,
                    header = FALSE, col.names = c("lon","lat","err"))
  
  spdf <- data.frame(lon = xyz$lon, lat = xyz$lat, err = xyz$err)
  coordinates(spdf) <- c("lon", "lat")
  crs <- CRS("+proj=longlat +datum=WGS84")
  
  # Assign the CRS to the spdf object
  proj4string(spdf) <- crs
  
  variogram_model <- variogram(err ~ 1, data = spdf)
  # plot(variogram_model, xlab = "distance (km)")
  
  vfit <- fit.variogram(variogram_model, 
                        model = vgm("Exp", nugget = 0))
  L <- vfit$range[2] # km
  
  # plot(vfit)
  
  return(L)
  
}


