# This routine computes Finite-Time Lyapunov Exponents (FTLE) fields from a
# matrix of lagrangian particles positions evolving FORWARD in time.
#
# http://mmae.iit.edu/shadden/LCS-tutorial/computation.html
#
# MapsF 2012


# Load necessary libraries
library(ncdf4)
library(matrixStats)
library(RColorBrewer)
library(heatmaply)

# Start timer
tic()

# Load input NetCDF files
infile <- readline(prompt = 'Input file name root? ')
dt <- as.numeric(readline(prompt = 'Time step for computing FTLE (in days, even)? '))

ncid <- nc_open(paste0('/Users/frederic/Documents/JOB/Projets/Krill/Modelling/FTLE/', infile, '_FTLE.nc'), write = FALSE)

pdimid <- nc_inq_dimid(ncid, 'part')
tdimid <- nc_inq_dimid(ncid, 'time')

plen <- nc_inq_dim(ncid, pdimid)$size
tlen <- nc_inq_dim(ncid, tdimid)$size

lonid <- nc_inq_varid(ncid, 'lon')
latid <- nc_inq_varid(ncid, 'lat')
dxid <- nc_inq_varid(ncid, 'dx')
dyid <- nc_inq_varid(ncid, 'dy')
xid <- nc_inq_varid(ncid, 'xt')
yid <- nc_inq_varid(ncid, 'yt')
tid <- nc_inq_varid(ncid, 'time')

lon <- ncvar_get(ncid, lonid)
lat <- ncvar_get(ncid, latid)
dx <- ncvar_get(ncid, dxid)
dy <- ncvar_get(ncid, dyid)
time <- ncvar_get(ncid, tid)

m <- 197
n <- 234
mi <- 0
ni <- 0
mf <- m
nf <- n

start <- c(1, 1)
count <- c(m, n)

x <- ncvar_get(ncid, xid, start = start, count = count)
x[x > 1e9] <- NA

y <- ncvar_get(ncid, yid, start = start, count = count)
y[y > 1e9] <- NA

# Get the correct timesteps & convert coordinates in distances
X <- array(NA, dim = c(dim(x), 2))
Y <- array(NA, dim = c(dim(y), 2))

for (i in 1:dim(x)[3]) {
  X[, , i, 1] <- x[, , i, i]
  X[, , i, 2] <- x[, , i, i + round(dt * 0.5)]
  Y[, , i, 1] <- y[, , i, i]
  Y[, , i, 2] <- y[, , i, i + round(dt * 0.5)]
}

X[X < 1] <- 1
Y[Y < 1] <- 1
X[X > m] <- m
Y[Y > n] <- n

# Interpolate position in km according to the dx & dy matrices
DX <- matrix(0, nrow = nrow(dx) + 1, ncol = ncol(dx))
DX[2:nrow(DX), ] <- cumsum(dx, 1)

DY <- matrix(0, nrow = nrow(dy), ncol = ncol(dy) + 1)
DY[, 2:ncol(DY)] <- cumsum(dy, 2)

idx <- which(!is.na(X[, , 1, 1]), arr.ind = TRUE)

XX <- array(NA, dim = dim(X))
YY <- array(NA, dim = dim(Y))

for (k in 1:nrow(idx)) {
  i <- idx[k, 1]
  j <- idx[k, 2]
  for (l in 1:(dim(X)[3] * dim(X)[4])) {
    XX[i, j, l] <- DX[ceiling(X[i, j, l, 1]), ceiling(Y[i, j, l, 1])] +
      dx[ceiling(X[i, j, l, 1]), ceiling(Y[i, j, l, 1])] *
      (X[i, j, l, 1] - floor(X[i, j, l, 1]))
    YY[i, j, l] <- DY[ceiling(X[i, j, l, 1]), ceiling(Y[i, j, l, 1])] +
      dy[ceiling(X[i, j, l, 1]), ceiling(Y[i, j, l, 1])] *
      (Y[i, j, l, 1] - floor(Y[i, j, l, 1]))
  }
}

iX <- array(c(NA, XX[-nrow(XX), , , ]), dim = dim(XX))
Xi <- array(c(XX[-1, , , ], NA), dim = dim(XX))
jX <- array(c(NA, XX[, -ncol(XX), , ]), dim = dim(XX))
Xj <- array(c(XX[, -1, , ], NA), dim = dim(XX))

iY <- array(c(NA, YY[-nrow(YY), , , ]), dim = dim(YY))
Yi <- array(c(YY[-1, , , ], NA), dim = dim(YY))
jY <- array(c(NA, YY[, -ncol(YY), , ]), dim = dim(YY))
Yj <- array(c(YY[, -1, , ], NA), dim = dim(YY))

dXY <- array(NA, dim = c(dim(X), 2, 2))
dXY[, , , 1, 1] <- (Xi[, , , 2] - iX[, , , 2]) / (Xi[, , , 1] - iX[, , , 1])
dXY[, , , 1, 2] <- (Xj[, , , 2] - jX[, , , 2]) / (Yj[, , , 1] - jY[, , , 1])
dXY[, , , 2, 1] <- (Yi[, , , 2] - iY[, , , 2]) / (Xi[, , , 1] - iX[, , , 1])
dXY[, , , 2, 2] <- (Yj[, , , 2] - jY[, , , 2]) / (Yj[, , , 1] - jY[, , , 1])

idx <- which(!is.na(XX[, , 1, 1]), arr.ind = TRUE)

FTLE <- array(NA, dim = c(nrow(dx), ncol(dx), dim(XX)[3]))

for (k in 1:nrow(idx)) {
  i <- idx[k, 1]
  j <- idx[k, 2]
  for (t in 1:dim(XX)[3]) {
    if (!is.na(dXY[i, j, t, , ])) {
      FTLE[mi + i, ni + j, t] <- log(sqrt(max(eigen(t(dXY[i, j, t, , ]) %*% dXY[i, j, t, , ]))$values)) / dt
    }
  }
}

# Control plots
heatmaply(FTLE[, , 1], rev(unique(lat)), unique(lon), k_row = nrow(FTLE), k_col = ncol(FTLE),
          colors = colorRampPalette(brewer.pal(9, "RdBu"))(100),
          main = paste0("FTLE Field at Time ", time[1], " days"))

heatmaply(rowMeans(FTLE), rev(unique(lat)), unique(lon), k_row = nrow(FTLE), k_col = ncol(FTLE),
          colors = colorRampPalette(brewer.pal(9, "RdBu"))(100),
          main = "Mean FTLE Field")

# Save FTLE field
outfile <- paste0(infile, "_FTLE_", as.character(dt), "d")

saveRDS(list(infile = infile, dt = dt, lon = lon, lat = lat, time = time, FTLE = FTLE), file = outfile)

# Stop timer and print elapsed time
toc()
