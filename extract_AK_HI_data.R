# This code handles the special cases of the states of Alaska and Hawaii. These two have data for
# 25 models (only) at a different resolution.

library(ncdf4)
library(Rcpp)
library(dplyr)

netcdf.dir <- "/media//eitan/My Book/nex-gddp"
output.dir <- "/fast/ak_hi"

grid <- read.csv("/media/eitan/My Book/ak_hi_grid_cty.csv")

models <- c("ACCESS1-0", "BNU-ESM", "CCSM4", "CESM1", "CNRM-CM5", "CSIRO-Mk3-6-0", "CanESM2", "GFDL-CM3", "GFDL-ESM2G",
            "GFDL-ESM2M", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "MIROC-ESM-CHEM", "MIROC-ESM", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR",
            "MRI-CGCM3", "NorESM1-M", "bcc-csm1-1", "inmcm4")

read.file <- function(mtype, fn) {
  # Start by getting coordinates:
  origin <- paste0(substr(fn, nchar(fn) - 6, nchar(fn) - 3), "-01-01")
  print(paste("Starting type:", mtype, "from:", origin, "at time:", Sys.time()))
  
  nc <- nc_open(paste0(netcdf.dir, "/", fn))
  lon <- ncvar_get(nc, "lon");
  lat <- ncvar_get(nc, "lat");
  df <- data.frame(expand.grid(lon, lat))

  # Read in modeled variable and cleanup missing values:
  gvar <- ifelse(mtype == "surface_runoff", "surface runoff", ifelse(mtype == "total_runoff", "total runoff", mtype))
  var <- ncvar_get(nc, gvar);
  fillvalue <- ncatt_get(nc, "pr", "_FillValue")$value;
  var[var == fillvalue] <- NA;
  nc_close(nc);

  # Attach one column per day of modeled variable
  names(df) <- c("LON","LAT");
  for (i in 1:dim(var)[3]) {
    df <- cbind(df, as.numeric(var[,,i]))
    names(df)[i + 2] <- format(as.Date(i - 1, origin = origin), "%m%d%Y")
  }
  
  # Reduce to desired grid points:
  df$LON[df$LON > 180] <- df$LON[df$LON > 180] - 360
  df <- merge(grid, df)
  
  return(df)
}

# Convert all NetCDF files of a given type and model to a single RData file
extract.grid.data <- function(mtype, model) {
  setwd(netcdf.dir)
  df <- grid
  for (fn in list.files(pattern = paste0(mtype, ".*", model, "_[12][09][0-9][0-9].nc"))) {
    annual <- read.file(mtype, fn)
    df <- cbind(df, annual[3:ncol(annual)])
  }
  setwd(output.dir)
  save(df, file = paste(sep = ".", "all", mtype, model, "RData"))
}


# For a given subset of tmax, compute the number of times for each grid point where a 3-day
# run of temperatures exceeds a threshold for that grid point.
sum.heat.waves.str <- '
NumericVector fast_heat_waves(const NumericMatrix& tmax, const NumericVector& thresholds) {
  NumericVector sums(tmax.nrow());
  const int end = tmax.ncol();

  for (int i = 0; i < tmax.nrow(); ++i) {
    sums[i] = 0;
    for (int j = tmax.ncol() - 3; j >= 0; --j) {
      sums[i] += ((tmax(i, j) >= thresholds[i]) && (tmax(i, j + 1) >= thresholds[i]) && (tmax(i, j + 2) >= thresholds[i]))? 1 : 0;
    }
  }

  return(sums);
}
'
cppFunction(sum.heat.waves.str, plugins = c("cpp11"))

# For a given model, compute the mean annual no. of heatwave days for each grid point and
# for each of the three epochs.
compute.mean.heatwave.days <- function(model) {
  print(paste("Working on model:", model, "at time:", Sys.time()))
  load(paste0(output.dir, "/all.tasmax.", model, ".RData"))
  start.91 = which(names(df) == "01011991")
  end.10 = which(names(df) == "12312010")
  start.21 = which(names(df) == "01012021")
  end.40 = which(names(df) == "12312040")
  start.41 = which(names(df) == "01012041")
  end.60 = which(names(df) == "12312060")

  thresholds <- apply(df[,start.91:end.10], 1,
                      function(row) { quantile(as.numeric(row), 0.95, na.rm = TRUE)})
  
  return(data.frame(LON = df$LON, LAT = df$LAT,
                    mean.days.91 = fast_heat_waves(as.matrix(df[,start.91:end.10]), thresholds) / 20,
                    mean.days.21 = fast_heat_waves(as.matrix(df[,start.21:end.40]), thresholds) / 20,
                    mean.days.41 = fast_heat_waves(as.matrix(df[,start.41:end.60]), thresholds) / 20
  ))
}

# Helper code for aggregate.heatwave.days for repeated code: aggregate by county and state
aggregate.heatwave.aux <- function(metric) {
  county.means <- group_by(metric, GEOID) %>% summarise_each(funs(mean))
  county.pop <- countypop[as.character(county.means$GEOID),]
  county.sums <- cbind(as.integer(county.means$GEOID / 1000), county.pop * county.means[,2:ncol(county.means)])
  names(county.sums)[1] = "State"
  state.sums <- group_by(county.sums, State) %>% summarise_each(funs(sum))
  state.sums <- cbind(state.sums, t(apply(state.sums[2:ncol(state.sums)], 1, quantile)))
  return(state.sums)
}

# Compute mean annual heatwave days for each model, grid point, and epoch.
# Aggregate over counties and states and summarize across models for each epoch.
aggregate.heatwave.days <- function() {
  res.91 <- grid
  res.21 <- grid
  res.41 <- grid

  for (m in models) {
    days <- compute.mean.heatwave.days(m)
    res.91 <- cbind(res.91, days$mean.days.91)
    res.21 <- cbind(res.21, days$mean.days.21)
    res.41 <- cbind(res.41, days$mean.days.41)
  }

  names(res.91)[6:ncol(res.91)] = models
  names(res.21)[6:ncol(res.21)] = models
  names(res.41)[6:ncol(res.41)] = models

  # Get population per counties of interest:
  pop <- read.csv("/media/eitan/My Book/countypov.csv")
  countypop <- data.frame(pop = pop$vul.pop, row.names = pop$GEOID)

  write.csv(aggregate.heatwave.aux(res.91[,5:ncol(res.91)]),
            file = paste0(output.dir, "/aggregated-heatwave.1991.csv"), row.names = FALSE)
  write.csv(aggregate.heatwave.aux(res.21[,5:ncol(res.21)]),
            file = paste0(output.dir, "/aggregated-heatwave.2021.csv"), row.names = FALSE)
  write.csv(aggregate.heatwave.aux(res.41[,5:ncol(res.41)]),
            file = paste0(output.dir, "/aggregated-heatwave.2041.csv"), row.names = FALSE)
}
