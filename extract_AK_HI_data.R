# This code handles the special cases of the states of Alaska and Hawaii. These two have data for
# 25 models (only) at a different resolution.

library(ncdf4)

netcdf.dir <- "/media//eitan/My Book/nex-gddp"
output.dir <- "/fast/ak_hi"

grid <- read.csv("/media/eitan/My Book/grid_ak_hi.csv")

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