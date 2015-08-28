# This code handles the special cases of the states of Alaska and Hawaii. These two have data for
# 25 models (only) at a different resolution.

library(ncdf4)
library(Rcpp)
library(dplyr)
library(doParallel);

netcdf.dir <- "/media//eitan/My Book/nex-gddp"
output.dir <- "/fast/ak_hi"

grid <- read.csv("/media/eitan/My Book/ak_hi_grid.csv")

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
  nc_close(nc)

  # Attach one column per day of modeled variable
  names(df) <- c("LON","LAT");
  for (i in 1:dim(var)[3]) {
    df <- cbind(df, as.numeric(var[,,i]))
    names(df)[i + 2] <- format(as.Date(i - 1, origin = origin), "%m%d%Y")
  }
  
  # Reduce to desired grid points:
  df$LON[df$LON > 180] <- df$LON[df$LON > 180] - 360
  df <- merge(grid[,1:2], df)

  return(df)
}

# Convert all NetCDF files of a given type and model to a single RData file
extract.grid.data <- function(mtype, model) {
  setwd(netcdf.dir)
  df <- grid
  for (fn in list.files(pattern = paste0(mtype, ".*", model, "_[12][09][0-9][0-9].nc"))) {
    annual <- read.file(mtype, fn)
  }
  setwd(output.dir)
  save(df, file = paste(sep = ".", "all", mtype, model, "RData"))
}

############################## Heatwave hazard ##########################################

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

# Helper code for repeated use: aggregate by county and state given a metric and county population.
# Assumes the first column in metric is the GEOID, and all other columns contain metrics to aggregate.
aggregate.by.county.and.state <- function(metric, countypop) {
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

  write.csv(aggregate.by.county.and.state(res.91[,5:ncol(res.91)], countypop),
            file = paste0(output.dir, "/aggregated-heatwave.1991.csv"), row.names = FALSE)
  write.csv(aggregate.by.county.and.state(res.21[,5:ncol(res.21)], countypop),
            file = paste0(output.dir, "/aggregated-heatwave.2021.csv"), row.names = FALSE)
  write.csv(aggregate.by.county.and.state(res.41[,5:ncol(res.41)], countypop),
            file = paste0(output.dir, "/aggregated-heatwave.2041.csv"), row.names = FALSE)
}


# Given a vector of daily precipitation amounts in mm, compute for each day
# the total precip in the 7 days ending in this day, in inches.
# For the first 6 days, no value is computed.
compute.weekly.precip <- function(daily.precip) {
  ret <- rep(0, 6)
  ret[7] <- sum(daily.precip[1:7])

  for (i in 8:length(daily.precip)) {
    ret[i] = ret[i - 1] + daily.precip[i] - daily.precip[i - 7]
  }

  ret = ret * 0.0393701  # Convert mm to inches
  return(ret)
}

############################# Wildfire hazard #################################################

# Given a vector of daily precipitation amounts in mm/s, compute for each day
# the total precip in the 7 days ending in this day, in inches.
# For the first 6 days, no value is computed.
weekly.precip.str <- '
NumericVector fast_weekly_precip(NumericVector daily_precip) {
  auto sz = daily_precip.size();
NumericVector ret(sz);

ret[0] = ret[1] = ret[2] = ret[3] = ret[4] = ret[5] = 0.;
ret[6] = daily_precip[0] + daily_precip[1] + daily_precip[2] + daily_precip[3] +
daily_precip[4] + daily_precip[5] + daily_precip[6];

for (auto i = 7; i < sz; ++i) {
ret[i] = ret[i - 1] + daily_precip[i] - daily_precip[i - 7];
}

return(ret * 0.0393701 * 60 * 60 * 24);  // Convert from mm/sec to in/day
}
'
cppFunction(weekly.precip.str, plugins = c("cpp11"))

# For a given series of max temperatures and precipitation, find the latest
# index prior to start.date that can be considered as KBDI=0, based on
# having the maximum weekly rainfall or 7.874in of rain, whicever is higher.
find.KBDI.start.index <- function(weekly.precip) {
  capped <- pmin(weekly.precip, 7.874)
  return(max(which(capped == max(capped))))
}

# Given an array of daily maximum temperatures (in C) and daily precipitation (in inches),
# compute the KBDI value for each day starting from the date (array index) for which KBDI
# is defined to be 0.
# The formula follows Liu et al., "Trends in global wildfire potential in a changing climate", in
# Forest Ecology and Management, 2010 (259), pp. 685--697. The algorithm is further elaborated
# in Janis et al., http://climate.geog.udel.edu/~climate/publication_html/Pdf/JJF_IJWF_02.pdf
# "Near-real time mapping of Keetch-Byram drought index in the south-eastern United States",
# Int. J. of Wildland Fire, 2002 (11), pp. 281--289
compute.KBDI.str <- '
NumericMatrix fast_kbdi(const NumericMatrix& daily_precip, const NumericMatrix& daily_tmax, const NumericVector& kbdi0_index) {
  const auto nr = daily_precip.nrow();
  const auto nc = daily_precip.ncol();
  NumericMatrix ret(nr, nc);
  std::fill(ret.begin(), ret.end(), 0.);

  for (auto j = 0; j < nr; ++j) {
    const auto first = kbdi0_index[j];  // Start computing KBDI here. Remember C++ is 0-index based!
    double Q = 0.;     // Latest KBDI value
    double cumP = 1000.;  // Cumulative precipitation up to now during wet spell

    // Compute mean annual precipitation:
    double R = 0.;
    for (auto i = first; i < nc; ++i)  R = R + daily_precip(j, i);
    R = 365 * R / (nc - first);      // Divide by no. of years.

    // Now, for each subsequent day, compute KBDI as a function of the previous days KBDI:
    for (auto i = first; i < nc; ++i) {
      const auto temp = daily_tmax(j, i);
      const auto p = daily_precip(j, i);

      // Change in saturation from yesterday, if temperature is high enough:
      const auto dQ = (temp <= 50.)?
                       0.
                     : 0.001 * (800 - Q) * (0.968 * exp(0.0486 * temp) - 8.3) /
                       (1. + 10.88 * exp(-0.0441 * R));

      // Change in precipitation (adjusted downward by 0.2in after a dry period of < 0.2in):
      cumP = (p > 0.)? cumP + p : 0.;
      // Have we already accumulated more than 0.2in before today in spell? if not, substract it
      const auto dP = 100. * ((cumP - p > 0.2)? p : std::max(cumP - 0.2, 0.));

      //      if (i <= first + 10)
      //        Rcpp::Rcout << "j:" << j << "   i:" << i << "   Q:" << Q << "   p:" << p << "   T:" << temp << "   dQ:" << dQ << "   dP:" << dP << std::endl;

      Q = std::max(Q + dQ - dP, 0.);
      ret(j, i) = Q;
    }
  }

  return (ret);
}
'
cppFunction(compute.KBDI.str, plugins = c("cpp11"))

# For a given model, compute its KBDI initialization date, then follow through from
# that date into the future, computing KBDI for each day.
compute.model.KBDI <- function(model) {
  print(paste("Reading inputs for model", model, "at time", Sys.time()))
  load(paste0(output.dir, "/all.tasmax.", model, ".RData"))
  tmax <- df
  load(paste0(output.dir, "/all.pr.", model, ".RData"))
  precip <- df

  print(paste("Computing KBDI==0 indices for model", model, "at time", Sys.time()))
  baseline.start <- which(names(precip) == "01011991") # The time from which we record measurements
  # The last date to consider for KBDI=0 start is the prior date. Compute weekly precip up to then:
  first <- which(names(precip) == "01011950")
  weekly.precip = apply(precip[,first:(baseline.start - 1)], 1, fast_weekly_precip)
  kbdi0 <- apply(weekly.precip, 2, find.KBDI.start.index)

  start.values <- data.frame(LON = precip$LON, LAT = precip$LAT, date = names(precip)[kbdi0 + first - 1], index = kbdi0,
                             weekly.precip.in = sapply(1:length(kbdi0), function(i) weekly.precip[kbdi0[i], i]))
  write.csv(start.values, paste0(output.dir, "/wildfire/start_at.", model, ".csv"), row.names = FALSE)

  print(paste("Computing KBDI values for model", model, "at time", Sys.time()))
  pr <- as.matrix(precip[,5:ncol(precip)]) * 0.0393701 * 86400  # Convert mm/sec to inches/day
  tm <- (as.matrix(tmax[,5:ncol(tmax)]) - 273.15) * 1.8 + 32.   # Convert K to F
  kbdi <- fast_kbdi(pr, tm, kbdi0)

  print(paste("Saving KBDI values for model", model, "at time", Sys.time()))
  kbdi <- cbind(grid, kbdi)
  names(kbdi) = names(precip)
  save(kbdi, file = paste0(output.dir, "/wildfire/kbdi.", model, ".RData"))

  print(paste("Done at", Sys.time()))
}

compute.wildfire.by.state <- function(threshold = 600, start.date = "01011991", end.date = "12312010") {
  annual.high.kbdi <- grid
  
  for (model in models) {
    print(paste("Reading data for model", model, "at time:", Sys.time()))
    load(paste0(output.dir, "/wildfire/kbdi.", model, ".RData"))
    start <- which(names(kbdi) == start.date)
    end <- which(names(kbdi) == end.date)
    high.days <- apply(kbdi[,start:end], 1, function(row) { sum(row >= threshold)} / 20)
    annual.high.kbdi <- cbind(annual.high.kbdi, high.days)
  }
  names(annual.high.kbdi) = c(names(grid), models)

  pop <- read.csv("/media/eitan/My Book/countypov.csv")
  countypop <- data.frame(pop = pop$POP2010, row.names = pop$GEOID)
  aggregated <- aggregate.by.county.and.state(annual.high.kbdi[,c(3,5:ncol(annual.high.kbdi))], countypop)

  outname <- paste0(output.dir, "wildfire/aggregated-kbdi-over", threshold, substr(start.date, 5, 8), substr(end.date, 5, 8))
  write.csv(aggregated, outname, row.names = FALSE)
}


############################################### Compute danger days  ################################

yr <- c(1:365);
lyr <- c(1:366);
yrday <- c(yr,lyr,yr,yr,yr,lyr,yr,yr,yr,lyr,yr,yr,yr,lyr,yr,yr,yr,lyr,yr,yr)

# Create all the INI input files to MTCLIM for AK/HI
create.ini.files <- function() {
  registerDoParallel(cores = 50);
  for (model in models) {
    for (startyr in c("1991", "2021", "2041")) {
      print(paste("Started writing ini files for model", model, "for start year", startyr, " at:", Sys.time()))
      foreach (i = 1:nrow(grid), .combine=c) %dopar% {
        dir <- paste0(output.dir, "/danger_days/mtclimata/", model, "/", startyr, "/")
        dir.create(dir,recursive = TRUE)
        fn <- paste0(startyr, "_", grid[i,"LAT"], "_" , grid[i, "LON"])
        fname <- paste0(dir, fn, ".ini");

        if (grid$elev[i] != -500) {
          cat(paste0("MTCLIM Initialization file for LAT ", grid[i, "LAT"], " LON ", grid[i, "LON"]), file = fname, sep = "\n");
          cat(paste0("Reference row in elev: ", i, " / model: ", model, " / start year: ", startyr), file = fname, sep = "\n", append = TRUE);
          cat("", file = fname, sep = "\n", append = TRUE)
          cat("IOFILES", file = fname, sep = "\n", append = TRUE)
          cat(paste0(fn, ".dat"), file = fname, sep = "\n", append = TRUE)
          cat(fn, file = fname, sep = "\n", append = TRUE)
          cat("", file = fname, sep = "\n", append = TRUE)
          cat("CONTROL", file = fname, sep = "\n", append = TRUE)
          cat("2", file = fname, sep = "\n", append = TRUE)
          cat("7305", file = fname, sep = "\n", append = TRUE)
          cat("0", file = fname, sep = "\n", append = TRUE)
          cat("1", file = fname, sep = "\n", append = TRUE)
          cat("1", file = fname, sep = "\n", append = TRUE)
          cat("", file = fname, sep = "\n", append = TRUE)
          cat("PARAMETERS", file = fname, sep = "\n", append = TRUE)
          cat(grid[i, "elev"], file = fname, sep = "\n", append = TRUE)
          cat("1.0", file = fname, sep = "\n", append = TRUE)
          cat(grid[i, "LAT"], file = fname, sep = "\n", append = TRUE)
          cat(grid[i, "elev"], file = fname, sep = "\n", append = TRUE)
          cat("0", file = fname, sep = "\n", append = TRUE)
          cat("0", file = fname, sep = "\n", append = TRUE)
          cat("1.0", file = fname, sep = "\n", append = TRUE)
          cat("0", file = fname, sep = "\n", append = TRUE)
          cat("0", file = fname, sep = "\n", append = TRUE)
          cat("0", file = fname, sep = "\n", append = TRUE)
          cat("0", file = fname, sep = "\n", append = TRUE)
          cat("", file = fname, sep = "\n", append = TRUE)
          cat("END", file = fname, sep = "\n", append = TRUE)
        }
      }
      print(paste("Finished writing ini files for model",model,"for start year",startyr," at: ",Sys.time()))
    }
  }
}


# Create data files for MTCLIM
create.dat.files <- function() {
  for (model in models) {
    print(paste("Loading tasmax for", model, "at", Sys.time()))
    load(paste0(output.dir, "/all.tasmax.", model, ".RData"))
    tmax <- cbind(grid, df[,5:ncol(df)] - 273.15)  # Convert K to C
    print(paste("Loading tasmin for", model, "at", Sys.time()))
    load(paste0(output.dir, "/all.tasmin.", model, ".RData"))
    tmin <- cbind(grid, df[,5:ncol(df)] - 273.15)  # Convert K to C
    load(paste0(output.dir, "/all.pr.", model, ".RData"))
    precip <- cbind(grid, df[,5:ncol(df)] * 86400 / 10)  # Convert mm/sec to cm/day

    for (startyr in c("1991", "2021", "2041")) {
      start <- which(names(tmax) == paste0("0101", startyr))
      end <- which(names(tmax) == paste0("1231", as.integer(startyr) + 19))
      dir <- paste0(output.dir, "/danger_days/mtclimata/", model, "/", startyr, "/")
      dir.create(dir,recursive = TRUE)

      foreach (i = 1:nrow(grid), .combine=c) %dopar% {
        if (grid[i, "elev"] != -500) {
          fn <- paste0(startyr, "_", grid[i, "LAT"], "_" , grid[i, "LON"])
          fname <- paste0(dir, fn, ".dat");
          cat("Year", "YearDay\t", "Tmax", "Tmin", "Prcp", "\n", file = fname, sep="\t\t");
          cat("\t\t\t", "(deg C)", "(deg C)", "(cm)", "\n", file= fname, sep="\t\t", append = TRUE)
          mintemp <- as.numeric(tmin[i, start:end]);
          maxtemp <- as.numeric(tmax[i, start:end]);
          pr <- as.numeric(precip[i, start:end]);
          labels <- colnames(tmax)[start:end];
          year <- as.numeric(substr(labels, 5, 8));
          df <- as.data.frame(year);
          df <- cbind(df, yrday, maxtemp, mintemp, pr);
          write.table(df, file = fname, append = TRUE, row.names = FALSE, col.name = FALSE);
        }
      }
    }
  }
}

# At this point, you need to run mtclim on all the .ini files. Then run go() from combine_humidity.R
# for each appropriate directory. Then, the code from heat_hazard.R can be used to compute heat index.
# Finally, you can call aggregate.by.county.and.state from this file to combine all the results.
