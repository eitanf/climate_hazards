# Compute wildfire risk based on a KBDI estimate (http://www.srs.fs.usda.gov/pubs/rp/uncaptured/rp_se038.pdf)

library(dplyr)
library(Rcpp)

models.dir <- "/data/modeldata"
wildfire.dir <- "/data/wildfire"
grid <- read.csv("/fast/grid.classify.csv")[,c(1,2,3,7)]

models <- c("access1-0_rcp85_r1i1p1", "bcc-csm1-1-m_rcp85_r1i1p1", "bcc-csm1-1_rcp85_r1i1p1", "canesm2_rcp85_r1i1p1",
            "ccsm4_rcp85_r1i1p1", "cesm1-bgc_rcp85_r1i1p1", "cesm1-cam5_rcp85_r1i1p1", "cmcc-cm_rcp85_r1i1p1",
            "cnrm-cm5_rcp85_r1i1p1", "csiro-mk3-6-0_rcp85_r1i1p1", "fgoals-g2_rcp85_r1i1p1", "fio-esm_rcp85_r1i1p1",
            "gfdl-cm3_rcp85_r1i1p1", "gfdl-esm2g_rcp85_r1i1p1", "gfdl-esm2m_rcp85_r1i1p1", "giss-e2-r_rcp85_r1i1p1",
            "hadgem2-ao_rcp85_r1i1p1", "hadgem2-cc_rcp85_r1i1p1", "hadgem2-es_rcp85_r1i1p1", "inmcm4_rcp85_r1i1p1",
            "ipsl-cm5a-mr_rcp85_r1i1p1", "ipsl-cm5b-lr_rcp85_r1i1p1", "miroc-esm-chem_rcp85_r1i1p1", "miroc-esm_rcp85_r1i1p1",
            "miroc5_rcp85_r1i1p1", "mpi-esm-lr_rcp85_r1i1p1", "mpi-esm-mr_rcp85_r1i1p1", "mri-cgcm3_rcp85_r1i1p1", "noresm1-m_rcp85_r1i1p1")

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

# Same as above, but implemented in C++11 for speed
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

  return(ret * 0.0393701);
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

load.data <- function(grid, model, mtype) {
  tmax <- grid
  for (year in c(1951, 1971, 1991, 2011, 2021, 2041)) {
    print(year)
    load(paste0(models.dir, "/", model, "/conus_c5.", model, ".daily.", year, ".", mtype, ".Rdata"))
    df <- df[!is.na(df[,3]),]
    tmax <- merge(tmax, df)
  }
  return(tmax)
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
  NumericMatrix ret(daily_precip.nrow(), daily_precip.ncol());

  for (auto j = 0; j < nr; ++j) {
    const auto first = kbdi0_index[j];  // Start computing KBDI here. Remember C++ is 0-index based!
    double Q = 0.;     // Latest KBDI value
    double cumP = 1000.;  // Cummulative precipitation up to now during wet spell

    // Compute mean annual preciptation:
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

      Q = std::max(Q + dQ - dP, 0.);
      ret(j, i) = Q;
    }
  }

  return (ret);
}
'

cppFunction(compute.KBDI.str, plugins = c("cpp11"))

# An R implementation of fast_kbdi, for debugging purposes
compute.kbdi <- function(daily.precip, daily.tmax, kbdi0.index) {
  nr = nrow(daily.precip)
  nc = ncol(daily.precip)
  ret = matrix(rep(0, nr * nc), nrow = nr, ncol = nc)

  for (j in 1:nr) {
    first = kbdi0.index[j] + 1  # Start computing KBDI here.
    Q = 0.   # Latest KBDI value
    cumP = 1000.  # Cummulative precipitation up to now during wet spell

    # Compute mean annual preciptation:
    R = 0.
    for (i in first:nc) R = R + daily.precip[j, i]
    R = 365 * R / (nc - first)      # Divide by no. of years.

    # Now, for each subsequent day, compute KBDI as a function of the previous days KBDI:
    for (i in first:nc) {
      temp = daily.tmax[j, i];
      #Change in saturation from yesterday, if temperature is high enough:
      dQ = ifelse(temp <= 50, 0,
        0.001 * (800 - Q) * (0.968 * exp(0.0486 * temp) - 8.3) /
        (1. + 10.88 * exp(-0.0441 * R)))

      # Change in precipitation (adjusted downward by 0.2in after a dry period of < 0.2in):
      p = daily.precip[j, i]
      cumP = ifelse(p > 0., cumP + p, 0.)
      # Have we already accumulated more than 0.2in before today in spell? if not, substract it
      dP = 100. * ifelse(cumP - p > 0.2, p, max(cumP - 0.2, 0.))

      Q = max(Q + dQ - dP, 0.)
      ret[j, i] = Q
    }
  }

  return (ret);
}

compute.model.KBDI <- function(model) {
  print(paste("Reading inputs for model", model, "at time", Sys.time()))
  tmax <- load.data(grid, model, "tasmax")
  precip <- load.data(grid, model, "pr")

  print(paste("Computing KBDI==0 indices for model", model, "at time", Sys.time()))
  baseline.start <- which(names(precip) == "01011991") # The time from which we record measurements
  # The last date to consider for KBDI=0 start is the prior date. Compute weekly precip up to then:
  weekly.precip = apply(precip[,5:(baseline.start - 1)], 1, fast_weekly_precip)
  kbdi0 <- apply(weekly.precip, 2, find.KBDI.start.index)

  print(paste("Computing KBDI values for model", model, "at time", Sys.time()))
  pr <- as.matrix(precip[,5:ncol(precip)]) * 0.0393701  # Convert to inches
  tm <- as.matrix(tmax[,5:ncol(tmax)]) * 9./5. + 32.    # Convert to F
  kbdi <- fast_kbdi(pr, tm, kbdi0)

  print(paste("Saving KBDI values for model", model, "at time", Sys.time()))
  write.csv(cbind(precip[,1:4], kbdi),
            bzfile(paste0(kbdi.dir, "/kbdi.", model, ".csv.bz2")),
            row.names = FALSE)

  print(paste("Done at", Sys.time()))
}