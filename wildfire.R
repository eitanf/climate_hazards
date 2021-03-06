# Compute wildfire risk based on a KBDI estimate (http://www.srs.fs.usda.gov/pubs/rp/uncaptured/rp_se038.pdf)

library(plyr)
library(dplyr)
library(Rcpp)
library(ggplot2)
library(sp)
library(maptools)
library(maps)

models.dir <- "/data/modeldata"
wildfire.dir <- "/data/wildfire"
grid <- read.csv("/fast/grid.classify.csv")[,c(1,2,3,7)]

map.counties <- map_data("county")
map.counties$county <- with(map.counties, paste(region, subregion, sep = ","))
map.states <- map_data("state")
counties <- map('county', fill = TRUE, col = "transparent", plot = FALSE)
IDs <- sapply(strsplit(counties$names, ":"), function(x) x[1])
counties_sp <- map2SpatialPolygons(counties, IDs = IDs,
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
pointsSP <- SpatialPoints(as.data.frame(grid[,c('LON','LAT')]),
                          proj4string = CRS("+proj=longlat +datum=WGS84"))
indices <- over(pointsSP, counties_sp)
county.names <- sapply(counties_sp@polygons, function(x) x@ID)

# US counties that have no grid points:
missing.counties <- c(6075, 25019, 29510, 44001, 53055, 51720, 51520, 51640, 51750, 51775, 51580, 51530, 51790, 51660,
                      51540, 51600, 51013, 51683, 51685, 51630, 51670, 51570, 51830, 51735, 51595, 51610, 51678, 51690, 51710)
# The nearest respective counties to use for the missing.counties
replacement.counties <- c(6081, 25007, 29189, 44005, 53029, 51195, 51191, 51077, 51121, 51161, 51005, 51163, 51015,
                          51165, 51003, 51059, 51059, 51153, 51153, 51177, 51149, 51041, 51095, 51650, 51081, 51059,
                          51163, 51089, 51740)

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

# An R implementation of fast_kbdi, for debugging purposes
compute.kbdi <- function(daily.precip, daily.tmax, kbdi0.index) {
  nr = nrow(daily.precip)
  nc = ncol(daily.precip)
  ret = matrix(rep(0, nr * nc), nrow = nr, ncol = nc)

  for (j in 1:nr) {
    first = kbdi0.index[j] + 1  # Start computing KBDI here.
    Q = 0.   # Latest KBDI value
    cumP = 1000.  # Cumulative precipitation up to now during wet spell

    # Compute mean annual precipitation:
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

#      if (i <= first + 10)
#        print(paste("j:", j, "i:", i, "Q:", Q, "p:", p, "T:", temp, "dQ:", dQ, "dP:", dP))

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

  start.values <- data.frame(LON = precip$LON, LAT = precip$LAT, date = names(precip)[kbdi0 + 4], index = kbdi0,
                             weekly.precip.in = sapply(1:length(kbdi0), function(i) weekly.precip[kbdi0[i], i]))
  write.csv(start.values, paste0(wildfire.dir, "/start_at.", model, ".csv"), row.names = FALSE)

  print(paste("Computing KBDI values for model", model, "at time", Sys.time()))
  pr <- as.matrix(precip[,5:ncol(precip)]) * 0.0393701  # Convert to inches
  tm <- as.matrix(tmax[,5:ncol(tmax)]) * 9./5. + 32.    # Convert to F
  kbdi <- fast_kbdi(pr, tm, kbdi0)

  print(paste("Saving KBDI values for model", model, "at time", Sys.time()))
  kbdi <- cbind(precip[,1:4], kbdi)
  names(kbdi)[5:ncol(kbdi)] = names(precip)[5:ncol(kbdi)]
  save(kbdi, file = paste0(wildfire.dir, "/kbdi.", model, ".Rdata"))

  print(paste("Done at", Sys.time()))
}

# For a given model and grid point, this function avearages the KBDI values over three
# period of 20 years, during June, July, and August only. Returns a dataframe with the
# three means per grid point, for all grid points.
aggregate.years <- function(model) {
  print(paste("Reading", model, "at", Sys.time()))
  load(file = paste0(wildfire.dir, "/kbdi.", model, ".Rdata"))

  summers.91 = (substr(names(kbdi), 1, 2) >= "06"   & substr(names(kbdi), 1, 2) <= "08" &
                substr(names(kbdi), 5, 8) >= "1991" & substr(names(kbdi), 5, 8) <= "2010")
  summers.21 = (substr(names(kbdi), 1, 2) >= "06"   & substr(names(kbdi), 1, 2) <= "08" &
                substr(names(kbdi), 5, 8) >= "2021" & substr(names(kbdi), 5, 8) <= "2040")
  summers.41 = (substr(names(kbdi), 1, 2) >= "06"   & substr(names(kbdi), 1, 2) <= "08" &
                substr(names(kbdi), 5, 8) >= "2041" & substr(names(kbdi), 5, 8) <= "2060")

  means <- data.frame(mean.summers.91 = apply(kbdi[summers.91], 1, mean),
                      mean.summers.21 = apply(kbdi[summers.21], 1, mean),
                      mean.summers.41 = apply(kbdi[summers.41], 1, mean))

  ret <- cbind(kbdi[,1:4], means)
  write.csv(ret, paste0(wildfire.dir, "/summer-means.", model, ".csv"))
  return(ret)
}


start.date <- function(point.data) {
  ret <- c(LON = point.data[[1]]$LON,
           LAT = point.data[[1]]$LAT,
           median.idx = median(sapply(point.data, function(r) r$index)),
           min.year = as.integer(min(sapply(point.data, function(r) r$index)) / 365 + 1950))

  return(ret)
}

# For a geo given data (already in color buckets), plot it on a US map.
plot.on.map <- function(data, breaks, title, ylabel) {
  plot <- ggplot(data, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = colorBuckets)) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_fill_brewer(palette="PuRd", labels = breaks, name = ylabel) +
    coord_map(project="globular") +
    ggtitle(title) +
    geom_path(data = map.states, colour = "black", size = .3) +
    geom_path(data = map.counties, colour = "white", size = .5, alpha = .1)

  return(plot)
}

# This function computes for each grid point the median date for KBDI=0 start across all
# the models, the uses this data to plot it as a heatmap overlayed on a US county map.
plot.start.dates <- function() {
  files <- lapply(models, function(m) paste0(wildfire.dir, "/start_at.", m, ".csv"))
  data <- lapply(files, read.csv)

  dates <- data.frame(LON = numeric(), LAT = numeric(), median.idx = integer(), min.year = integer())

  for (i in 1:nrow(data[[1]])) {
    dates <- rbind(dates, start.date(lapply(data, "[", i,)))
  }
  names(dates) = c("LON", "LAT", "median.start.day", "min.year")
  write.csv(dates, paste0(wildfire.dir, "/start_dates_across_models.csv"), row.names = FALSE)

  dates <- read.csv(paste0(wildfire.dir, "/start_dates_across_models.csv"))

  dates$county <- county.names[indices]
  dates.mean <- ddply(dates, c("county"), summarize, means = mean(median.start.day))
  dates.min <- ddply(dates, c("county"), summarize, mins = min(min.year))
  dates.un <- dates[!duplicated(dates[c("county")]),]
  dates.1 <- merge(dates.mean, dates.un, by = 'county')
  dates.2 <- merge(dates.min, dates.un, by = 'county')

  breaks <- c(1000, 5000, 7500, 10000, 12500, 15000)
  dates.1$colorBuckets <- as.factor(as.numeric(cut(dates.mean$means, breaks)))
  mergedata <- merge(map.counties, dates.1, by = "county")
  mergedata <- mergedata[order(mergedata$means),]

  pdf(paste0(wildfire.dir, "/start_dates.pdf"))
  print(plot.on.map(mergedata, breaks, "Median start day, averaged per county", "KBDI=0 start day since 1950/01/01"))
  dev.off()
}

# Create a plot (with 29 charts, one per model) mapping the quantity of precipitation
# that was used to start the KBDI computation (KBDI=0) for each grid point.
plot.start.precip <- function() {
  files <- lapply(models, function(m) paste0(wildfire.dir, "/start_at.", m, ".csv"))
  data <- lapply(files, read.csv)

  pdf(paste0(wildfire.dir, "/start_precip.pdf"), onefile = TRUE)
  for (i in 1:29) {
    cur <- data[[i]]
    cur$county <- county.names[indices]
    rain.mean <- ddply(cur, c("county"), summarize, means = mean(weekly.precip.in))
    rain.un <- cur[!duplicated(cur[c("county")]),]
    rain <- merge(rain.mean, rain.un, by = 'county')
  
    breaks <- c(2, 3, 4, 5, 6, 7, 8, 10, 12, 14)
    rain$colorBuckets <- as.factor(as.numeric(cut(rain.mean$means, breaks)))
    mergedata <- merge(map.counties, rain, by = "county")
    mergedata <- mergedata[order(mergedata$means),]

    print(plot.on.map(mergedata, breaks, paste("Weekly precip averaged per county for model", models[i]), "Inch"))
  }
  dev.off()
}

# Plot the total no. of days above a given KBDI threshold for a given period.
# Plots three maps, one for each of 25%, 50%, 75% percentile results across models.
plot.high.KBDI.days <- function(threshold = 600, start.date = "01011991", end.date = "12312010") {
  annual.high.kbdi <- grid

  for (model in models) {
    print(paste("Reading data for model", model))
    load(paste0(wildfire.dir, "/kbdi.", model, ".Rdata"))
    start <- which(names(kbdi) == start.date)
    end <- which(names(kbdi) == end.date)
    high.days <- apply(kbdi[,start:end], 1, function(row) { sum(row >= threshold)})
    annual.high.kbdi <- cbind(annual.high.kbdi, high.days)
    rm(kbdi)
    gc()
  }

  basename <- paste(sep = "-", "/kbdi-over", threshold, substr(start.date, 5, 8), substr(end.date, 5, 8))
  names(annual.high.kbdi) = c(names(grid), models)
  write.csv(annual.high.kbdi, paste0(wildfire.dir, basename, ".csv"), row.names = FALSE)

  annual.high.kbdi <- read.csv(paste0(wildfire.dir, basename, ".csv"))
  tper <- apply(annual.high.kbdi[,5:ncol(annual.high.kbdi)], 1, quantile)
  percentiles <- cbind(grid[,1:2], t(tper)[,2:4])
  names(percentiles)[3:5] = c("p25", "p50", "p75")

  pdf(paste0(wildfire.dir, basename, ".pdf"), onefile = TRUE)
  breaks <- c(100, 200, 300, 400, 500, 750, 1000, 2000, 4000, 5000)
  percentiles$county <- county.names[indices]

  per.mean <- ddply(percentiles, c("county"), summarize, means = mean(p25))
  per.un <- percentiles[!duplicated(percentiles[c("county")]),]
  per <- merge(per.mean, per.un, by = 'county')
  per$colorBuckets <- as.factor(as.numeric(cut(per.mean$means, breaks)))
  mergedata <- merge(map.counties, per, by = "county")
  mergedata <- mergedata[order(mergedata$means),]
  print(plot.on.map(mergedata, breaks, paste0("KBDI>=", threshold, ", p25 across models"), "Days [1991-2010]"))

  per.mean <- ddply(percentiles, c("county"), summarize, means = mean(p50))
  per.un <- percentiles[!duplicated(percentiles[c("county")]),]
  per <- merge(per.mean, per.un, by = 'county')
  per$colorBuckets <- as.factor(as.numeric(cut(per.mean$means, breaks)))
  mergedata <- merge(map.counties, per, by = "county")
  mergedata <- mergedata[order(mergedata$means),]
  print(plot.on.map(mergedata, breaks, paste0("KBDI>=", threshold, ", p50 across models"), "Days [1991-2010]"))

  per.mean <- ddply(percentiles, c("county"), summarize, means = mean(p75))
  per.un <- percentiles[!duplicated(percentiles[c("county")]),]
  per <- merge(per.mean, per.un, by = 'county')
  per$colorBuckets <- as.factor(as.numeric(cut(per.mean$means, breaks)))
  mergedata <- merge(map.counties, per, by = "county")
  mergedata <- mergedata[order(mergedata$means),]
  print(plot.on.map(mergedata, breaks, paste0("KBDI>=", threshold, ", p75 across models"), "Days [1991-2010]"))

  dev.off()
}

# This function reads in the no. of days exceeding a KBDI threshold for a given epoch,
# Then aggregates each model by county first and then by state. Aggregation is defined as follows:
# Multiply each county's annual mean by the WUI population of that county, and add up across the counties in the state.
aggregate.by.state <- function(threshold = 600, start.date = "01011991", end.date = "12312010") {
  basename <- paste(sep = "-", "/kbdi-over", threshold, substr(start.date, 5, 8), substr(end.date, 5, 8))
  kbdi.days <- read.csv(paste0(wildfire.dir, basename, ".csv"))

  # First, we'll aggregate by county, multiplying our metric by WUI population

  # Get mean metric per county:
  county.means <- group_by(kbdi.days[,-c(1,2,4)], GEOID) %>% summarise_each(funs(mean))

  # Copy over data from some counties to counties with no grid points, so we have all US counties:
  for (i in 1:length(missing.counties)) {
    county.means <- rbind(county.means,
                          cbind(data.frame(GEOID = missing.counties[i]),
                                county.means[county.means$GEOID == replacement.counties[i], 2:ncol(county.means)]))
  }

  # Divide by 20 to look at annual mean:
  annual.means <- county.means[2:ncol(county.means)] / 20

  # Get population per counties of interest:
  pop <- read.csv("/media/eitan/My Book/wui.cty.all.csv")
  countypop <- data.frame(pop = pop$WUIPOP, row.names = pop$ID)
  county.pop <- countypop[as.character(county.means$GEOID),]

  # Multiply by county WUI population and convert GEOID to state ID:
  county.sums <- cbind(as.integer(county.means$GEOID / 1000), county.pop * annual.means)
  names(county.sums)[1] = "State"

  # Now aggregate by state:
  state.sums <- group_by(county.sums, State) %>% summarise_each(funs(sum))
  statepop <- read.csv("/media/eitan/My Book/statepop.csv")

  # Create final table for saving:
  final <- data.frame(STFIPS = state.sums$State, POP2010 = statepop$POP10)
  final <- cbind(final, state.sums[2:ncol(state.sums)])
  final <- cbind(final, t(apply(state.sums[2:ncol(state.sums)], 1, quantile)))

  outname <- paste(sep = "-", "/aggregated-kbdi-over", threshold, substr(start.date, 5, 8), substr(end.date, 5, 8))
  write.csv(final, paste0(wildfire.dir, outname), row.names = FALSE)
}
