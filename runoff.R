# Compute total water runoff in future (modeled) periods based on 95th-percentile data modeled on the past.

library(dplyr)

root.dir <- "/data/runoff"

models <- c("access1-0_rcp85_r1i1p1", "bcc-csm1-1-m_rcp85_r1i1p1", "bcc-csm1-1_rcp85_r1i1p1", "canesm2_rcp85_r1i1p1",
            "ccsm4_rcp85_r1i1p1", "cesm1-bgc_rcp85_r1i1p1", "cesm1-cam5_rcp85_r1i1p1", "cmcc-cm_rcp85_r1i1p1",
            "cnrm-cm5_rcp85_r1i1p1", "csiro-mk3-6-0_rcp85_r1i1p1", "fgoals-g2_rcp85_r1i1p1", "fio-esm_rcp85_r1i1p1",
            "gfdl-cm3_rcp85_r1i1p1", "gfdl-esm2g_rcp85_r1i1p1", "gfdl-esm2m_rcp85_r1i1p1", "giss-e2-r_rcp85_r1i1p1",
            "hadgem2-ao_rcp85_r1i1p1", "hadgem2-cc_rcp85_r1i1p1", "hadgem2-es_rcp85_r1i1p1", "inmcm4_rcp85_r1i1p1",
            "ipsl-cm5a-mr_rcp85_r1i1p1", "ipsl-cm5b-lr_rcp85_r1i1p1", "miroc-esm-chem_rcp85_r1i1p1", "miroc-esm_rcp85_r1i1p1",
            "miroc5_rcp85_r1i1p1", "mpi-esm-lr_rcp85_r1i1p1", "mpi-esm-mr_rcp85_r1i1p1", "mri-cgcm3_rcp85_r1i1p1", "noresm1-m_rcp85_r1i1p1")


# Given a date string of the foramt "MMDDYYYY" return a 4-letter string representing
# the water year, running from Oct. 1st of the previous year to Sep. 30th of this one.
water.year <- function(date) {
  yr <- substr(date, 5, 8)
  if (substr(date, 1, 2) < "10") {
    return (yr)
  } else {
    return (as.character(as.numeric(yr) + 1))
  }
}

# For a given model, compute the total future runoff per water year that exceeds a threshold.
# The threshold is the 95th percentile value of the past period.
total_runoff <- function(model) {
  # Read in modeled data water years 1991-2010, 2021-2040, 2041-2060:
  print(paste("Reading in runoff files for model", model, "at time", Sys.time()))
  runoff.91 <- read.csv(paste0(root.dir, "/total_runoff.", model, ".1991.csv"), check.names = FALSE)
  runoff.21 <- read.csv(paste0(root.dir, "/total_runoff.", model, ".2021.csv"), check.names = FALSE)
  runoff.41 <- read.csv(paste0(root.dir, "/total_runoff.", model, ".2041.csv"), check.names = FALSE)

  # Save first five columns and delete them, to compare only runoff values:
  coords <- runoff.91[,1:5]
  runoff.91 <- runoff.91[,6:ncol(runoff.91)]
  runoff.21 <- runoff.21[,6:ncol(runoff.21)]
  runoff.41 <- runoff.41[,6:ncol(runoff.41)]
  
  # Compute 95th quantile of runoff in baseline period, per grid point:
  print(paste("Copmuting thresholds at time", Sys.time()))
  thresholds <- apply(runoff.91, 1,
                      function(row) { quantile(as.numeric(row), 0.95, na.rm = TRUE)})

  # For each grid point in all periods, compute total runoff that exceeds threshold:
  print(paste("Computing excess runoff at time", Sys.time()))
  excess.91 <- numeric(nrow(runoff.91))
  excess.21 <- numeric(nrow(runoff.21))
  excess.41 <- numeric(nrow(runoff.41))
  
  for (i in 1:length(thresholds)) {
    excess.91[i] = sum(ifelse(as.numeric(runoff.91[i,]) > thresholds[i], as.numeric(runoff.91[i,]) - thresholds[i], 0))
    excess.21[i] = sum(ifelse(as.numeric(runoff.21[i,]) > thresholds[i], as.numeric(runoff.21[i,]) - thresholds[i], 0))
    excess.41[i] = sum(ifelse(as.numeric(runoff.41[i,]) > thresholds[i], as.numeric(runoff.41[i,]) - thresholds[i], 0))
    if (i %% 525 == 0) {
      print(paste("Completed", round(100 * i / length(thresholds), 1), "percent at time", Sys.time()))
    }
  }
  
  save.91 <- coords
  save.91$total_runoff <- excess.91
  save.21 <- coords
  save.21$total_runoff <- excess.21
  save.41 <- coords
  save.41$total_runoff <- excess.41
  write.csv(save.91, file = paste0(root.dir, "/tmp_save_excess_runoff.", model, ".1991.csv"), row.names = FALSE)
  write.csv(save.21, file = paste0(root.dir, "/tmp_save_excess_runoff.", model, ".2021.csv"), row.names = FALSE)
  write.csv(save.41, file = paste0(root.dir, "/tmp_save_excess_runoff.", model, ".2041.csv"), row.names = FALSE)

  print(paste("Done at time", Sys.time()))
}

# Temporary crutch to move from a bad grid merge to a better one
replace.grid.data <- function(grid, model) {
  for (yr in c(1991, 2021, 2041)) {
    fn <- paste0(root.dir, "/tmp_save_excess_runoff.", model, ".", yr, ".csv")
    orig <- read.csv(fn)
    merged <- merge(grid[,c(1,2,3,7,8)], orig[,c(1,2,6)])
    merged$state = as.integer(merged$GEOID / 1000)
    write.csv(merged, file = fn, row.names = FALSE)
  }
}


# (Obsolote) function to aggregate data by watershed and state, ignoring population
average.by.watershed.state <- function(model) {
  for (year in c(1991, 2021, 2041)) {
    print (paste("Processing model:", model, "year:", year))
    data <- read.csv(paste0(root.dir, "/tmp_save_excess_runoff.", model, ".", year, ".csv"))
    by.state.huc <- group_by(data, HUC8, state) %>% summarise(mean.runoff = mean(total_runoff) / 20) # Twenty years
    names(by.state.huc) <- c("HUC8", "State", paste0("epoch-", year))
    if (year == 1991) {
      mdata <- by.state.huc
    } else {
      mdata <- merge(mdata, by.state.huc)
    }
  }
  return (mdata)
}


plot.runoff.map <- function(year) {
  library(ggplot2)
  library(sp)
  library(maptools)
  library(plyr)
  library(maps)

  runoff <- read.csv(paste0(root.dir, "/median_excess_runoff.", year, ".csv"))
  map.counties <- map_data("county")
  map.states <- map_data("state")

  counties <- map('county', fill = TRUE, col = "transparent", plot = FALSE)
  IDs <- sapply(strsplit(counties$names, ":"), function(x) x[1])
  counties_sp <- map2SpatialPolygons(counties, IDs = IDs,
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))
  pointsSP <- SpatialPoints(as.data.frame(runoff[,c('LON','LAT')]),
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  indices <- over(pointsSP, counties_sp)

  county.names <- sapply(counties_sp@polygons, function(x) x@ID)
  runoff$county<-county.names[indices]
  mapcounties <- map_data("county")
  mapstates <- map_data("state")
  mapcounties$county <- with(mapcounties, paste(region, subregion, sep = ","))
  runoff.mean <- ddply(runoff, c("county"), summarize, means = mean(median.excess.runoff))
  runoff.un <- runoff[!duplicated(runoff[c("county")]),]
  runoff.final <- merge(runoff.mean, runoff.un, by = 'county')
  breaks <- c(1, 100,500,1000,2000,5000,10000)
  runoff.final$colorBuckets <- as.factor(as.numeric(cut(runoff.mean$means, breaks)))
  mergedata <- merge(mapcounties, runoff.final, by = "county")
  mergedata <- mergedata[order(mergedata$means),]

  pdf(paste0(root.dir, "/county-excess-runoff-", year, ".pdf"))
  ggplot(mergedata, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = colorBuckets)) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_fill_brewer(palette="PuRd", labels = breaks, name = "Excess runoff (mm/m^2)") +
    coord_map(project="globular") +
    ggtitle(paste0("Median excess runoff ", year, "-", year + 19, ", averaged per county")) +
    geom_path(data = map.states, colour = "black", size = .3) +
    geom_path(data = mapcounties, colour = "white", size = .5, alpha = .1)
  dev.off()
}

# For each state, aggregate excess runoff by watershed, and average weighted by population
compute.state.runoff <- function(watersheds, model) {
  setwd(root.dir)
  model.91 <- read.csv(paste0("tmp_save_excess_runoff.", model, ".1991.csv"))
  model.21 <- read.csv(paste0("tmp_save_excess_runoff.", model, ".2021.csv"))
  model.41 <- read.csv(paste0("tmp_save_excess_runoff.", model, ".2041.csv"))

  runoff.91 <- group_by(model.91, HUC8, state) %>% summarise(annual.mean.runoff = mean(total_runoff) / 20)
  runoff.21 <- group_by(model.21, HUC8, state) %>% summarise(annual.mean.runoff = mean(total_runoff) / 20)
  runoff.41 <- group_by(model.41, HUC8, state) %>% summarise(annual.mean.runoff = mean(total_runoff) / 20)

  merged.91 <- merge(watersheds, runoff.91)
  merged.21 <- merge(watersheds, runoff.21)
  merged.41 <- merge(watersheds, runoff.41)

  agg.91 <- group_by(merged.91, state) %>% summarise(wrunoff = sum(annual.mean.runoff * POP), pop = sum(POP), n = n())
  agg.21 <- group_by(merged.21, state) %>% summarise(wrunoff = sum(annual.mean.runoff * POP), pop = sum(POP), n = n())
  agg.41 <- group_by(merged.41, state) %>% summarise(wrunoff = sum(annual.mean.runoff * POP), pop = sum(POP), n = n())

  ret <- data.frame(state = agg.91$state,
                    weighted.runoff.1991 = agg.91$wrunoff / agg.91$pop,
                    weighted.runoff.2021 = agg.21$wrunoff / agg.21$pop,
                    weighted.runoff.2041 = agg.41$wrunoff / agg.41$pop
  )

  write.csv(ret, paste0("weighted_runoff_summary.", model, ".csv"), row.names = FALSE)
  ret
}

# For a given list of data for one state (each list item containing aggregated
# excess runoff for one model and three epochs), compute the median change from
# the 1991 epoch to the 2021 and 2041. Median change is only computed if 20 or
# more of the models agree on the direction (sign) of the change. Otherwise,
# it is defined NA.
median.change <- function(state.data) {
  ratios <- sapply(state.data, function (row)
    (row$weighted.runoff.2021 - row$weighted.runoff.1991) / row$weighted.runoff.1991)
  if (sum(ratios > 0) >= 20) {
    med.21 <-  median(ratios[ratios > 0])
  } else if (sum(ratios < 0) >= 20) {
    med.21 <-  median(ratios[ratios < 0])
  } else {
    med.21 <- NA
  }

  ratios <- sapply(state.data, function (row)
    (row$weighted.runoff.2041 - row$weighted.runoff.1991) / row$weighted.runoff.1991)
  if (sum(ratios > 0) >= 20) {
    med.41 <-  median(ratios[ratios > 0])
  } else if (sum(ratios < 0) >= 20) {
    med.41 <-  median(ratios[ratios < 0])
  } else {
    med.41 <- NA
  }

  return(c(round(state.data[[1]]$state, 0), round(med.21, 4), round(med.41, 4)))
}

# For all states, compute the median change from 2021 to 1991 and 1991 to 2041,
# if a large majority of models agree on the sign of the change.
median.runoff.agreement <- function(models) {
  files <- lapply(models, function(m) paste0("weighted_runoff_summary.", m, ".csv"))
  data <- lapply(files, read.csv)
  ret <- data.frame(state = numeric(), median.2021 = numeric(), median.2041 = numeric())

  for (i in 1:nrow(data[[1]])) {
    ret <- rbind(ret, median.change(lapply(data, "[", i,)))
  }

  names(ret) = c("state", "median.change.ratio.2021", "median.change.ratio.2041")
  return(ret)
}

for (model in models) {
  total_runoff(model)
}

watersheds <- read.csv(paste0(root.dir, "/huc8pops.csv"))
watersheds <- rename(watersheds, state = STFIPS)

for (model in models) {
  compute.state.runoff(watersheds, model)
}

final.results <- median.runoff.agreement(models)
write.csv(final.results, "weighted_runoff_summary.csv", row.names = FALSE)
