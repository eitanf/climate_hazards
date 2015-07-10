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
    merged$state = substr(merged$GEOID, 1, 2)
    write.csv(merged, file = fn, row.names = FALSE)
  }
}

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

for (model in models) {
  total_runoff(model)
}

for (model in models) {
  mdata <- average.by.watershed.state(model)
  write.csv(mdata, paste0(root.dir, "/annual.means.", model, ".csv", row.names = FALSE))
}