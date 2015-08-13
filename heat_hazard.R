# Combine humidity and teperature to come up with a heat risk category
library(stringr)
library(dplyr)

# Preliminaries: set globals and read in geo grid and heat index lookup table
root_dir <- "/fast"
models <- c("access1-0_rcp85_r1i1p1", "bcc-csm1-1-m_rcp85_r1i1p1", "bcc-csm1-1_rcp85_r1i1p1", "canesm2_rcp85_r1i1p1",
            "ccsm4_rcp85_r1i1p1", "cesm1-bgc_rcp85_r1i1p1", "cesm1-cam5_rcp85_r1i1p1", "cmcc-cm_rcp85_r1i1p1",
            "cnrm-cm5_rcp85_r1i1p1", "csiro-mk3-6-0_rcp85_r1i1p1", "fgoals-g2_rcp85_r1i1p1", "fio-esm_rcp85_r1i1p1",
            "gfdl-cm3_rcp85_r1i1p1", "gfdl-esm2g_rcp85_r1i1p1", "gfdl-esm2m_rcp85_r1i1p1", "giss-e2-r_rcp85_r1i1p1",
            "hadgem2-ao_rcp85_r1i1p1", "hadgem2-cc_rcp85_r1i1p1", "hadgem2-es_rcp85_r1i1p1", "inmcm4_rcp85_r1i1p1",
            "ipsl-cm5a-mr_rcp85_r1i1p1", "ipsl-cm5b-lr_rcp85_r1i1p1", "miroc-esm-chem_rcp85_r1i1p1", "miroc-esm_rcp85_r1i1p1",
            "miroc5_rcp85_r1i1p1", "mpi-esm-lr_rcp85_r1i1p1", "mpi-esm-mr_rcp85_r1i1p1", "mri-cgcm3_rcp85_r1i1p1", "noresm1-m_rcp85_r1i1p1")
years <- c("2021", "2041")

grid <- read.csv(paste(root_dir, "grid.classify.csv", sep = "/"))
danger_index <- read.csv(paste(root_dir, "heat_index_celsius.csv", sep = "/"), skip = 1)
rownames(danger_index) <- danger_index$Humidity / 100

# US counties that have no grid points:
missing.counties <- c(6075, 25019, 29510, 44001, 53055, 51720, 51520, 51640, 51750, 51775, 51580, 51530, 51790, 51660,
                      51540, 51600, 51013, 51683, 51685, 51630, 51670, 51570, 51830, 51735, 51595, 51610, 51678, 51690, 51710)
# The nearest respective counties to use for the missing.counties
replacement.counties <- c(6081, 25007, 29189, 44005, 53029, 51195, 51191, 51077, 51121, 51161, 51005, 51163, 51015,
                          51165, 51003, 51059, 51059, 51153, 51153, 51177, 51149, 51041, 51095, 51650, 51081, 51059,
                          51163, 51089, 51740)

statepop <- read.csv("/media/eitan/My Book/statepop.csv")

# Read in humidity data from current directory, merge it with geo and sort it by lon/lat:
read_humidity <- function() {
  print ("Reading humidities...")
  hum <- read.csv("humidity.csv", check.names = FALSE)
  names(hum)[3:7307] = str_pad(names(hum)[3:7307], 8, pad = "0")
  hum <- merge(hum, grid[,1:2])
  hum <- hum[with(hum, order(LON, LAT)),]
  print ("Done!")
  hum
}

# Read in humidity data from current directory. Ensure it has exactly the same geo points
# and ordering as the humidity data, for index matching
read_temps <- function(hum) {
  print ("Reading temperatures...")
  load("tmax.Rdata")
  tmax <- merge(tmax, hum[,1:2])
  tmax <- tmax[with(tmax, order(LON, LAT)),]
  print ("Done!")
  tmax
}

# Compute a level of heat danger for a given humidity an temperature pair.
# Humidity is first rounded to the nearst 5% value to match lookup table.
heat_index <- function(humidity, temp) {
  h <- as.character(floor(humidity * 20) / 20)
  if (temp < danger_index[h,"T1"])  return(0)
  if (temp < danger_index[h,"T2"])  return(1)
  if (temp < danger_index[h,"T3"])  return(2)
  if (temp < danger_index[h,"T4"])  return(3)
  return (4)
}

# Given a row with humidity values and a row with temperature values, compute
# row with the appropriate heat index values
compute_row <- function(hums, temps) {
  ret <- hums
  for (j in 3:ncol(temps)) {
    ret[,j] <- heat_index(hums[,j], temps[,j])
  }
  ret
}

# Compute the head index for an entire model/year set of coordinates and days
process_model <- function(model, year) {
  print (paste("Working on model:", model, "year:", year, "at:", Sys.time()))
  setwd(paste(root_dir, "mtclimdata", model, year, sep = "/"))
  hum <- read_humidity()
  tmax <- read_temps(hum)
  write.csv(hum[,1:7307], "sorted-hum.csv", row.names = FALSE)
  write.csv(tmax[,1:7307], "sorted-tmax.csv", row.names = FALSE)
  
  # Return matrix starts from hum, replacing all humidities with heat index
#  ret <- data.frame()
#  for (i in 1:nrow(hum)) {
#    ret <- rbind(ret, compute_row(hum[i,], tmax[i,]))
#  }
}


# This function receives a vector of metrics (with GEOID for keys) and aggregates the metric
# by county and then by state. Aggregation is defined as follows:
# Multiply each county's annual mean by the population of that county, and add up across the counties in the state.
aggregate.by.state <- function(metrics) {
  # Get mean metric per county:
  data <- data.frame(ID = names(metrics), value = metrics / 20)  # Divide by 20 to get annual mean
  grouped <- group_by(data, ID) %>% summarise(mean.value = mean(value))
  county.means <- data.frame(ID = as.numeric(as.character(grouped$ID)),
                             county.mean = grouped$mean.value)

  # Copy over data from some counties to counties with no grid points, so we have all US counties:
  for (i in 1:length(missing.counties)) {
    county.means <- rbind(county.means,
                          c(missing.counties[i], county.means[county.means$ID == replacement.counties[i], 2]))
  }

  # Get population per counties of interest:
  pop <- read.csv("/media/eitan/My Book/countypov.csv")
  countypop <- data.frame(pop = pop$vul.pop, row.names = pop$GEOID)
  county.pop <- countypop[as.character(county.means$ID),]

  # Multiply by county population and convert FIPS ID to state ID:
  county.sums <- data.frame(STFIPS = as.integer(county.means$ID / 1000),
                            county.sum = county.pop * county.means$county.mean)

  # Now aggregate by state:
  state.sums <- group_by(county.sums, STFIPS) %>% summarise(state.sum = sum(county.sum))
  return(state.sums)
}

# For a given model and start year, return the weighted no. of days that exceed a heat index
# threshold, averaged per year and aggregated by states (see aggregate.by.state)
aggregate.by.model <- function(model, startyr, threshold = 3) {
  fn = paste(sep="/", "/fast/mtclimdata", model, startyr, "heat_index.csv")
  hi <- read.csv(fn, check.names = FALSE)
  gridded <- merge(grid[,1:3], hi, by=c("LON", "LAT"))
  high.days <- apply(gridded[,4:ncol(gridded)], 1, function(row) { sum(row >= threshold)})
  names(high.days) = gridded$GEOID
  return(aggregate.by.state(high.days))
}

# Output a summary table of the weighted no. of danger days for the observed period (1991-2010)
output.observed.aggregates <- function() {
  df <- statepop[,c(1,4)]
  means <- aggregate.by.model("baseline", "1991")
  df <- merge(df, means)
  names(df)[ncol(df)] = "baseline"
  write.csv(df, paste0(root_dir, "/mtclimdata/aggregated-1991.csv"), row.names = FALSE)
}

# Output a summary table of the weighted no. of danger days for future periods, with quantiles
output.future.aggregates <- function(startyr = "2021") {
  df <- statepop[,c(1,4)]
  for (m in models) {
    print(paste("Working on model", m, Sys.time()))
    df <- merge(df, aggregate.by.model(m, startyr))
    names(df)[ncol(df)] = m
  }
  df <- cbind(df, t(apply(df[3:ncol(df)], 1, quantile)))
  write.csv(df, paste0(root_dir, "/mtclimdata/aggregated-", startyr, ".csv"), row.names = FALSE)
}
