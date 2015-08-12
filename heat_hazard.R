# Combine humidity and teperature to come up with a heat risk category
library(stringr)

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

system.time(
for (model in models) {
  for (year in years) {
    process_model(model, year)
  }
}
)