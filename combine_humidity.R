library(doParallel);
registerDoParallel(cores=6);
root_dir <- "/fast/mtclimdata/"

satpres <- function(t) {
	es <- 100*6.11*(10 ^ ((7.5*t)/(237.3+t)))
	return(es)
}



onefile <- function(fn) {
	#get lat and lon
	line1 <- readLines(fn,1);
	lat <- strsplit(line1, " ")[[1]][6]
	lon <- strsplit(line1, " ")[[1]][8]
	pt <- read.table(file=fn,skip=4,stringsAsFactors = FALSE);
	columns <- c("YR","DAY","TMAX","TMIN","TDAY","PRECIP","VP","SRAD","DAYLEN");
	colnames(pt) <- columns;
	temps <- pt$TDAY;
	sats <- satpres(temps);
	vp <- pt$VP;
	relhum <- as.numeric(vp) / as.numeric(sats);
	relhum <- as.numeric(round(relhum,digits=4));
	row <- c(lat,lon,relhum)
	return(row)
}



go <- function(model, year) {
#	load("column.names.Rdata");
  setwd(model)
  setwd(year)
	relhum <- c()
	out <- c()
	i <- 1
  all.files <- list.files(pattern='.mtc43')
	for (out.file in all.files) {
	  relhum <- rbind(relhum,onefile(out.file))
		i <- i + 1
		if (i%%512 == 0) {
		  print(paste("Processed", round(100*i/length(all.files),1), "pct at", Sys.time()))
		  out <- rbind(out,relhum);
			relhum <- c()
		}

	}
	relhum <- as.data.frame(relhum)
	#colnames(relhum) <- names
	out <- as.data.frame(out)
	print(dim(out))
#	colnames(out) <- names
	relhum <- rbind(out,relhum)
	write.csv(relhum,paste0("humidity.csv"), row.names = FALSE)
  setwd("../..")
}

transform_all <- function() {
  models <- c("access1-0_rcp85_r1i1p1", "bcc-csm1-1-m_rcp85_r1i1p1", "bcc-csm1-1_rcp85_r1i1p1", "canesm2_rcp85_r1i1p1",
              "ccsm4_rcp85_r1i1p1", "cesm1-bgc_rcp85_r1i1p1", "cesm1-cam5_rcp85_r1i1p1", "cmcc-cm_rcp85_r1i1p1",
              "cnrm-cm5_rcp85_r1i1p1", "csiro-mk3-6-0_rcp85_r1i1p1", "fgoals-g2_rcp85_r1i1p1", "fio-esm_rcp85_r1i1p1",
              "gfdl-cm3_rcp85_r1i1p1", "gfdl-esm2g_rcp85_r1i1p1", "gfdl-esm2m_rcp85_r1i1p1", "giss-e2-r_rcp85_r1i1p1",
              "hadgem2-ao_rcp85_r1i1p1", "hadgem2-cc_rcp85_r1i1p1", "hadgem2-es_rcp85_r1i1p1", "inmcm4_rcp85_r1i1p1",
              "ipsl-cm5a-mr_rcp85_r1i1p1", "ipsl-cm5b-lr_rcp85_r1i1p1", "miroc-esm-chem_rcp85_r1i1p1", "miroc-esm_rcp85_r1i1p1",
              "miroc5_rcp85_r1i1p1", "mpi-esm-lr_rcp85_r1i1p1", "mpi-esm-mr_rcp85_r1i1p1", "mri-cgcm3_rcp85_r1i1p1", "noresm1-m_rcp85_r1i1p1")
  years <- c("2021", "2041")
  setwd(root_dir)
  
  foreach (m = 17:29, .combine=c) %dopar% {
  model = models[m]
    foreach (y = 1:length(years), .combine=c) %dopar% {
      year = years[y]
      print (paste("Working on model", model, "starting year:", year))
      go(model, year)
    }
  }
}