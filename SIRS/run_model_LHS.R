# Script to run SIR model with waning with Latin Hypercube Sampling of parameters

# load packages
library("FME")
library("doParallel")
library("dplyr")

# change to false if don't want to plot or save individual simulation results
save_plots = TRUE
save_output_file = TRUE

# set number of cores to use for parallel processing
cl = makeCluster(6, outfile="")
registerDoParallel(cl)

# read default parameters file
parms = read.table(file="parameters.dat",header=F,row.names=1)
parms = setNames(as.numeric(t(parms)),rownames(parms))

# parameters to vary and their ranges
log_base = 10
parms_to_vary = t(data.frame(imm_min=c(0,0.75),
                            imm_shape=c(log(0.5,base=log_base),log(10,base=log_base)),
                            imm_t=c(25,250)))
nparms = nrow(parms_to_vary)
# names of parameters that are being varied for plotting
parm_names = c("Long maintained immunity",
                "Steepness of waning immunity function",
                "Mean time immune")
# whather the parameter is being sampled over a log scale or not
parms_log = c(FALSE,
              TRUE,
              FALSE)

# number of simulations to run
number_simulations = 10

# sample parameters using Latin Hypercube sampling
parms_LHS = Latinhyper(parms_to_vary,num=number_simulations)
# convert log scale parameters back
parms_LHS[,parms_log] = log_base^parms_LHS[,parms_log]

# list to save simulation results
simulations_out = vector(mode='list', length=number_simulations)

# create directory to save files
# name directory with timestamp
my_dir = paste0(getwd(),"/output/",format(Sys.time(), "%y%m%d_%H%M%S"))
dir.create(my_dir,recursive=TRUE)
if (save_plots | save_output_file){
  dir.create(paste0(my_dir,"/simulations"))

  # copy parameter file to output directory
  file.copy("parameters.dat",paste0(my_dir,"/parameters.dat"))
}

# run simulations in parallel
simulations_out = foreach(i = 1:number_simulations) %dopar% {
  # load required package
  suppressWarnings(suppressMessages(require(zoo)))

  # load model functions
  source("peaks_function.R")
  source("SIRS_with_waning.R")

  # if saving plots or output files create subdirectory
  if (save_plots | save_output_file){
    # create directory for current simulation run
    dir_sub = paste0(my_dir,"/simulations/",i)
    dir.create(dir_sub)

    # save parameter file to simulation subdirectory
    write.table(parms,file=paste0(dir_sub,"/parmeters.dat"),quote=FALSE,col.names=FALSE,sep="        ")
  }

  # loop over parameters we are varying and change their value
  for (p in 1:nparms){
    parms[row.names(parms_to_vary)[p]] = parms_LHS[i,p]
  }

  ##### Run simulation
  sim = SIRS_with_waning(parms,save_plots,save_output_file,dir_sub)

  # get the peaks in the infection curve
  peaks = sim$peaks
  n_peaks = length(peaks)

  # time elapsed since the previous peak
  peaks_dt = rep(NA,n_peaks)
  if (n_peaks>1) peaks_dt[2:n_peaks] = peaks[2:n_peaks] - peaks[1:(n_peaks-1)]

  # return simulation output
  if (n_peaks > 0){
    out = as.data.frame(t(matrix(rep(parms,n_peaks),length(parms),n_peaks)))
    names(out) = names(parms)
    out$peaks = peaks
    out$peaks_dt = peaks_dt
    out$simulation_id = rep(i,n_peaks)
    out$simulation_dir = rep(my_dir,n_peaks)
    out
  }
}

# convert list of individual simulations output to data frame
all_peaks = bind_rows(simulations_out)

# save LHS outputs (the time of the peaks in the infection curve for each simulation)
write.csv(all_peaks,file=paste0(my_dir,"/LHS_output.csv"),quote=FALSE,row.names=FALSE)
