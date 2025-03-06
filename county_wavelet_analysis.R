#-----------------------------------------------------------#
# Run a wavelet analysis of NY Times Covid-19 cases and     #
# correlate wavelets to climatic and demographic variables. #
# Author: Ilan N. Rubin                                     #
# Date: February 27, 2025                                   #
#-----------------------------------------------------------#

# load packages
library("dplyr");library("tidyr")
library("ggplot2");library("RColorBrewer");library("ggtext")
library("stats")
library("jsonlite")
library("tigris");library("sf")
library("FME")
library("doParallel")
library("WaveletComp") # wavelet analysis
# library("zoo") # for peak algorithm

# set working directory
setwd("~/Documents/Periodicity/code/wavelets")

# load the wavelet analysis script
source("peaks_function.R")
source("wavelet_function.R")
source("plot_function.R")

# load scripts to import case data
source("nytimes_case_data.R")

# import google mobility data
source("google_mobility_data.R")

# create output directory
if (!dir.exists("output/states")) dir.create("output/states",recursive=TRUE)

# loop over fips
no_county = nytimes_case_data$county=="Unknown" # cut out cases that are not assigned to a specific county
all_fips = sort(unique(nytimes_case_data$fips[!no_county]))
n_fips = length(all_fips) # unique fips to run wavelet analysis

# initiate parallel processing
cl = makeCluster(7, outfile="")
registerDoParallel(cl)

# create simulation log
if (!dir.exists("logs")) dir.create("logs")
log_file = paste0("logs/log_",format(Sys.time(), "%Y%m%d_%H%M"),".txt")
file.create(log_file)
write(paste(format(Sys.time(), "%Y-%m-%d %H:%M"),"Started"),log_file,append=TRUE)

# function to return output from wavelet even if it is NA
return_output =  function(x,output) if(any(is.na(x))){x}else{x[[output]]}

# Loop over all fips in parallel to run wavelet analysis
wavelet_out = foreach(cnt = 1:n_fips, .inorder=FALSE, .errorhandling="pass", .verbose=FALSE) %dopar% {

  # packages and scripts need to be reloaded inside parallel loop
  #load dplyr, waveletcomp, and zoo
  r = require("dplyr")
  r = require("WaveletComp")
  r = require("zoo")

  # load the wavelet analysis script
  source("wavelet_function.R")

  # extract cases for just this region
  this_fips = all_fips[cnt]
  temp_cases = nytimes_case_data[nytimes_case_data$fips==this_fips,]

  # save state and county
  this_state = temp_cases$state[1]
  this_county = temp_cases$county[1]

  # update console and log file
  print(paste(this_state,"-",this_county))
  write(paste(format(Sys.time(), "%Y-%m-%d %H:%M"),":",cnt,this_state,this_county),log_file,append=TRUE)

  # get county mobility data
  temp_mobility = mobility[mobility$fips==this_fips,]

  # create directory for plots if it does not already exist and set directory
  if (!dir.exists(paste0("output/states/",this_state))) dir.create(paste0("output/states/",this_state))
  if (!dir.exists(paste0("output/states/",this_state,"/",this_county))) dir.create(paste0("output/states/",this_state,"/",this_county))
  dir = paste0("output/states/",this_state,"/",this_county)

  # plot case counts
  plot_cases(temp_cases,dir)
  plot_cases(temp_cases,dir,column="log_cases_100k")

  # only run wavelet is the county had recorded cases
  if (sum(temp_cases$cases) > 1){
    # run wavelet analysis
    wavelet_log = run_wavelet(temp_cases,"log_cases_100k",directory=dir,identifier="log_cases",save_plots=TRUE)

    # get the wavelet ridge for plotting
    ridge_log = which(wavelet_log$Ridge==1,arr.ind=TRUE)

    # find peaks in the spectrum
    peaks_log = findpeaks(wavelet_log$Period, wavelet_log$Power.avg, prom=0.01, w=1)
  }else{
    wavelet_log = NA
    ridge_log = NA
    peaks_log = NA
  }

  # run mobility wavelets as long as they have at least 500 days of data
  if (sum(!is.na(temp_mobility[mobility_column[1]]))>500){
    wavelet_workplace = run_wavelet(temp_mobility,mobility_column[1],identifier="workplace",directory=dir,save_plots=FALSE)
  }else{
    wavelet_workplace = NA
  }
  if (sum(!is.na(temp_mobility[mobility_column[2]]))>500){
    wavelet_retail = run_wavelet(temp_mobility,mobility_column[2],identifier="retail",directory=dir,save_plots=FALSE)
  }else{
    wavelet_retail = NA
  }

  # return county output as a list
  list(cases = temp_cases$cases_avg_per_100k,
        date = temp_cases$date,
        wavelet_log = return_output(wavelet_log,"Power.avg"),
        wavelet_date = return_output(wavelet_log,"Period"),
        ridge_log = ridge_log,
        peaks_log = peaks_log$x,
        wavelet_workplace = return_output(wavelet_workplace,"Power.avg"),
        wavelet_retail = return_output(wavelet_retail,"Power.avg"),
        region = temp_cases[1,c("state","county","fips")])
}
write(paste(format(Sys.time(), "%Y-%m-%d %H:%M"),"Finished"),log_file,append=TRUE)

# stop parallel processing
stopCluster(cl)

# convert output list from parallel processing to dataframes
all_region =  bind_rows(lapply(wavelet_out,function(x){ if(is.null(x$region)){c(state=NA,county=NA,fips=NA)}else{ x$region}}))
wavelet_date =  wavelet_out[[1]]$wavelet_date
all_cases =  bind_rows(lapply(wavelet_out,function(x){n=length(x$date);data.frame(fips=rep(x$region$fips,n),state=rep(x$region$state,n),county=rep(x$region$county,n),x[c("date","cases")])}))
all_wavelet_log =  t(sapply(wavelet_out,function(x){ if(all(is.na(x$wavelet_log))){rep(NA,length(wavelet_date))}else{x$wavelet_log}}))
all_peaks_log =  bind_rows(lapply(wavelet_out,function(x){n=length(x$peaks_log);data.frame(fips=rep(x$region$fips,n),state=rep(x$region$state,n),county=rep(x$region$county,n),x["peaks_log"])}))
all_wavelet_workplace =  t(sapply(wavelet_out,function(x){ if(all(is.na(x$wavelet_workplace))){rep(NA,length(wavelet_date))}else{x$wavelet_workplace}}))
all_wavelet_retail =  t(sapply(wavelet_out,function(x){ if(all(is.na(x$wavelet_retail))){rep(NA,length(wavelet_date))}else{x$wavelet_retail}}))
all_ridge_log =  bind_rows(lapply(wavelet_out,function(x){n=dim(x$ridge_log)[1];if(!is.null(n)){data.frame(fips=rep(x$region$fips,n),state=rep(x$region$state,n),county=rep(x$region$county,n),x["ridge_log"])}}))

# write output to files
write.csv(all_cases,file="output/all_counties_cases.csv",row.names=FALSE,quote=FALSE)
write.csv(cbind(all_region,all_wavelet_log),file="output/all_counties_wavelet_log.csv",row.names=FALSE,quote=FALSE)
write.csv(cbind(all_region,all_wavelet_retail),file="output/all_counties_wavelet_retail.csv",row.names=FALSE,quote=FALSE)
write.csv(cbind(all_region,all_wavelet_workplace),file="output/all_counties_wavelet_workplace.csv",row.names=FALSE,quote=FALSE)
write.csv(all_ridge_log,file="output/all_counties_ridge_log.csv",row.names=FALSE,quote=FALSE)

# calculate annual periodicity of cases and mobility and write to file
wavelet_365 = sort(abs(wavelet_date-365),index.return=TRUE)$ix[1] # find the period closest to exactly 365 days
all_region$annual_periodicity = all_wavelet_log[,wavelet_365]
all_region$retail = all_wavelet_retail[,wavelet_365]
all_region$workplace = all_wavelet_workplace[,wavelet_365]
write.csv(all_region,file="output/all_counties_periodicity.csv",row.names=FALSE,quote=FALSE)





##### Plot results
# load data for correlates
# some of the data sources from NOAA and the Census are large and this can take a little while
source("climate_data.R")
source("election_data.R")
source("census_data.R")
source("mapping_data.R")

# get the total number of cases in each county
all_cases = summarize(group_by(nytimes_case_data,fips),total_cases=sum(cases))

# add all the data to the shape files
all_data = full_join(usgeo,all_region,by="fips")
all_data = left_join(all_data,election_2020[c("fips","perc_repub")],by="fips")
all_data = left_join(all_data,climate,by="fips")
all_data = left_join(all_data,census[c("fips","pop","area","perc_insured","perc_poverty","over_65","density")],by="fips")
names(HHS)[names(HHS)=="STATE_NAME"] = "state"
all_data = left_join(all_data,HHS,by="state")
all_data = left_join(all_data,all_cases,by="fips")
all_data$total_cases_100k = all_data$total_cases/all_data$pop * 100000

# plot annual periodicity of COVID-19 cases
p = ggplot(all_data,aes(fill=annual_periodicity)) + 
  geom_sf(color="#00000020") +
  scale_fill_gradient(low="lightgoldenrod1",high="orchid4") +
  theme_void() +
  theme(legend.position=c(0.075,0.8),legend.title=element_blank())
ggsave("output/map_annual_periodicity.pdf",p,height=4.5,width=7.5)

# add borders of the HHS regions to the map as well
p = p +
    geom_sf(data=usgeo_HHS,color="#000000",fill=NA,linewidth=1.5)
ggsave("output/map_annual_periodicity_with_HHS.pdf",p,height=4.5,width=7.5)

# list of varibales to plot, whether they are logscale or not, and their axis titles
var_to_plot = c("pop","area","density","perc_repub","max_temp","min_temp","avg_percip","diff_temp","perc_insured","perc_poverty","over_65","retail","workplace","total_cases_100k")
var_log =     c(TRUE, TRUE,  TRUE,     FALSE,       FALSE,     FALSE,     FALSE,       FALSE,      FALSE,         FALSE,         FALSE,    FALSE,   FALSE,      TRUE)
var_name =    c("Log population size","Log geographic area (sq miles)","Log population per sq mile","Republican vote share, 2020 presidential election","Maximum temperature (F)","Minimum temperature (F)","Average percipitation","Temperature variability (F)","% of population with health insurance","% of population under 200% of the poverty line","% of population 65+","Annual periodicity of Google mobility (retail)","Annual periodicity of google mobility (workplace)","Log total number of COVID-19 cases per 100k")

# Loop over all variables, map and plot
for (v in 1:length(var_to_plot)){
  var = var_to_plot[v]

  # map variable and save
  p = ggplot(all_data,aes(fill=.data[[var]])) + 
    geom_sf(color="#00000020") +
    theme_void() +
    theme(legend.position=c(0.075,0.8),legend.title=element_blank())
  p = p + ylab(var_name[v])
  # add logscale coloring if relevant
  if (var_log[v]){
    p = p + scale_fill_gradient(low="lightgoldenrod1",high="orchid4",trans = "log")
  }else{
    p = p + scale_fill_gradient(low="lightgoldenrod1",high="orchid4")
  }
  ggsave(paste0("output/map_",var,".pdf"),p,height=4.5,width=7.5)

  # Plot variable versus COVID-19 annual periodicity
  # only if that region had more than 500 cases to minimize weirdenss in the wavelets
  p = ggplot(filter(all_data,total_cases>500),aes(x=.data[[var]],y=annual_periodicity)) +
    geom_smooth(method="lm",formula="y~x",color="black",alpha=0.8) +
    geom_point(size=0.1,color="black") +
    theme_classic() 
  p = p + xlab(var_name[v]) + ylab("Annual periodicity of COVID-19 cases")
  # add log scale to axis if relevant
  if (var_log[v]){
    p = p + scale_x_log10()
  }
  ggsave(paste0("output/covid_periodicity_vs_",var,".pdf"),p,height=4.5,width=7.5)
}

# Plot boxplot for HHS regions versus annual periodicity
p = ggplot(filter(all_data,!is.na(HHS),total_cases>500), aes(x=HHS,y=annual_periodicity)) +
  geom_boxplot() +
  xlab("HHS region") + ylab("Annual periodicity of COVID-19 cases") +
  theme_classic()
ggsave("output/covid_periodicity_vs_HHS_region.pdf",p,height=4.5,width=7.5)




### Calculate correlations and regression
# take the log of any variables that are compared on a logscale
for (v in 1:length(var_to_plot)){
  if(var_log[v]){
    all_data[[var_to_plot[v]]] = log(all_data[[var_to_plot[v]]],10)
  }
}

# calculate multivariable regression of all correlatees except temperature variability and population density
# with annual COVID-19 periodicity
var_to_model = var_to_plot[var_to_plot!="diff_temp"&var_to_plot!="density"]
multivariate_model = lm(paste("annual_periodicity ~ ",paste(var_to_model,collapse="+")),data=na.omit(all_data))

# number of correlates
n_var = length(var_to_plot)

# data frame to save regression results
# regression_results = data.frame(var=var_to_model,
#                         univariable=rep("",n_var),
#                         univariable_r2=rep("",n_var),
#                         univariable_pearsons = rep("",n_var),
#                         multivariable_max_temp=rep("",n_var),
#                         multivariable_diff_temp=rep("",n_var))
regression_results = data.frame(var=var_to_plot,
                                univariable_coefs=rep("",n_var),
                                univariable_r2=rep("",n_var),
                                univariable_pearsons = rep("",n_var),
                                multivariable_coefs=rep("",n_var),
                                latex = rep("",n_var))


# loop over all variables and calculate univariable regression
# save univariable regression and multivariable regression results to table
# cut out places with 500 or fewer total cases to control for weirdness in the wavelet analysis
for (i in 1:n_var){
  v = var_to_plot[i]

  # calculate univariable regression
  temp_model = lm(paste("annual_periodicity ~ ",v),data=filter(all_data,total_cases>500))
  # save regression coefficients and 95% confidence intervals
  regression_results$univariable_coefs[i] = paste0(sprintf("%0.2e",temp_model$coefficients[2])," (",paste(sprintf("%0.2e",confint(temp_model)[2,]),collapse=", "),")")
  # save univariable regression r^2
  regression_results$univariable_r2[i] = sprintf("%0.2e",summary(temp_model)$r.squared)
  # save Pearson's Correlation Coefficient
  regression_results$univariable_pearsons[i] = sprintf("%0.2e",cor.test(all_data[[v]],all_data$annual_periodicity)$estimate)

  # save the coefficients and 95% confidence intervals for this variable from the multivariable model
  if (v %in% row.names(confint(multivariate_model))){
    w = which(row.names(confint(multivariate_model))==v)
    regression_results$multivariable_coefs[i] = paste0(sprintf("%0.2e",multivariate_model$coefficients[w])," (",paste(sprintf("%0.2e",confint(multivariate_model)[w,]),collapse=", "),")")
    multivariate_latex = paste0("$",sprintf("%0.2e",multivariate_model$coefficients[w]),"$ & $(",paste(sprintf("%0.2e",confint(multivariate_model)[w,]),collapse=", "),")$")
  }else{
    multivariate_latex=""
  }
  # if (v %in% row.names(confint(diff_model))){
  #   w = which(row.names(confint(diff_model))==v)
  #   regression_results[i,4] = paste0(sprintf("%0.2e",diff_model$coefficients[w])," (",paste(sprintf("%0.2e",confint(diff_model)[w,]),collapse=", "),")")
  # }

  regression_results$latex[i] = multivariate_latex
}
# add column for ANOVA for hhs regions
HHS_anova = lm(annual_periodicity ~ HHS,data=filter(all_data,!is.na(HHS)))
regression_results = rbind(regression_results,c("HHS","",sprintf("%0.2e",summary(HHS_anova)$r.squared),"","",""))
write.table(regression_results,"regression_results.tsv",quote=FALSE,row.names=FALSE,sep="\t")

# write the r^2 of the multivariable regression to the console
print(paste("Multivariable regression r^2:",sprintf("%0.2e",summary(multivariate_model)$r.squared),
            "and adjusted r^2:",sprintf("%0.2e",summary(multivariate_model)$adj.r.squared)))



# calculate the timing of the peak of each wave in log COVID-19 cases
# loop over all FIPS in parallel and local maxima wiht a minimum prominance of 10^3
all_peaks_in_cases = foreach(cnt = 1:n_fips, .inorder=FALSE, .errorhandling="pass", .verbose=FALSE) %dopar% {
  r = require("zoo")
  this_fips = all_fips[cnt]
  
  temp_cases = nytimes_case_data[nytimes_case_data$fips==this_fips,]
  temp_peaks = findpeaks(temp_cases$date,temp_cases$log_cases_100k, prom=3, w=1)$x
  list(peaks = temp_peaks,
      fips = this_fips)
}

all_case_peaks =  bind_rows(lapply(all_peaks_in_cases,function(x){n=length(x$peaks);data.frame(fips=rep(x$fips,n),x["peaks"])}))
all_case_peaks = left_join(all_case_peaks,all_region,by="fips")
all_case_peaks = left_join(all_case_peaks,climate,by="fips")
all_case_peaks$day = as.numeric(strftime(all_case_peaks$peaks, format = "%j"))
all_case_peaks$date = as.Date(format(all_case_peaks$peaks, "%m-%d"),"%m-%d")

# plot histogram of peak timing
p = ggplot(all_case_peaks,aes(x=date)) +
  geom_histogram(binwidth=7) +
  ylab("# of peaks in log COVID-19 cases per week of the year") +
  theme_classic() +
  theme(axis.title.x=element_blank()) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b",)
ggsave("output/histogram_peak_timing.pdf",p,height=5,width=5)
