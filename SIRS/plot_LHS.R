# script to plot LHS output

# load packages
library("dplyr")
library("tidyr")
library("ggplot2")
library("latex2exp")

# location of directories with LHS results that will be plotted
results_dirs = c("./output/250306_232741",
                  "./output/250306_234052")

# directory to save plots
plot_dir = "./output"

# read in all output files into one dataframe
all_output = bind_rows(lapply(results_dirs,function(x) read.csv(paste0(x,"/LHS_output.csv"))))

# find the number of peaks and interepidemic period for each simulation
cols_to_save = names(all_output)[!names(all_output) %in% c("peaks","peaks_dt")]
output = all_output %>%
        summarise(mean_interepidemic = mean(peaks_dt,na.rm=TRUE),
                  median_interepidemic = median(peaks_dt,na.rm=TRUE),
                  n_peaks = n(),
                  .by=all_of(cols_to_save))

# which parameters are going to be plotted
parms_to_plot = c("imm_min",
                  "imm_shape",
                  "imm_t")
# LHS ranges for each (copied from the run_model_LHS.R file)
parms_range = t(data.frame(imm_min=c(0,0.75),
                            imm_shape=c(log(0.5,base=10),log(10,base=10)),
                            imm_t=c(25,250)))
parms_log = c(FALSE,TRUE,FALSE)

# parameter names for plotting
parameter_titles = c(TeX(r"( Long-term immunity ($\omega_{\infty}$) )"),
                TeX(r"( Immunity function shape ($k$) )"),
                TeX(r"( Immunity half-life ($\lambda$) )"))

# bin the parameters into groups to be able to calculate the proportion of simulations
# with similar parameterizations that result in annual or sub-annual waves
# the number of bins to split the ranges into
n_breaks = 25
for (i in 1:length(parms_to_plot)){
  p = parms_to_plot[i]

  temp_breaks = seq(parms_range[p,1],parms_range[p,2],length=n_breaks)
  if (parms_log[i]) temp_breaks = 10^temp_breaks
  temp_bin = cut(output[[p]], breaks=temp_breaks,labels=FALSE)

  # define the bin by its center and save that to the output data frame
  output[[paste0(p,"_bin")]] = ((temp_breaks[2:n_breaks] - temp_breaks[1:(n_breaks-1)])/2 + temp_breaks[1:(n_breaks-1)])[temp_bin]
}

# pivot output data table to have one row per parameter we will plot
output_long = pivot_longer(output,cols=paste0(parms_to_plot,"_bin"),names_to="parameter",values_to="parameter_value")

# summarize output to find the proportion of simulations with a given immunity parameter returns a periodicity category
output_summary = summarize(group_by(output_long,parameter,parameter_value,seas_amp),
      p_no_period = sum(is.na(median_interepidemic))/n(),
      p_subannual_period = sum(median_interepidemic<=0.75*365 & n_peaks>=5,na.rm=T)/n(),
      p_annual_period = sum(median_interepidemic>0.75*365 & median_interepidemic<=1.5*365 & n_peaks>=5,na.rm=T)/n(),
      p_subannual_damped = sum(median_interepidemic<=0.75*365 & n_peaks<5,na.rm=T)/n(), # damped oscillations (e.g., subannual period for some amount of time that approaches an equiblirium)
      p_annual_damped = sum(median_interepidemic>0.75*365 & median_interepidemic<=1.5*365 & n_peaks<5,na.rm=T)/n(),
      p_superannual_period = sum(median_interepidemic>1.25*365,na.rm=T)/n(),
      mean_interepidemic = mean(mean_interepidemic,na.rm=T),
      n=n())

# set the parameters to factors with LaTeX names
output_summary$parameter = factor(output_summary$parameter)
levels(output_summary$parameter) = parameter_titles

# pivot output summary data table to have one row per periodicity category we are interested in
output_summary_long = pivot_longer(output_summary,cols=c("p_annual_period","p_subannual_period"),names_to="annual_sub",values_to="proportion")

# rename periodicity levels
output_summary_long$annual_sub = factor(output_summary_long$annual_sub)
levels(output_summary_long$annual_sub) = c("'Annual periodicity'","'Subannual periodicity'")

# set colors
my_colors = c("grey30","dodgerblue3","orchid4","firebrick3","goldenrod")

# plot the proportion of simulation that result in annual or sub-annual periodicity
# color it by the amount of seasonal forcing (seas_amp)
# and show plots for each of the parameters in parms_to_plot
p = ggplot(output_summary_long,aes(x=parameter_value,
                                y=proportion,
                                # color=factor(seas_amp,labels=c("0%","2.5%","5%","10%","20%")),
                                color=factor(seas_amp),
                                group=as.factor(seas_amp))) +
  geom_path(size=1.5,alpha=0.75) + 
  scale_color_manual(values=my_colors,na.value="black")
p = p + facet_grid(vars(annual_sub),vars(parameter),scales="free_x",
                    switch = 'both', label = "label_parsed")
p = p + xlab("") +ylab("Proportion of simulations resulting in")
p = p + theme(legend.position = "top",
              legend.background = element_rect(fill = "transparent"),
              legend.title = element_blank(),
              legend.key = element_blank(),
              legend.margin=margin(t=-5,b=-10))
p = p + theme(axis.title.x=element_blank(),
              strip.background = element_rect(colour="#FFFFFF", fill="#FFFFFF"),
              strip.placement = "outside",
              panel.background = element_rect(fill = "white", colour = "black", linewidth=1.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
ggsave(paste0(plot_dir,"/proportion_annual_sunannual.pdf"),p,width=6.5,height=4)

# add logscale and vertical line at k = 1 for the classical SIRS
p = p + geom_vline(xintercept = 1,linetype="dotted",size=0.75)
p = p + scale_x_log10(breaks = c(0.5,1,5,10)) +
              annotation_logticks(sides="b")
ggsave(paste0(plot_dir,"/proportion_annual_sunannual_logscale.pdf"),p,width=6.5,height=4)


# make the same plot but with super annual and damped periodicity
# pivot output summary data table to have one row per seasonality parameter
output_summary_long = pivot_longer(output_summary,cols=c("p_superannual_period","p_annual_period","p_subannual_period","p_no_period","p_annual_damped","p_subannual_damped"),names_to="annual_sub",values_to="proportion")

# rename periodicity levels
output_summary_long$annual_sub = factor(output_summary_long$annual_sub,levels=c("p_no_period","p_superannual_period","p_annual_period","p_subannual_period","p_annual_damped","p_subannual_damped"))
levels(output_summary_long$annual_sub) = c("'stable\nequilibrium'","'super-annual\nperiodicity'","'annual\nperiodicity'","'sub-annual\nperiodicity'","'damped\nannual\nperiodicity'","'damped\nsub-annual\nperiodicity'")

# plot with all levels of periodicity
p = ggplot(output_summary_long,aes(x=parameter_value,
                                y=proportion,
                                # color=factor(seas_amp,labels=c("0%","2.5%","5%","10%","20%")),
                                group=as.factor(seas_amp))) +
  geom_path(size=1.5,alpha=0.75) + 
  scale_color_manual(values=my_colors,na.value="black")
p = p + facet_grid(vars(annual_sub),vars(parameter),scales="free",
                    switch = 'both', label = "label_parsed")
p = p + xlab("") +ylab("Proportion of simulations resulting in")
p = p + theme(legend.position = "top",
              legend.background = element_rect(fill = "transparent"),
              legend.title = element_blank(),
              legend.key = element_blank(),
              legend.margin=margin(t=-5,b=-10))
p = p + theme(axis.title.x=element_blank(),
              strip.background = element_rect(colour="#FFFFFF", fill="#FFFFFF"),
              strip.placement = "outside",
              panel.background = element_rect(fill = "white", colour = "black", linewidth=1.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

ggsave(paste0(plot_dir,"/proportion_all_periodicity.pdf"),p,width=6.5,height=6.5)

# add logscale and vertical line at k = 1 for the classical SIRS
p = p + geom_vline(xintercept = 1,linetype="dotted",size=0.75)
p = p + scale_x_log10(breaks = c(0.5,1,5,10)) +
              annotation_logticks(sides="b")
ggsave(paste0(plot_dir,"/proportion_all_periodicity_logscale.pdf"),p,width=6.5,height=6.5)
