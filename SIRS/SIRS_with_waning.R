### immune waning function
waning = function(x,imm_min,imm_shape,imm_t) {
  return((1-imm_min) * (exp(-(x/(imm_t) * ((-log(0.5))^(1/imm_shape)))^imm_shape)) + imm_min)
}

### seasonality function
seasonal = function(x,amp,period,shift) {
  seasonal = 1 + amp * cos((x*2/period + shift) * pi)
  return(seasonal)
}

### function with IBM implementation of SIRS model with waning immunity
SIRS_with_waning = function(parms,
                            save_plots=TRUE,
                            save_output_file=TRUE,
                            output_dir="."){
  # start timer
  start_time <- Sys.time()

  # turn the parameters list into variables
  for (i in 1:length(parms)){
    assign(names(parms[i]),parms[[i]])
  }

  # set seed
  if (seed>0) set.seed(seed)

  # set time vector that the simulation will run over
  time <- seq(1,t_max,by=dt)
  time_steps = length(time)

  ### seasonality vector
  seas = seasonal(time,seas_amp,seas_period,seas_shift)

  # create population arrays
  infected = rep(FALSE,N) # infection status
  infection_end_time = rep(NA,N)
  infection_start_time = rep(NA,N)

  # create array to store history of the number of infected individuals
  I_hist = rep(0,time_steps)

  # initialize initially infected individuals
  I = sample(N, I0) # randomly choose initial infected
  I_hist[1] = I0
  num_I = I0 # counter of current number of infected individuals

  # draw time infected from exponential distribution for each infected individual
  for (i in I){
    temp_time = round(rexp(1,gamma*dt))
    if (temp_time > 0){
      infected[i] = TRUE
    }
    infection_start_time[i] = 1
    infection_end_time[i] = temp_time + 1
  }

  ### run simulation
  for(t in 2:time_steps){
    # calculate binomially distributed number of exposed based on number infected and beta
    # 'exposed' = will become infected if completely susceptible
    n_exposed = rbinom(n = 1, size = N, prob = 1 - exp(-(num_I+0)*R0*gamma*dt*seas[t]/N)) #with continual re-seeding (+1)
    if (n_exposed == 0 & num_I == 0) n_exposed = 1 # ensure that the epidemic does that end
    which_exposed = sample(N, size = n_exposed, replace = F)

    # filter exposed for only those that are not already infected
    which_suscept_exposed = which_exposed[!infected[which_exposed]] 
    n_suscept_exposed = length(which_suscept_exposed)

    # calculate who is infected based on who is exposed and immunity level
    if (n_suscept_exposed > 0){
      # set immunity for all exposed to 0, then we will recalulate immunity based on time since last infectious
      immunity = rep(0,n_suscept_exposed)

      # exposed individuals with immunity
      which_previously_infected = which_suscept_exposed[!is.na(infection_end_time[which_suscept_exposed])]
      immunity[!is.na(infection_end_time[which_suscept_exposed])] = waning(t-infection_end_time[which_previously_infected],imm_min,imm_shape,imm_t)

      # binomially distributed infection success based on amount of individual immunity
      infection_success = rbinom(n_suscept_exposed, size = 1, prob = 1-immunity)
      which_infected = which_suscept_exposed[infection_success == 1]

      # exponentially distributed infectious times
      if (length(which_infected)>0){
        infection_times = round(rexp(length(which_infected),gamma*dt))
        infection_start_time[which_infected] = t
        infection_end_time[which_infected] = t + infection_times
        infected[which_infected] = TRUE
      }
    }

    # individuals recover based on pre-calculated infection end tim
    infected[which(infection_end_time==t)] = FALSE
    num_I = sum(infected)
    I_hist[t] = num_I
  }

  # calculate peaks in the number of infections
  peaks = findpeaks(time,I_hist/N,prom=0.05,span=0.05,w=1)

  if (save_output_file){
    ### file to save output
    dummy=file.create(paste0(output_dir,"/simulation_output.dat"),showWarnings=FALSE)
    out_file = file(paste0(output_dir,"/simulation_output.dat"))
    writeLines(paste(time,(N-I_hist)/N,I_hist/N,sep="\t\t"),out_file)
    close(out_file)

    # file with the peaks in cases
    write.table(peaks$i,file=paste0(output_dir,"/peaks.dat"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="")
  }

  if (save_plots){
    ### save individual immunity curve
    pdf(paste0(output_dir,"/immunity_curve.pdf"),height=2,width=3)
    par(mar=c(4,4,1,1))
    plot(time,waning(time,imm_min,imm_shape,imm_t),
          type = "l", lwd = 4,
          ylim = c(0, 1),
          xlim = c(0,365),
          cex=0.5,
          yaxt='n',
          ylab="Immunity",xlab="Time (days)")
    lines(c(imm_t,imm_t),c(0,1),col=rgb(0,0,0),lwd=2,lty='dotted')
    axis(2,at=c(0,0.5,1),las=2)
    dummy=dev.off()

    # plot with proportion of population infected, amount of population immunity,
    # and the amount of seasonal forcing
    pdf(paste0(output_dir,"/infections_and_forcing.pdf"),width=5,height=5)
    par(mar=c(4,4,1,2))
    plot(time,I_hist/N,
          type='l',lwd=3,col="firebrick",
          xaxt='n',
          ylim=c(0,ceiling(max(I_hist/N)*10)/10),
          ylab="Proportion of poopulation infected (red)\n and seasonal forcing (grey)", xlab="Time (days)")
    axis(1, at=seq(0,t_max,by=365),labels=round(seq(0,t_max,by=365)/365),las=2)

    par(new = TRUE)
    plot(time, seas, type = "l", lwd=2, axes = FALSE, bty = "n", ylim=c(0,2),xlab = "", ylab = "", col='#00000060')
    axis(side=4)
    dummy=dev.off()

    # plot just the proportion of the population infected
    # with circles showing the peaks in the infection curve
    pdf(paste0(output_dir,"/infections.pdf"),width=5,height=5)
    par(mar=c(4,4,1,2))
    plot(time,I_hist/N,
          type='l',lwd=3,
          xaxt="n",
          ylab="Proportion of poopulation infected", xlab="Time (days)")
    axis(1, at=seq(0,t_max,by=365),labels=round(seq(0,t_max,by=365)/365),las=2)
    points(time[peaks$i],(I_hist/N)[peaks$i],col='red')
    dummy=dev.off()
  }
  
  end_time <- Sys.time()
  print(paste0(output_dir,": ",end_time-start_time," to run sim"))

  # return output
  out = list(IBM = data.frame(time=time,I=I_hist/N),
            peaks=peaks$i)
  return(out)
}
