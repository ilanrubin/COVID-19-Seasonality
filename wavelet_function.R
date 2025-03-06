# function to run wavelet analysis and save figures and the peaks of the period of the wave decomposition
run_wavelet = function(data, incidence_col,
                        time_col="date", time_format="date",
                        identifier="",
                        save_plots=TRUE, write_files=TRUE,
                        invisible_TF = TRUE,
                        directory=".",...){
  # add underbar to identifier if necessary
  if (identifier != ""){
    identifier = paste0("_",identifier)
  }

  # ensure date column exists and is in the right format
  if (time_col != "date"){
    data$date = data[[time_col]]
  }
  if (time_format == "date"){
    data$date = as.Date(data$date)
    show_date = TRUE
  }else{
    data$date = as.numeric(data$date)
    show_date = FALSE
  }

  # filter out any times that don't exist
  data = data[!is.na(data[[incidence_col]]),]

  # calculate dt
  dt = as.numeric(data$date[2] - data$date[1])

  # run analysis
  if (invisible_TF){
    wavelet = invisible(analyze.wavelet(data, incidence_col,
                             loess.span = 0,
                             dt = dt, dj = 1/250,
                             lowerPeriod = 64,
                             upperPeriod = 512,
                             make.pval = TRUE, n.sim = 10,
                             verbose=FALSE))
  }else{
    wavelet = analyze.wavelet(data, incidence_col,
                             loess.span = 0,
                             dt = dt, dj = 1/250,
                             lowerPeriod = 64,
                             upperPeriod = 512,
                             make.pval = TRUE, n.sim = 10)
  }
  


    # find peaks in period reconstruction
    # peaks = argmax(wavelet$Period, wavelet$Power.avg, w=1, span=0.01)
    peaks = findpeaks(wavelet$Period, wavelet$Power.avg, prom=0.01, w=1)

  # plot and save figures
  if (save_plots){

    png(paste0(directory,"/wavelet_power_spectrum_period",identifier,".png"),width=7,height=7,units="in",res=600)
    par(mar=c(3,4,1,1))
    wavelet_image = wavelet$Power
    wavelet_image = log(wavelet_image,base=2)
    wavelet_image[wavelet_image<log(0.1/(2^4),base=2)] = log(0.1/(2^4),base=2)
    image(wavelet$series$date,wavelet$Period,t(wavelet_image),
          xlab="",ylab="Wave Period (days)",
          # log='y',
          col=hcl.colors(1000,palette='Geyser'))
          # col=hcl.colors(1000,palette='Geyser')[round(log(1:250)/log(250)*1000)])
    ridge = which(wavelet$Ridge==1,arr.ind=TRUE)
    points(wavelet$series$date[ridge[,2]],wavelet$Period[ridge[,1]],pch='.',cex=2,col="#000000")
    lines(c(min(wavelet$series$date),max(wavelet$series$date)),c(365,365),lty="dashed")
    dev.off()

    pdf(paste0(directory,"/data_with_reconstructed_wave",identifier,".pdf"),width=5,height=5)
    par(mar=c(4,3,1,1))
    reconstruct(wavelet, plot.waves = FALSE, only.ridge = F,lwd = c(2,4),
              show.date=show_date,
              col=c("#000000","#b22222b0"),
              show.legend=FALSE,
              verbose=FALSE)
    dummy=dev.off()

    pdf(paste0(directory,"/wavelet_period",identifier,".pdf"),width=7,height=7)
    par(mar=c(4,4,1,1))
    y_max = max(wavelet$Power.avg)
    y_min = min(wavelet$Power.avg)
    plot(wavelet$Period, wavelet$Power.avg,type='l',lwd=3,xlab="Wavelet Period",ylab="Wavelet Power")
    lines(c(365,365),c(y_min,y_max),lty='dashed',lwd=2,col="#00000099")
    if (length(peaks$i)>0){
      for (i in 1:length(peaks$i)){
        # lines(c(peaks$x[i],peaks$x[i]),c(y_min,y_max),lwd=2,col="#00000099")
        if (wavelet$Power.avg.pval[peaks$i[i]] < 0.05){
          color = "#b22222aa"
        }else{
          color = "#000000aa"
        }
        lines(c(peaks$x[i],peaks$x[i]),c(y_min,wavelet$Power.avg[peaks$i[i]]),lwd=2,col=color)
        points(peaks$x[i],wavelet$Power.avg[peaks$i[i]],cex=2,pch=16,col=color)
      }
      text(peaks$x+6,y_min,labels=round(peaks$x,2),srt=-90,adj=1)
    }
    text(365+6,y_min,labels=365,srt=-90,adj=1)
    dummy=dev.off()

  }

  # save peaks and full period in file
  if (write_files){
    write.table(data.frame(time=wavelet$Period,power=wavelet$Power.avg,p.val=wavelet$Power.avg.pval),
                file=paste0(directory,"/wavelet_period",identifier,".dat"),
                quote=FALSE,row.names=FALSE)
    write.table(data.frame(time=peaks$x,power=wavelet$Power.avg[peaks$i],p.val=wavelet$Power.avg.pval[peaks$i]),
                file=paste0(directory,"/wavelet_peaks",identifier,".dat"),
                quote=FALSE,row.names=FALSE)
  }
  
  # print(paste(directory,"Finished with wavelet"))
  return(wavelet)
}

