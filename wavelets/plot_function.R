# function to plot case counts
plot_cases = function(df,dir,column="cases_avg_per_100k"){
  pdf(paste0(dir,paste0("/",column,".pdf")),height=7,width=7)
  par(mar=c(3,4,1,1))
  if (length(grep("log",column))>0){
    plot(df$date,df[[column]],type='l',lwd=2,col='firebrick',
      log='y',xlab="",ylab="Log of COVID-19 cases per 100k")
  }else{
    plot(df$date,df[[column]],type='l',lwd=2,col='firebrick',
      xlab="",ylab="COVID-19 cases per 100k")
  }
  dev.off()
}
