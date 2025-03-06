# download NOAA Monthly U.S. Climate Gridded Dataset (NClimGrid)
# https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00332

# URL
nclimgrid_link = "https://noaa-nclimgrid-daily-pds.s3.amazonaws.com/EpiNOAA/v1-0-0/csv/cty/"

# years relevant
min_year = 2020
max_year = 2023

# loop over years and download
for (y in min_year:max_year){
  for (m in 1:12){
    d = paste0(y,sprintf("%02d",m))
    
    # download and read monthly data
    temp = read.csv(paste0(nclimgrid_link,d,"-scaled.csv"),
                    colClasses=c("date"="Date"))
    
    # filter out rows without dates
    temp = filter(temp,!is.na(date))
    if (!exists("nclimgrid")){
      nclimgrid = temp
    }else{
      nclimgrid = rbind(nclimgrid,temp)
    }
  }
}

# pull year out of the date
nclimgrid$year = substring(nclimgrid$date,1,4)

# format fips
nclimgrid$fips = sprintf("%05d",nclimgrid$fips)

# calculate average monthly temperature highs, lows, and average percipitation
climate_yearly = summarize(group_by(nclimgrid,fips,year),
          max_temp = max(tmax),
          min_temp = min(tmin),
          avg_percip = mean(prcp),
          .groups="drop_last")
climate = summarize(climate_yearly,
          max_temp = mean(max_temp),
          min_temp = mean(min_temp),
          avg_percip = mean(avg_percip))
climate$diff_temp = climate$max_temp - climate$min_temp


