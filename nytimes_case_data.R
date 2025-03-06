# download NY Times COVID-19 case data
# https://github.com/nytimes/covid-19-data

# URL
nytimes_link = "https://raw.githubusercontent.com/nytimes/covid-19-data/refs/heads/master/rolling-averages/us-counties-"

# years relevant
min_year = 2020
max_year = 2023

# loop over years and download
for (y in min_year:max_year){
  temp = read.csv(paste0(nytimes_link,y,".csv"),
                  colClasses=c("date"="Date"))
  if (!exists("nytimes_case_data")){
    nytimes_case_data = temp
  }else{
    nytimes_case_data = rbind(nytimes_case_data,temp)
  }
}

# pull out FIPS
nytimes_case_data$fips = substr(nytimes_case_data$geoid,5,10)

# calculate log cases
nytimes_case_data$log_cases_100k = nytimes_case_data$cases_avg_per_100k
nytimes_case_data$log_cases_100k[nytimes_case_data$log_cases_100k<0] = 0
nytimes_case_data$log_cases_100k = log(nytimes_case_data$log_cases_100k)
nytimes_case_data$log_cases_100k[nytimes_case_data$log_cases_100k<0] = NA