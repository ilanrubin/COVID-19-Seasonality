# download Google Mobility Reports
# https://www.google.com/covid19/mobility/
# ZIP file: https://www.gstatic.com/covid19/mobility/Region_Mobility_Report_CSVs.zip

# location of Google Mobility data files
mobility_dir = "../data/Region_Mobility_Report_CSVs/"

# years relevant
min_year = 2020
max_year = 2022

# loop over years and read
for (y in min_year:max_year){
  temp = read.csv(paste0(mobility_dir,y,"_US_Region_Mobility_Report.csv"),
                  colClasses=c("date"="Date"))
  if (!exists("mobility")){
    mobility = temp
  }else{
    mobility = rbind(mobility,temp)
  }
}

# filter out any rows without FIPS data
mobility = filter(mobility,!is.na(census_fips_code))

# format FIPS
mobility$fips = sprintf("%05d",mobility$census_fips_code)

# which column to use for mobility periodicity ananlysis
mobility_column = c("workplaces_percent_change_from_baseline","retail_and_recreation_percent_change_from_baseline")
