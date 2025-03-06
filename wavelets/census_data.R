# download census data from the 2020 Census, American Community Survey, and Gazetteer Files

# function to download census json data and convert it to a data.frame
download.census = function(link){
  json = jsonlite::fromJSON(link,flatten=TRUE) # download JSON
  df = as.data.frame(json[-1,]) # convert to a data frame
  names(df) = json[1,] # get column names from the first row the the JSON
  df$fips = paste0(df$state,df$county) # find FIPS from state and county
  return(df)
}

# api for 2020 decennial census total population (P1_001N) from the 2020 Census Demographic and Housing Characteristics File (DHC)
# https://www.census.gov/data/tables/2023/dec/2020-census-dhc.html
census_link = "https://api.census.gov/data/2020/dec/dhc?get=NAME,COUNTY,P1_001N&for=county:*"
census = download.census(census_link)
names(census)[names(census)=="P1_001N"] = "pop" # rename population estimate to pop
census$pop = as.numeric(census$pop) # convert population estimate to numeric
census$fips = paste0(census$state,census$county) # find FIPS by combining state and county

# Get county gepgraphic size estimates from Gazetteer Files
# https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.2020.html#list-tab-264479560
census_geography_file = "../data/2020_Gaz_counties_national.txt"
census_geography = read.table(census_geography_file,sep="\t",header=TRUE, quote = "\"")
census_geography$fips = sprintf("%05d",census_geography$GEOID)
census_geography$area = as.numeric(census_geography$ALAND_SQMI)

# read in % insured and income estimates from the American Community Survery Small Area Health Insurance Estimates (SAHIE) Program
# https://www.census.gov/programs-surveys/sahie.html
sahie_link = "https://api.census.gov/data/timeseries/healthins/sahie?get=NAME,PCTIC_PT,NIPR_PT,IPRCAT&for=county:*&time=2020"
sahie = download.census(sahie_link)
sahie = summarize(group_by(sahie,fips),
      perc_insured = as.numeric(PCTIC_PT[IPRCAT==0]),
      perc_poverty = as.numeric(NIPR_PT[IPRCAT==1])/as.numeric(NIPR_PT[IPRCAT==0]))

# demography estimates from the American Community Survey 5-Year Data (2009-2023) Profiles
# https://www.census.gov/data/developers/data-sets/acs-5year.2020.html#list-tab-1806015614
demography_link = "https://api.census.gov/data/2020/acs/acs5/profile?get=NAME,DP05_0001E,DP05_0024E&for=county:*"
demography = download.census(demography_link)
demography$over_65 = as.numeric(demography$DP05_0024E) / as.numeric(demography$DP05_0001E)

# combine all data into the census data frame
census = left_join(census,census_geography[c("fips","area")],by="fips")
census = left_join(census,sahie[c("fips","perc_insured","perc_poverty")],by="fips")
census = left_join(census,demography[c("fips","over_65")],by="fips")
census$density = census$pop/census$area
