# download mapping data for US counties and HHS regions

# a list of of the states and territories to maps (DC and Puerto Rico)
states = unique(nytimes_case_data$state)
states = states[states!='Northern Mariana Islands' &
                states!='American Samoa'&
                states!='Guam' &
                states!="Virgin Islands"]

# get us county shape files from the counties function in the tigris package
usgeo = counties(states,cb=TRUE,year=2020)

# shift alaska, hawaii, and puerto rico to fit on the plot
usgeo = shift_geometry(usgeo,position = "below",preserve_area = FALSE)

# find FIPS from state and county
usgeo = mutate(usgeo,fips = paste0(STATEFP, COUNTYFP))

# read in HHS region boundaries for plotting
HHS = read.csv("../data/HHS_regions_by_state.csv")
# change HHS regions to factor
HHS$HHS = as.factor(HHS$HHS)
# note which states and territories are in the continential US
HHS$continent = !HHS$STATE_NAME %in% c('Northern Mariana Islands','American Samoa','Guam',"Virgin Islands","Puerto Rico",'Hawaii','Alaska')

# calculate HHS regions map by combining all counties in each HHS region into one shape
usgeo_HHS = left_join(usgeo,HHS,by="STATE_NAME")
usgeo_HHS = summarize(group_by(usgeo_HHS,HHS,continent),geometry=st_union(geometry))
usgeo_HHS$center = st_centroid(usgeo_HHS$geometry)
