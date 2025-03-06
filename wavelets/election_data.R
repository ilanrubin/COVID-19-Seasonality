# download presidential election results from the MIT Election Data and Science Lab
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VOQCHQ

# location of election data file
election_file = "../data/countypres_2000-2020.csv"

# read in data
county_election = read.csv(election_file)

# format FIPS
county_election$fips = sprintf("%05d",county_election$county_fips)

# pull out 2020 presidential election results for each county and summarize by major party
election_2020 = summarize(group_by(filter(county_election,year==2020,office=="US PRESIDENT"),fips),
                          total_votes=sum(candidatevotes),
                          democrat=sum(candidatevotes[party=="DEMOCRAT"]),
                          repub=sum(candidatevotes[party=="REPUBLICAN"]))

# calculate the percentage of votes that were for a Republican Party candidate
election_2020$perc_repub = election_2020$repub/election_2020$total_votes
