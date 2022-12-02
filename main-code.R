library(tidyverse)
library(lubridate)

dataset1 <- read_csv("dataset/site9991799.csv")

dataset1$Date <- mdy(dataset1$Date)

metaData <- read_csv("dataset/locations_metadata.txt")

#Lakes with no names in the metadata, name was assigned with NM and Pour_lat value

metaData <- metaData %>% 
  mutate(Lake_name = coalesce(Lake_name, as.character(Pour_lat)))

#Add lake name to the reflectance data