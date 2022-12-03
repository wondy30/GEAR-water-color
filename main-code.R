library(tidyverse)
library(lubridate)

#Data was downloaded in two batches, so read the batches
dataset1 <- read_csv("dataset/site1998.csv")
dataset2 <- read_csv("dataset/GEAR-points2-HLSL30-020-results.csv")

#combine
dataset <- dataset1 %>% 
  bind_rows(dataset2)

#read the metadata for the lakes, sampling location
metaData <- read_csv("dataset/locations_metadata.txt")

#Add lake name to the reluctance data
##Lakes with no names in the metadata, name was assigned with Pour_lat value
metaData <- metaData %>% 
  mutate(Lake_name = coalesce(Lake_name, as.character(Pour_lat)))

fidLakeName <- metaData %>% 
  select(FID, Lake_name) 

dataset <- left_join(dataset, fidLakeName, by=c("ID"="FID"))

#aggregate to annual average
#dataset$Date <- mdy(dataset1$Date) #Not required as the data column read in date format 

#dataset$Year <- year(dataset$Date)  #this can be done later

#Decided to calculate color DW first and then aggregate per year annual

dataset <- dataset %>% 
  filter(HLSL30_020_B01 > 0, HLSL30_020_B02 > 0, HLSL30_020_B03 > 0, HLSL30_020_B04 > 0) #drop all negative reflectance values

summary(dataset)

dim(dataset)

#Run XYZ calculation, dominant wavelength using Lehman equation








