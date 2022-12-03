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

write_csv(dataset, "dataset/cleanData.csv")

cleanData <- read_csv("dataset/cleanData.csv")          #row data

#Run XYZ calculation, dominant wavelength using Lehman equation

source("functions.r")       #from Xiao et al
require(sf)
require(tidyverse)
require(ggspectra)
require(colorscience)
require(ggridges)

colnames(cleanData)

cleanData <- rename(cleanData, "ultraBlue" = "HLSL30_020_B01", "Blue" = "HLSL30_020_B02", "Green" = "HLSL30_020_B03", "Red" = "HLSL30_020_B04")       #rename columns to match symbol in the functions

windows()
hist(cleanData$Green)
cleanData <- cleanData %>% 
  filter(Green < 0.6)       #Clean cloud signals, reflectance higher than 0.6

hist(cleanData$Green)

cleanData <- cleanData %>% 
  mutate(dwPost = chroma(Red, Green, Blue))   #first method using chroma

cleanData <- cleanData %>% 
  mutate(dwChroma = chroma2(Red, Green, Blue), ## calculate dw with the mean reflectance with RGB
         dwLehmann = chromaLehmann(ultraBlue, Blue, Green, Red)$wl, ## calculate dw with the mean reflectances with RGB and uBlue
         s = chromaLehmann(ultraBlue, Blue, Green, Red)$s, ## calculate purity with the mean reflectance with RGB and uBlue
         sp = chromaLehmann(ultraBlue, Blue, Green, Red)$sp,
         brightness = sqrt((Red^2 + Green^2 + Blue^2) / 3))

par(mfrow=c(1,2))
hist(cleanData$dwLehmann, breaks=20)
hist(cleanData$dwChroma, breaks=20)

summary(cleanData)
cleanData <- cleanData %>% 
  drop_na(dwLehmann)


## density plot
### label values
d <- density(cleanData$dwLehmann)

max_mode <- tibble(x = d$x, y = d$y) %>% 
  mutate(grp = x < 529) %>% 
  group_by(grp) %>% 
  summarise(ymax = max(y),
            x = x[which(y == ymax)]) %>% 
  ungroup() %>% 
  rename(y = ymax) %>% 
  select(-grp)

min_mode <- tibble(x = d$x, y = d$y) %>% 
  filter(x >= 495 & x <= 574) %>% 
  slice_min(order_by = y)

mode_to_label <- bind_rows(max_mode, min_mode)

cleanData$Year <- year(cleanData$Date)

cleanData %>%
  ggplot() +
  geom_vline(data = mode_to_label, aes(xintercept = x), color = "darkgrey") + 
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", 
                                aes(x = dwLehmann, y = 0, fill = stat(x)), 
                                color = NA, show.legend = F, bandwidth = 3, scale = 1) +
  scale_fill_gradientn(limits = c(450, 600), colours = dw2FUI(450:600), na.value='#FFFFFF00') +
  facet_wrap(facets = vars(Year)) +
  scale_x_continuous(breaks = seq(470, 600, by = 30), limits = c(470, 600), 
                     sec.axis = dup_axis(name = NULL, breaks = mode_to_label$x, 
                                         labels = paste(format(mode_to_label$x, digits = 2)))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), 
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(plot.margin = unit(x = c(0, 0.05, 0, 0.05), units = "in"),
    text = element_text(size = 7),
    axis.text.x = element_text(angle = 0),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank())




