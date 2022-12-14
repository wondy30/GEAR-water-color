library(tidyverse)
library(lubridate)
library("sf")
library("rnaturalearth")      #provides a map of countries of the entire world
library("rnaturalearthdata")
require(ggspectra)
source("functions.r")       #from Xiao et al
require(colorscience)
require(ggridges)

setwd("/water_color/GEAR-water-color")

lakes3Data <- read_csv("dataset/threelakesData.csv")
#format and calculate month, year, will be used for aggregation
lakes3Data$Date <- mdy(lakes3Data$Date)
lakes3Data$month <- month(lakes3Data$Date)
lakes3Data$Year <- year(lakes3Data$Date)

lakes3Data <- lakes3Data %>% 
  filter(HLSL30_020_B01 > 0, HLSL30_020_B02 > 0, HLSL30_020_B03 > 0, HLSL30_020_B04 > 0) #drop all negative reflectance values

lakes3Data <- rename(lakes3Data, "ultraBlue" = "HLSL30_020_B01", "Blue" = "HLSL30_020_B02", "Green" = "HLSL30_020_B03", "Red" = "HLSL30_020_B04")       #rename columns to match symbol in the functions
lakes3Data <- lakes3Data %>% 
  filter(Green < 0.99, Red < 0.99, Blue < 0.99, ultraBlue < 0.99) #remove artifacts R > 1
  
#fliter based on water mask. see readme here: https://lpdaac.usgs.gov/documents/1326/HLS_User_Guide_V2.pdf

lakes3Data <- lakes3Data %>% 
  filter(HLSL30_020_Fmask_Water_Description=="Yes")

write_csv(lakes3Data, "dataset/watermaskfiltered.csv")

#read watermaskfiltered.csv from here to start using the clean data
lakes3Data <- read_csv("dataset/watermaskfiltered.csv")

lakes3Data <- lakes3Data %>% 
  mutate(dwChroma = chroma2(Red, Green, Blue), ## calculate dw with the mean reflectance with RGB
         dwLehmann = chromaLehmann(ultraBlue, Blue, Green, Red)$wl, ## calculate dw with the mean reflectance with RGB and uBlue
         s = chromaLehmann(ultraBlue, Blue, Green, Red)$s, ## calculate purity with the mean reflectance with RGB and uBlue
         sp = chromaLehmann(ultraBlue, Blue, Green, Red)$sp,
         brightness = sqrt((Red^2 + Green^2 + Blue^2) / 3))

par(mfrow=c(1,2))
hist(lakes3Data$dwLehmann, breaks=20)
hist(lakes3Data$dwChroma, breaks=20)

summary(lakes3Data)

#Calculate annual mean dw
lakes3Data <- lakes3Data %>% 
  group_by(ID, Year) %>% 
  mutate(meanDwLehmann = mean(dwLehmann)) %>% 
  distinct(ID, Year, .keep_all=TRUE)

#density plot
d <- density(lakes3Data$meanDwLehmann)

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

windows()
lakes3Data %>%
  ggplot() +
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", aes(x = meanDwLehmann, y = 0, fill = stat(x)), color = NA, show.legend = F, bandwidth = 3, scale = 1, ) + 
  scale_fill_gradientn(limits = c(450, 600), colours = dw2FUI(450:600), na.value='#FFFFFF00') +
  facet_wrap(~Category, ncol=1) + 
  scale_x_continuous(breaks = seq(470, 600, by = 30), limits = c(470, 599)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Dominant wavelength (nm)", y = "Density") +
  theme_bw()

#create map
library("sf")
library("rnaturalearth")      #provides a map of countries of the entire world
library("rnaturalearthdata")
require(ggspectra)

unique(lakes3Data$Category)
theme_set(theme_bw())     #classic dark-on-light theme
world <- ne_countries(scale = "medium", returnclass = "sf")   #pull country data and choose the scale
class(world)

victoria <- lakes3Data %>% 
  filter (Category=="Victoria")  
ziway <- lakes3Data %>% 
  filter (Category=="Ziway")  
tanga <- lakes3Data %>% 
  filter (Category=="Tanganyika")  
ggplot() +
  geom_sf(data=world, fill = "white") +
  #geom_text(data = world, aes(label = admin), size = 4) +
  geom_point(data=victoria, size=5, aes(x=Longitude, y=Latitude, color=dwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "DW (nm)") +
  coord_sf(xlim = c(31, 35), ylim = c(-3.0, 0.5), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")

ggplot() +
  geom_sf(data=world, fill = "white") +
  #geom_text(data = world, aes(label = admin), size = 4) +
  geom_point(data=ziway, size=5, aes(x=Longitude, y=Latitude, color=dwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "DW (nm)") +
  coord_sf(xlim = c(38.7, 38.95), ylim = c(7.85, 8.12), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")

ggplot() +
  geom_sf(data=world, fill = "white") +
  #geom_text(data = world, aes(label = admin), size = 4) +
  geom_point(data=tanga, size=4, aes(x=Longitude, y=Latitude, color=dwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "DW (nm)") +
  coord_sf(xlim = c(29.0, 31.3), ylim = c(-8.9, -3.2), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")
