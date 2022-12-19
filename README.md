# GEAR-water-color

Draft code for Lake Erie

setwd("D:/Research/Projects/Lakes_color/Code")

library(tidyverse)
library(lubridate)
library("sf")
library("rnaturalearth")      #provides a map of countries of the entire world
library("rnaturalearthdata")
require(ggspectra)
source("functions.r")       #from Xiao et al
require(colorscience)
require(ggridges)

lakeErieData <- read_csv("LakeErie-HLSL30-020-results.csv")

lakeErieData$month <- month(lakeErieData$Date)
lakeErieData$Year <- year(lakeErieData$Date)

lakeErieData <- lakeErieData %>% 
  filter(HLSL30_020_B01 > 0, HLSL30_020_B02 > 0, HLSL30_020_B03 > 0, HLSL30_020_B04 > 0) #drop all negative reflectance values

lakeErieData <- rename(lakeErieData, "ultraBlue" = "HLSL30_020_B01", "Blue" = "HLSL30_020_B02", "Green" = "HLSL30_020_B03", "Red" = "HLSL30_020_B04")       #rename columns to match symbol in the functions
lakeErieData <- lakeErieData %>% 
  filter(Green < 0.99, Red < 0.99, Blue < 0.99, ultraBlue < 0.99) #remove artifacts R > 1

#fliter based on water mask. see readme here: https://lpdaac.usgs.gov/documents/1326/HLS_User_Guide_V2.pdf

lakeErieData <- lakeErieData %>% 
  filter(HLSL30_020_Fmask_Water_Description=="Yes")

write_csv(lakeErieData, "watermaskfilteredErie.csv")

lakeErieData <- lakeErieData %>% 
  mutate(dwChroma = chroma2(Red, Green, Blue), ## calculate dw with the mean reflectance with RGB
         dwLehmann = chromaLehmann(ultraBlue, Blue, Green, Red)$wl, ## calculate dw with the mean reflectance with RGB and uBlue
         s = chromaLehmann(ultraBlue, Blue, Green, Red)$s, ## calculate purity with the mean reflectance with RGB and uBlue
         sp = chromaLehmann(ultraBlue, Blue, Green, Red)$sp,
         brightness = sqrt((Red^2 + Green^2 + Blue^2) / 3))

par(mfrow=c(1,2))
hist(lakeErieData$dwLehmann, breaks=20)
hist(lakeErieData$dwChroma, breaks=20)
summary(lakeErieData)

windows()
lakeErieData %>%
  ggplot() +
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", aes(x = dwLehmann, y = 0, fill = stat(x)), color = NA, show.legend = F, bandwidth = 3, scale = 1, ) + 
  scale_fill_gradientn(limits = c(450, 600), colours = dw2FUI(450:600), na.value='#FFFFFF00') +
  facet_wrap(~Category, ncol=2) + 
  scale_x_continuous(breaks = seq(470, 600, by = 30), limits = c(470, 599)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Dominant wavelength (nm)", y = "Density") +
  theme_bw()

#create a map
theme_set(theme_bw())     #classic dark-on-light theme
world <- ne_countries(scale = "medium", returnclass = "sf")   #pull country data and choose the scale
class(world)
  
#read lakes shapefile
shp <- readOGR(dsn = "D:/Research/Projects/Lakes_color/lake_erie_data/hydro_p_LakeErie", layer = "hydro_p_LakeErie")
summary(shp@data)

ggplot() +
  geom_sf(fill = "white") +
  geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "blue", 
               fill = "white") +
  #geom_text(data = world, aes(label = admin), size = 4) +
  geom_point(data=lakeErieData[lakeErieData$Year==2019,], size=6, 
             aes(x=Longitude, y=Latitude, color=dwLehmann)) +
  #ggrepel::geom_text_repel(data = lakeErieData, aes(x = Longitude, y = Latitude, label = Category, fontface = "bold", nudge_x = c(1, -1.5, 2, 2, -1), nudge_y = c(0.25, -0.25, 0.5, 0.5, -0.5)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), 
                        na.value='#FFFFFF00', name = "Dominant\nWavelength\n2019\n(nm)") +
  coord_sf(xlim = c(-83.1, -83.5), ylim = c(41.6, 41.9), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")

#Monthly average
monthlyMeanDW <- lakeErieData %>% 
  group_by(Category, Year, month) %>% 
  mutate(meanDwLehmann = mean(dwLehmann)) %>% 
  distinct(Category, Year, .keep_all=TRUE)

aug2016 <- monthlyMeanDW %>% 
  filter(month==8, Year==2016)

ggplot() +
  geom_sf(fill = "white") +
  geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "blue", 
               fill = "white") +
  geom_text(data = aug2016, aes(x=Longitude, y=Latitude, label = Category), size = 4, nudge_x = 0.02) +
  geom_point(data=aug2016, size=6, 
             aes(x=Longitude, y=Latitude, color=meanDwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), 
                        na.value='#FFFFFF00', name = "Dominant\nWavelength\nAug-2016\n(nm)") +
  coord_sf(xlim = c(-83.1, -83.5), ylim = c(41.6, 41.9), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")



monthlyMeanDW$Month <- month(monthlyMeanDW$Date, label = TRUE)

monthlyMeanDW %>% 
  filter(Category == "WE16") %>% 
  ggplot(aes(x= Month, y= meanDwLehmann, fill = meanDwLehmann)) +
  scale_fill_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), 
                        na.value='#FFFFFF00', name = "Dominant\nWavelength\nWE16\n(nm)") + 
  geom_bar(stat = "identity") +
  facet_wrap(~Year, ncol = 2) +
  coord_cartesian(ylim=c(490,590))

#Specific month over the years

monthlyMeanDW %>% 
  filter(Month == "Oct") %>% 
  ggplot(aes(x= Year, y= meanDwLehmann, fill = meanDwLehmann)) +
  scale_fill_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), 
                       na.value='#FFFFFF00', name = "October\nDominant\nWavelength\n(nm)") + 
  geom_bar(stat = "identity") +
  facet_wrap(~Category, ncol = 2) +
  coord_cartesian(ylim=c(490,590)) +
  labs(x = "Year", y = "Dominant wavelength (nm)") +
  scale_x_continuous(breaks = seq(2013, 2022, 1),labels = seq(2013, 2022, 1))

#long-term average, monthly
longtermMeanDW <- lakeErieData %>% 
  group_by(Category, month) %>% 
  mutate(longtermDwLehmann = mean(dwLehmann)) %>% 
  distinct(Category, .keep_all=TRUE)

longtermMeanDW %>% 
  ggplot(aes(x= Month, y= longtermDwLehmann, fill = longtermDwLehmann)) +
  scale_fill_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), 
                       na.value='#FFFFFF00', name = "Dominant\nWavelength\n(nm)") + 
  geom_bar(stat = "identity") +
  facet_wrap(~Category, ncol = 2) +
  coord_cartesian(ylim=c(530,580)) +
  labs(x = "Months", y = "Dominant wavelength (nm)") 
