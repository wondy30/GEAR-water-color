library(tidyverse)
library(lubridate)
library("sf")
library("rnaturalearth")      #provides a map of countries of the entire world
library("rnaturalearthdata")
require(ggspectra)
source("functions.r")       #from Xiao et al
require(colorscience)
require(ggridges)

#This code is used to test the best way to summarize water color using Lake Victoria as an example

dwDataClean <- read_csv("dataset/noNaDwLehmann.csv")        #clean dw data

victoria <- dwDataClean %>% 
  filter(Lake_name=="Victoria")
victoria$Year <- year(victoria$Date)

world <- ne_countries(scale = "medium", returnclass = "sf")   #pull country data and choose the scale
class(world)
colnames(victoria)
windows()
theme_set(theme_bw())     #classic dark-on-light theme


vic18 <- victoria %>% 
  filter (Year==2022) 
ggplot() +
  geom_sf() +
  #geom_text(data = world, aes(label = admin), size = 4) +
  geom_point(data=vic18, size=5, aes(x=Longitude, y=Latitude, color=dwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "Water Color\n2022 (nm)") +
  coord_sf(xlim = c(31.5, 34.6), ylim = c(0, -3), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")
hist(victoria$dwLehmann)



#Round to the nearest 2
library(plyr)
victoria <- victoria %>% 
  mutate(dwLehmannR2= round_any(dwLehmann, 2, f=round))
#calculate modal
library(modeest)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

modal_victoria <- victoria %>% 
  group_by(Year) %>% 
  summarise(modalColor=getmode(dwLehmannR2), medianColor=median(dwLehmann), meanColor=mean(dwLehmann)) 

sumvictoria <- modal_victoria %>% pivot_longer(-Year, names_to = "summary", values_to = "dwLehmann")
  
# Scatter plot

ggplot() +
  geom_point(data=sumvictoria, size=6, aes(x=Year, y=dwLehmann, color=dwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "Water color\n(nm)") +
  scale_x_continuous(breaks = seq(2013, 2022, by = 2)) +
  facet_wrap( ~summary)
