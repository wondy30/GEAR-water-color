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
         dwLehmann = chromaLehmann(ultraBlue, Blue, Green, Red)$wl, ## calculate dw with the mean reflectance with RGB and uBlue
         s = chromaLehmann(ultraBlue, Blue, Green, Red)$s, ## calculate purity with the mean reflectance with RGB and uBlue
         sp = chromaLehmann(ultraBlue, Blue, Green, Red)$sp,
         brightness = sqrt((Red^2 + Green^2 + Blue^2) / 3))

par(mfrow=c(1,2))
hist(cleanData$dwLehmann, breaks=20)
hist(cleanData$dwChroma, breaks=20)

summary(cleanData)
cleanData <- cleanData %>% 
  drop_na(dwLehmann)

write_csv(cleanData, "dataset/noNaDwLehmann.csv")

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
  #facet_wrap(facets = vars(Year)) +
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

#--------------------------------------##----------------------------------------#
#Calculate modal color by lake
##let us round the dw to the nearst 5

cleanData$dwLehmann <- round(cleanData$dwLehmann,0)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

cleanData$Year <- year(cleanData$Date)

modalColor <- cleanData %>% 
  group_by(Lake_name, Year) %>% 
  summarise(modalDwLehmann = getmode(dwLehmann))

#density
d <- density(modalColor$modalDwLehmann)

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
modalColor %>%
  ggplot() +
  geom_vline(data = mode_to_label, aes(xintercept = x), color = "darkgrey") + 
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", 
                                aes(x = modalDwLehmann, y = 0, fill = stat(x)), 
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

#Visual classification
#modalColor <- modalColor %>% 
#mutate(class=case_when(modalDwLehmann <= 480 ~ "B", modalDwLehmann <= 530 ~ "G", modalDwLehmann <= 580 ~ "Y", modalDwLehmann > 580 ~ "O"))


modalColor %>%
  ggplot() +
  geom_vline(data = mode_to_label, aes(xintercept = x), color = "darkgrey") + 
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", 
                                aes(x = modalDwLehmann, y = 0, fill = stat(x)), 
                                color = NA, show.legend = F, bandwidth = 3, scale = 1) +
  scale_fill_gradientn(limits = c(450, 600), colours = dw2FUI(450:600), na.value='#FFFFFF00') +
  #facet_wrap(facets = vars(class), ncol=1) +
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

#Let us classify using kmeans clustering and see if we can have groupings
library("cluster")
library("factoextra") # clustering algorithms & visualization

kmeansClass <- modalColor

res.km <- kmeans(kmeansClass[,c(-1,-2)], 4)

kmeansClass$Cluster <- res.km$cluster

windows()
kmeansClass %>%
  ggplot() +
  geom_vline(data = mode_to_label, aes(xintercept = x), color = "darkgrey") + 
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", 
                                aes(x = modalDwLehmann, y = 0, fill = stat(x)), 
                                color = NA, show.legend = F, bandwidth = 3, scale = 1) +
  scale_fill_gradientn(limits = c(450, 600), colours = dw2FUI(450:600), na.value='#FFFFFF00') +
  facet_wrap(facets = vars(Cluster), ncol=1) +
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


groups <- kmeansClass %>% 
  group_by(Cluster, Lake_name) %>% 
  summarise(n())

groups_wide <- groups %>% 
  pivot_wider(names_from = Cluster, values_from = 'n()')

plot(kmeansClass$modalDwLehmann, col = res.km$cluster)
points(res.km$cluster, col = 1:2, pch = 8, cex = 2)

#--------------------------##--------------------------------#
#modal color, in this case, for the entire time period
modalColorByLake <- cleanData %>% 
  group_by(Lake_name) %>% 
  summarise(modalDwLehmann = getmode(dwLehmann)) 
  

#recaluclate density
d <- density(cleanData$modalDwLehmann)

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

#classify
lake.km <- kmeans(modalColorByLake[,-1], 4)

modalColorByLake$Cluster <- lake.km$cluster

modalColorByLake %>%
  ggplot() +
  geom_vline(data = mode_to_label, aes(xintercept = x), color = "darkgrey") + 
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", 
                                aes(x = modalDwLehmann, y = 0, fill = stat(x)), 
                                color = NA, show.legend = F, bandwidth = 3, scale = 1) +
  scale_fill_gradientn(limits = c(450, 600), colours = dw2FUI(450:600), na.value='#FFFFFF00') +
  facet_wrap(facets = vars(Cluster), ncol=1) +
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

groups <- modalColorByLake [,-2]

test <- kmeansClass %>%  count(Cluster, Lake_name, sort = TRUE)

test2 <- test %>% 
  group_by(Lake_name, Cluster) %>% 
  summarise(max=max(n)) 


groups_wide <- groups %>% 
  pivot_wider(names_from = Cluster, values_from = `Lake_name` )


#----------------------------------##---------------------------------#
#merge the modal color with the metadata, plot on a map, classify...

library("sf")
library("rnaturalearth")      #provides a map of countries of the entire world
library("rnaturalearthdata")
require(ggspectra)

meta <- metaData[,c(6:25)]
meta <- meta %>% distinct() 
meta <- meta[-63,]  #tana kenya removed

write_csv(modalColorMeta, "dataset/modalColorMetaData.csv") #read and start from here to create a map

modalColorMeta <- left_join(modalColorByLake, meta, by="Lake_name")

world <- ne_countries(scale = "medium", returnclass = "sf")   #pull country data and choose the scale
class(world)

windows()
theme_set(theme_bw())     #classic dark-on-light theme
ggplot() +
  geom_sf(data=world) +
  #geom_text(data = world, aes(label = admin), size = 4) +
  geom_point(data=modalColorMeta, size=5, aes(x=Pour_long, y=Pour_lat, color=modalDwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "Modal color") +
  coord_sf(xlim = c(25, 45), ylim = c(-20, 20), expand = FALSE)

#-------------------------------------##---------------------------------#
#median, sd, I think the getmode function has flaw

MedianSdData <- cleanData %>% 
  group_by(Lake_name) %>% 
  summarise(medianDwLehmann=median(dwLehmann), sdDwLehmann=std(dwLehmann))
summary(MedianSdData)
#scatter plot
windows()
MedianSdData %>% 
  ggplot(aes(x = medianDwLehmann, y = sdDwLehmann)) +
  geom_point(size = 6, aes(color = medianDwLehmann)) +
  xlim(490, 590) +
  ylim(0, 40) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "Median color") +
  #geom_smooth(alpha = 0.7, method = "loess") +
  labs(
    x = "Median color (nm)",
    y = "Temporal color standard deviation (nm)",
    title = "Median color variation for each lake"
    #shape = "Reservoir"
  ) +
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=14),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust=1.5),
        legend.text=element_text(size=11),
        legend.title=element_text(size=14)) 

#on a map
meta <- read_csv("dataset/modalColorMetaData.csv")

medianColorMeta <- left_join(MedianSdData, meta, by="Lake_name")
world <- ne_countries(scale = "medium", returnclass = "sf")   #pull country data and choose the scale
class(world)

windows()
theme_set(theme_bw())     #classic dark-on-light theme
ggplot() +
  geom_sf(data=world) +
  #geom_text(data = world, aes(label = admin), size = 4) +
  geom_point(data=medianColorMeta, size=5, aes(x=Pour_long, y=Pour_lat, color=medianDwLehmann)) +
  scale_color_gradientn(limits = c(472, 588), colours = dw2FUI(472:588), na.value='#FFFFFF00', name = "Median color") +
  coord_sf(xlim = c(25, 45), ylim = c(-20, 20), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude")
