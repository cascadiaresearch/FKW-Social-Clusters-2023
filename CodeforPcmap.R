# Code for FKW effort and sightings map for Pc manuscript 2022 November #
# written by AEH, Cascadia Research Collective #

library(sf)
library(tidyverse)
library(ggplot2)
library(ggspatial)
library(lubridate)



# start loading up effort files first #
# CRC effort #
CRCeffort <- read.csv("Effort/Effort_thruAug2022_2022NOVv2.csv")
CRCeffort$Date_ <- as.POSIXct(CRCeffort$Date_, format="%m/%d/%Y")
CRCeffort <- filter(CRCeffort, Year<2022)

# PWF effort #
PWFeffort <- read.csv("Effort/ALLPWFTracklinepoints_singlepointswaypointsremoved_clean.csv")

# PIFSC large boat effort #
PIFSC_lb_effort <- read.csv("Effort/pifsc_all_survey_effort_tracks.csv")
PIFSC_lb_effort$DateTime <- as.POSIXct(PIFSC_lb_effort$DateTime, format="%Y-%m-%d %H:%M:%S")
PIFSC_lb_effort$date <- date(PIFSC_lb_effort$DateTime)
PIFSC_lb_effort$Lon[PIFSC_lb_effort$Lon > 0] <- PIFSC_lb_effort$Lon[PIFSC_lb_effort$Lon > 0]*(-1)

# PIFSC small boat effort #
PIFSC_sb_effort <- read.csv("Effort/PIFSCsmallboattracks_HIdate_forcompiledmap_2022NOV16.csv")
PIFSC_sb_effort <- na.omit(PIFSC_sb_effort)








# build sf objects from effort files #
# CRC effort tracklines
CRCeffort_sf <- st_as_sf(CRCeffort, coords=c("Long_dd", "Lat_dd"), crs=4326)
CRCeffort_tracklines <- CRCeffort_sf %>%
  group_by(Date_, Vessel) %>%
  summarise(do_union=F) %>%
  st_cast("LINESTRING")
plot(CRCeffort_tracklines)

# PWF effort tracklines #
PWFeffort_sf <- st_as_sf(PWFeffort, coords=c("lon", "lat"), crs=4326)
PWFeffort_tracklines <- PWFeffort_sf %>%
  group_by(Year, trksegID) %>%
  summarise(do_union=F) %>%
  st_cast("LINESTRING")
plot(PWFeffort_tracklines)

# PIFSC large boat effort tracklines #
PIFSC_lb_effort_sf <- st_as_sf(PIFSC_lb_effort, coords=c("Lon", "Lat"), crs=4326)
PIFSC_lb_effort_tracklines <- PIFSC_lb_effort_sf %>%
  group_by(Cruise, date) %>%
  summarise(do_union=F) %>%
  st_cast("LINESTRING")
plot(PIFSC_lb_effort_tracklines)

# PIFSC small boat effort tracklines #
PIFSC_sb_effort_sf <- st_as_sf(PIFSC_sb_effort, coords=c("long", "lat"), crs=4326)
PIFSC_sb_effort_tracklines <- PIFSC_sb_effort_sf %>%
  group_by(date, cruise) %>%
  summarise(do_union=F) %>%
  st_cast("LINESTRING")
plot(PIFSC_sb_effort_tracklines)










# load up sightings #
sights <- read.csv("Sightings/Pc_allsightings_CRC_PWF_PIFSC_opportunistic.csv")
sights <- select(sights, Date, Group, Island, Lat, Long)
sights <- na.omit(sights)
sights_sf <- st_as_sf(sights, coords=c("Long", "Lat"), crs=4326)
plot(sights_sf)






# compile map objects #
theme_map <- function() {theme_bw() + 
    theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 0.75), axis.text = element_text(colour = 'black'), plot.title = element_text(colour = 'black', face = 'bold'))}

esri_ocean <- paste0('https://services.arcgisonline.com/arcgis/rest/services/',
                     'Ocean/World_Ocean_Base/MapServer/tile/${z}/${y}/${x}.jpeg')

# labels

label <- data.frame(label=c(paste0("Kaua\U02BBi"), paste0("Ni\U02BBihau"), paste0("O\u02BBahu"), "Maui", paste0("Kaho\u02BBolawe"), paste0("L\u0101na\u02BBi"), paste0("Moloka\U02BBi"), paste0("Hawai\U02BBi")), 
                    lat=c(22.05, 21.8, 21.5, 20.81, 20.45, 20.9, 21.3, 19.6), 
                    lon=c(-159.55, -159.9, -158.02, -156.3, -156.4, -157.3, -157.02, -155.5))
label_sf <- st_as_sf(label, coords=c("lon", "lat"), crs=4326)

# coastlines for non-esri map

coast <- st_read("Shapefiles", layer = "Coastline") # read in coastline shapefiles
coastr <- st_transform(coast, crs = 4326)# transform to projection of data 
st_crs(coastr) # check that coastline is in correct crs






# Build Map - all tracklines + sightings # 
ggplot() + 
  geom_sf(data=PWFeffort_tracklines, lwd=0.1, color="blue", alpha=0.75) +
  geom_sf(data=CRCeffort_tracklines, lwd=0.1, color="red", alpha=0.75) + 
  geom_sf(data=PIFSC_lb_effort_tracklines, lwd=0.15, color="darkgrey", alpha=0.75) +
  geom_sf(data=PIFSC_sb_effort_tracklines, lwd=0.15, color="darkgrey", alpha=0.75) + 
  geom_sf(data=coastr, fill="gray", color="black", lwd=0.25) +
  geom_sf(data=sights_sf, shape=21, fill="white", color="black", size=2) +
  coord_sf(xlim=c(-160.5, -154.75), ylim=c(18.75, 22.5), crs=4326) + 
  theme_map() + xlab("Longitude") + ylab("Latitude") +
  geom_sf_text(data=label_sf, aes(label=label, fontface='bold'), size=2.5)
ggsave("Fig1_PcMap_19992021_ALLEFFORT_tracklines_sightings_v6.tiff", dpi=450)


# map with all tracklines colored the same #

# lighest gray is just "darkgrey"#
# medium gray is "#999999"#
# darkest gray is "#848484" #
ggplot() + 
  geom_sf(data=PWFeffort_tracklines, lwd=0.1, color="#848484", alpha=0.75) +
  geom_sf(data=CRCeffort_tracklines, lwd=0.1, color="#848484", alpha=0.75) + 
  geom_sf(data=PIFSC_lb_effort_tracklines, lwd=0.15, color="#848484", alpha=0.75) +
  geom_sf(data=PIFSC_sb_effort_tracklines, lwd=0.15, color="#848484", alpha=0.75) + 
  geom_sf(data=coastr, fill="white", color="black", lwd=0.75) +
  geom_sf(data=sights_sf, shape=21, fill="white", color="black", size=2) +
  coord_sf(xlim=c(-160.5, -154.75), ylim=c(18.75, 22.5), crs=4326) + 
  theme_map() + xlab("Longitude") + ylab("Latitude") +
  annotation_north_arrow(location = "tr", which_north="true", style = north_arrow_fancy_orienteering()) + 
  annotation_scale(location="bl", text_cex = unit(1, "cm")) +
  geom_sf_text(data=label_sf, aes(label=label, fontface='bold'), size=2.5)
ggsave("PcMap_19992021_ALLEFFORT_tracklines_sightings_v8.jpg", dpi=1000)



