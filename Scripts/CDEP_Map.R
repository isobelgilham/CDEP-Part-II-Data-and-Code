install.packages('rnaturalearth')
install.packages('rnaturalearthdata')
install.packages('rnaturalearth')
install.packages('ggrepel')
install.packages('svglite')

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggrepel)
library(svglite)

# C4, C5, A2, A4, A5, C7, C8
plots <- data.frame(longitude = c(-3.7060884,-3.7028265,-3.6857965,-3.7292048,-3.6794880,-3.6703808,-3.7013495, -3.6919503, -3.6806094, -3.6781641, -3.7770081, -3.7915448), 
                    latitude = c(57.248576,57.252960,57.225248,57.252708,57.255538,57.255040,57.257239, 57.261957, 57.168733, 57.169152, 57.148249, 57.149758 ), 
                    treatment= c('older treatment','older treatment', 'intermediate treatment', 'intermediate treatment', 
                                 'intermediate treatment', 'intermediate treatment', 'recent treatment', 'recent treatment'
                                 , 'old growth', 'mature', 'old growth', 'mature'),
                    dwc= c('dwc','dwc','dwc','dwc','dwc','dwc','dwc', 'dwc', 'control', 'control','control','control'))


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf() +
  geom_point(data = plots, aes(x = longitude, y = latitude, fill= treatment), size = 4, 
             shape = 21) +
  coord_sf(xlim = c(-3.75, -3.65), ylim = c(57.1, 57.3), expand = FALSE)


ggplot(data = world) +
  geom_sf() +
  geom_point(data = plots, aes(x = longitude, y = latitude, fill = treatment), 
             size = 4, shape = 21, color = "black") +
  coord_sf(xlim = c(-3.75, -3.65), ylim = c(57.1, 57.3), expand = FALSE) +
  theme_minimal()


install.packages("ggspatial")
install.packages("rosm")
install.packages("raster") 
install.packages("prettymapr")

library('prettymapr')
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(ggspatial)
library(rosm)

# Load in map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Establish dataframe with coordinates of plots
plots <- data.frame(
  longitude = c(-3.7060884, -3.7028265, -3.6857965, -3.7292048, -3.6794880, -3.6703808, -3.7013495, 
                -3.6919503, -3.6806094, -3.6781641, -3.7770081, -3.7915448), 
  latitude = c(57.248576, 57.252960, 57.225248, 57.252708, 57.255538, 57.255040, 57.257239, 
               57.261957, 57.168733, 57.169152, 57.148249, 57.149758), 
  treatment= c('older treatment','older treatment', 'intermediate treatment', 'intermediate treatment', 
               'intermediate treatment', 'intermediate treatment', 'recent treatment', 'recent treatment',
               'old growth', 'mature', 'old growth', 'mature')
)


                                          
##### Plot #####
library(ggmap)
library(ggplot2)
library(sf)

register_google(key = "key")  

# satellite map
satellite_map <- get_map(location = c(lon = -3.72, lat = 57.205), 
                         zoom = 12, 
                         source = "google", 
                         maptype = "satellite")

# plots to sf 
plots_sf <- st_as_sf(plots, coords = c("longitude", "latitude"), crs = 4326)

# colours
dwc_colors <- c("dwc" = "red", "control" = "cyan")

# Plot
p <- ggmap(satellite_map) +
  geom_point(data = plots, aes(x = longitude, y = latitude, fill = dwc), 
             size = 5, shape = 22, color = "black") +
  scale_fill_manual(values = dwc_colors, 
                    name = "Treatment",  # Legend title
                    labels = c("Control", "Deadwood Creation")) +   theme_minimal() +
  
  # axis label
  labs(x = "Longitude (°)", y = "Latitude (°)") +
  theme(
    axis.title = element_text(size = 16),  # Bigger axis titles
    axis.text = element_text(size = 14),  # Bigger axis numbers
    legend.text = element_text(size = 14),  # Bigger legend text
    legend.title = element_text(size = 16),  # Bigger legend title
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5)  # Boxed legend
  )

# save as svg 
ggsave("maplotnew.svg", plot = p, width = 8, height = 6, dpi = 300)


