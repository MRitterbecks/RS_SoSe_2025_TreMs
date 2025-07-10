##############################
########## Workflow 02 ##########
##############################
# download lidar data, cleaning, pre-processing, and calculated forest structural
# metrics as predictors for further modelling
##############################

setwd("/data/")
rm(list = ls())

##########  Packages ##########
#install.packages(c("sf", "lidR", "dplyr"))
library(sf)
library(lidR)
library(dplyr)
##############################
##############################


##########  Data Input ##########
# Reference: Ammer, Christian; Juchheim, Julia (2021):
# Three-dimensional trees, Hainich 2014, using terrestrial laser scanning, SHAPE.
# Version 2. Biodiversity Exploratories Information System. Dataset.
# https://www.bexis.uni-jena.de/ddm/data/Showdata/20044?version=2

download.file("https://www.bexis.uni-jena.de/ddm/Data/DownloadZip/20044?version=3649",
              destfile="data/three-dimensional_trees.zip")
unzip("data/three-dimensional_trees.zip")
unzip("data/three-dimensional_trees/Hainich_trees.zip")
##############################
##############################


########## Data Pre-Processing ##########
xyz_trees <- list.files("data/three-dimensional_trees/Hainich_trees/Hainich_trees",
                        pattern = "\\.xyz$", full.names = TRUE)

lidar_trees <- data.frame()
for (file in xyz_trees) {
    temp_data <- read.table(file, header = FALSE)
    plot_id <- tools::file_path_sans_ext(basename(file))
    temp_data$PlotID <- plot_id
    lidar_trees <- rbind(lidar_trees, temp_data)
}

head(lidar_trees)
colnames(lidar_trees) <- c("X", "Y", "Z", "PlotID")



# transforming xyz file into las file
lidar_trees <- lidR::LAS(lidar_trees)
# filtering outlier
hist(lidar_trees@data$Z)
Z_sd <- sd(lidar_trees@data$Z)
Z_mean <- mean(lidar_trees@data$Z)
lower_thresh <- Z_mean - 6 * Z_sd
upper_thresh <- Z_mean + 6 * Z_sd
lidar_trees <-  lidR::filter_poi(lidar_trees, Z >= lower_thresh & Z <= upper_thresh)
# filtering only vegetation
lidar_trees <- lidR::filter_poi(lidar_trees, Z >= 0 & Z <= 50)
# filtering noise
lidar_trees <- lidR::classify_noise(lidar_trees, algorithm = sor(5, 1.5))
# transforming las file into xyz
lidar_trees <- lidar_trees@data[, .(X, Y, Z, PlotID)]
hist(lidar_trees$Z)
##############################
##############################

##########  forest structural diversity metrics ##########
######################
#####################
# metrics from cloud_metrics
# metrics from Z vector
lidar_metrics <- lidar_trees %>%
  group_by(PlotID) %>%
  summarise(
    vert_sd = sd(Z, na.rm = TRUE),  # vertical standard deviation
    mean_height = mean(Z, na.rm = TRUE), # mean height
  )
print(lidar_metrics)

# aggregated values throughout mean per plot
lidar_metrics <- lidar_metrics %>%
  mutate(PlotID = sub("_[0-9]+$", "", PlotID))

lidar_metrics <- lidar_metrics %>%
  group_by(PlotID) %>%
  summarise(across(c(vert_sd, mean_height), mean, na.rm = TRUE))


write.csv(lidar_metrics, "output/predictors/lidar_metrics.csv")
##############################
##### End of Workflow 02 #####
##############################



