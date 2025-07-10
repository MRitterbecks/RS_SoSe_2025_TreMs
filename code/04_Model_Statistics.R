##############################
########## Workflow 04 ##########
##############################
# download TreM data, cleaning,
# pre-processing, and modelling lm model with predictors,
# conducting statistical analysis, visualization of results
##############################
##############################

setwd("/data/")
rm(list = ls())

##########  Packages ##########
#install.packages(c("dplyr", "readr", "sf", "ggplot2", "ggpmisc", "gridExtra"))
library(dplyr)
library(readr)
library(sf)
library(ggplot2)
library(ggpmisc)
library(gridExtra)
##############################
##############################


##########  Data Input ##########
# Reference: Schall, Peter (2023): Tree related microhabitat data from 43 forest plots,
# 2017, used in "Among stand heterogeneity is key for biodiversity in managed beech
# forests but does not question the value of unmanaged forests:
# Response to Bruun & Heilmann Clausen (2021)", JAPPL. Version 6.
# Biodiversity Exploratories Information System. Dataset.
# https://www.bexis.uni-jena.de/ddm/data/Showdata/30980?version=6

# ground truth data/ field data
download.file("https://www.bexis.uni-jena.de/ddm/Data/DownloadZip/30980?format=text%2Fcsv&version=10079",
              destfile="data/TreMs.zip")
unzip("data/TreMs.zip")

TreMs <- readr::read_csv("data/TreMs/30980_6_data.csv")

########## Data Pre-Processing ##########
# caclculating mean TreMs per sample plot
TreMs <- TreMs %>%
  group_by(PlotID) %>%
  summarise(
    mean_n_TreMs = mean(Abu, na.rm = TRUE)
  )

HEW_Plots <- sf::read_sf("data/experiementalplots.gpkg")
HEW_Plots <- HEW_Plots %>%
  filter(grepl("HAI", explrtr))
HEW_Plots <- HEW_Plots %>%
  filter(grepl("forest", type))
HEW_Plots <- HEW_Plots %>%
  rename(
    PlotID = ep,
  )

TreMs_data <- TreMs
TreMs_data$PlotID <- sprintf("HEW%02d", as.numeric(sub("HEW", "", TreMs_data$PlotID)))
TreMs_data <- merge(
  TreMs_data,
  HEW_Plots[, c("PlotID", "geom")],
  by = "PlotID",
  all.x = TRUE
)

########## Predictors ##########
biodiv_indices <- readr::read_csv("output/predictors/biodiv_indices.csv")
lidar_metrics <- readr::read_csv("output/predictors/lidar_metrics.csv")
forest_structure <- readr::read_csv("output/predictors/forest_structure.csv")

biodiv_indices$PlotID <- sprintf("HEW%02d", as.numeric(sub("HEW", "", biodiv_indices$PlotID)))
lidar_metrics$PlotID <- toupper(lidar_metrics$PlotID)
lidar_metrics$PlotID <- sprintf("HEW%02d", as.numeric(sub("HEW", "", lidar_metrics$PlotID)))
forest_structure$PlotID <- sprintf("HEW%02d", as.numeric(sub("HEW", "", forest_structure$PlotID)))

# merge predictors
merge1 <- biodiv_indices %>%
  dplyr::left_join(HEW_Plots[, c("PlotID", "geom")], by = "PlotID")
merge2 <- lidar_metrics %>%
  dplyr::left_join(HEW_Plots[, c("PlotID", "geom")], by = "PlotID")
merge3 <- forest_structure %>%
  dplyr::left_join(HEW_Plots[, c("PlotID", "geom")], by = "PlotID")
merge4 <- TreMs_data %>%
  dplyr::left_join(HEW_Plots[, c("PlotID", "geom")], by = "PlotID")

# merge all values to one dataset
predictors <- merge1 %>%
  left_join(merge2, by = "PlotID") %>%
  left_join(merge3, by = "PlotID") %>%
  left_join(merge4, by = "PlotID")

# deleting columns
predictors <- subset(predictors, select = -c(...1.x, ...1.y, geom.y.x, ...1,
                                             geom.x.y, geom.y.y, geom.x.x,
                                             n_Species, n_Individuals))

write.csv(predictors, file = "output/predictors/predictors_all.csv")


########## Model ##########
# creating lm model
# simpson-diversity
model_lm_simpson <- lm(Simpson_diversity_a ~
                         vert_sd +
                         mean_dbh +
                         n_Layer +
                         mean_n_TreMs,
                       data = predictors)
summary(model_lm_simpson)

# evenness
model_lm_evenness <- lm(Evenness_simpson_a ~
                         vert_sd +
                         mean_dbh +
                         n_Layer +
                          mean_n_TreMs,
                       data = predictors)
summary(model_lm_evenness)
##########################


########## Linear Regression Analysis & Visualization ##########
model <- lm(Simpson_diversity_a ~ vert_sd, data = predictors)
summary(model)
plot1 <- ggplot(predictors, aes(x = vert_sd, y = Simpson_diversity_a)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linetype = "dashed") +
  ggpmisc::stat_poly_eq(aes(label = paste(..eq.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               color = "black") +
  ggpmisc::stat_poly_eq(aes(label = paste(..rr.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               vjust = 2,
               color = "black") +
  labs(
    x = "Vertical Standard Deviation",
    y = "Simpson-Diversity-Index"
  ) +
  theme_minimal() +
  ggtitle("p-value = ns") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

model <- lm(Simpson_diversity_a ~ mean_dbh, data = predictors)
summary(model)
plot2 <- ggplot(predictors, aes(x = mean_dbh, y = Simpson_diversity_a)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               color = "black") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               vjust = 2,
               color = "black") +
  labs(
    x = "Mean Diameter at Breast Height [m]",
    y = "Simpson-Diversity-Index"
  ) +
  theme_minimal() +
  ggtitle("p-value = ns") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

model <- lm(Simpson_diversity_a ~ n_Layer, data = predictors)
summary(model)
plot3 <- ggplot(predictors, aes(x = n_Layer, y = Simpson_diversity_a)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               color = "black") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               vjust = 2,
               color = "black") +
  labs(
    x = "Effective Number of Layers [n]",
    y = "Simpson-Diversity-Index"
  ) +
  theme_minimal() +
  ggtitle("p-value ≤ 0.1 ") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

model <- lm(Simpson_diversity_a ~ mean_n_TreMs, data = predictors)
summary(model)
plot4 <- ggplot(predictors, aes(x = mean_n_TreMs, y = Simpson_diversity_a)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", linetype = "dashed") +
  stat_poly_eq(aes(label = paste(..eq.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               color = "black") +
  stat_poly_eq(aes(label = paste(..rr.label..)),
               formula = y ~ x,
               parse = TRUE,
               size = 4,
               vjust = 2,
               color = "black") +
  labs(
    x = "Mean Number of TreMs",
    y = "Simpson-Diversity-Index"
  ) +
  theme_minimal() +
  ggtitle("p-value ≤ 0.1") +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 4)
##############################
##### End of Workflow 04 #####
##############################
