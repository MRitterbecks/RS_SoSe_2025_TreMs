##############################
########## Workflow 03 ##########
##############################
# download forest structure data (DBH, effective number of layers), cleaning,
# pre-processing, and calculated forest structural
# metrics as predictors for further modelling
##############################
##############################

setwd("/data/")
rm(list = ls())

##########  Packages ##########
#install.packages("dplyr")
library(dplyr)
##############################
##############################


##########  Data Input ##########
# Reference: Schall, Peter; Ammer, Christian; Schulze, Ernst-Detlef (2017):
# 1st forest inventory on all forest EPs, single tree data, 2008 - 2014. Version 2.
# Biodiversity Exploratories Information System. Dataset.
# https://www.bexis.uni-jena.de/ddm/data/Showdata/18268?version=2

# Reference: Ehbrecht, Martin; Weißing, Kim (2024):
# Effective Number of Layers (ENL) - Forest plots - 2014, 2019, 2023. Version 5.
# Biodiversity Exploratories Information System. Dataset.
# https://www.bexis.uni-jena.de/ddm/data/Showdata/31653?version=5

download.file("https://www.bexis.uni-jena.de/ddm/Data/DownloadZip/18268?format=text%2Fcsv&version=3459",
              destfile="data/Forest_Inventory.zip")
unzip("data/forest_inventory.zip")

download.file("https://www.bexis.uni-jena.de/ddm/Data/DownloadZip/19986?format=text%2Fcsv&version=3626",
              destfile="data/Number_Effective_Layers.zip")
unzip("data/number_effective_layers.zip")

DBH <- readr::read_csv("data/forest_inventory/18268_2_data.csv")
n_Layer <- readr::read_csv("data/number_effective_layers/31653_5_data.csv")



########## Data Pre-Processing ##########
# filtering to sample plots
DBH <- DBH %>%
  filter(grepl("HEW", EP_Plotid))
n_Layer <- n_Layer %>%
  filter(grepl("HEW", EP))

# extracting values
n_Layer <- n_Layer[, -c(3, 4, 5, 6, 7, 8)]

# calculating mean dbh values
DBH_values <- DBH %>%
  group_by(EP_Plotid) %>%
  summarise(
    mean_dbh = mean(dbh, na.rm = TRUE)
    )

# rename columns
DBH_values <- DBH_values %>%
  rename(
    PlotID = EP_Plotid,
    )
n_Layer <- n_Layer %>%
  rename(
    PlotID = EP,
  )

# convert cm into m
DBH_values$mean_dbh <- DBH_values$mean_dbh / 100


# merge datasets to one
forest_structure <- merge(
  DBH_values,
  n_Layer[, c("PlotID", "enl_2014")],
  by = "PlotID",
  all.x = TRUE
)

forest_structure <- forest_structure %>%
  rename(
    n_Layer = enl_2014
  )

write.csv(forest_structure, "output/predictors/forest_structure.csv")
##############################
##### End of Workflow 03 #####
##############################
