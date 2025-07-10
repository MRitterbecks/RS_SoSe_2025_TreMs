##############################
########## Workflow 01 ##########
##############################
# download biodiversity data for fauna (species groups: birds and bats) cleaning,
# pre-processing and calculated biodiversity metrics as predictors for further modelling
##############################

setwd("/data/")
rm(list = ls())

##########  Packages ##########
#install.packages(c("readr", "dplyr", "utils", "tidyverse", "vegan"))
library(readr)
library(dplyr)
library(utils)
library(tidyverse)
library(vegan)
##############################
##############################


##########  Data Input ##########
# Reference: Heidrich, Lea (2024): Abundance data of bats, birds, spiders, and insects.
# Version 8. Biodiversity Exploratories Information System. Dataset.
# https://www.bexis.uni-jena.de/ddm/data/Showdata/31609?version=8

download.file("https://www.bexis.uni-jena.de/ddm/Data/DownloadZip/31609?format=text%2Fcsv&version=11698",
              destfile = "data/abundance_fauna.zip")
unzip("abundance_fauna.zip")

abundance_fauna <- readr::read_csv("data/abundance_fauna/31609_8_data.csv")
##############################
##############################


##########  Data Pre-Processing  ##########
# extracting data for the study site Hainich-DÜn (HEW/HAI)
abundance_fauna <- abundance_fauna %>%
  dplyr::filter(grepl("HEW", PlotID))

# filtering dataset to the species groups birds and bats
birds <- abundance_fauna %>%
  dplyr::filter(grepl("Birds", Species_group))
bats <- abundance_fauna %>%
  dplyr::filter(grepl("Bats", Species_group))

# save modified datasets
write.csv(birds, "data/abundance_fauna/birds.csv")
write.csv(bats, "data/abundance_fauna/bats.csv")
view(birds)
view(bats)

# transform data from regular format to matrix format
# Script Reference: Schäfer, Deborah; Ostrowski, Andreas; Fischer, Markus (2021):
# R Script to Transform Data from Regular Format to Matrix Format. Version 2.
# Biodiversity Exploratories Information System. Dataset.
# https://www.bexis.uni-jena.de/ddm/data/Showdata/20766?version=2


# input datasets for loop
datasets <- list(
  birds = "birds.csv",
  bats = "bats.csv"
)

# Loop workflow
for (name in names(datasets)) {
  data <- get(name)
  fileout <-  paste0(name, "_crossed.csv")
  crossed <- data %>%
    pivot_wider(names_from = Species, values_from = Abundance, values_fill = list(Abundance = 0))
  utils::write.table(
    crossed,
    file = file.path("data/abundance_fauna/", fileout),
    quote = FALSE,
    sep = ",",
    dec = ".",
    row.names = FALSE
  )
}
##############################
##############################

##########  Biodiversity-Indices ##########
# load modified datasets
birds_crossed <- readr::read_csv("data/abundance_fauna/birds_crossed.csv")
bats_crossed <- readr::read_csv("data/abundance_fauna/bats_crossed.csv")

# transforming datasets to numerical values/class
birds_crossed_num <- birds_crossed[, sapply(birds_crossed, is.numeric)]
bats_crossed_num <- bats_crossed[, sapply(bats_crossed, is.numeric)]


# define function for biodiversity indices
# calculating number of species
# calculating number of individuals
# calculating simpson diversity (alpha diversity)
# calculating (simpson)-evenness
biodiv_indices <- function(data_num) {
  species_number <- vegan::specnumber(data_num)
  data_num$Species_count <- species_number

  individuals_count <- rowSums(data_num[, -ncol(data_num)])
  data_num$Individuals_count <- individuals_count

  simpson <- vegan::diversity(data_num, "simpson")
  data_num$Simpson_diversity_a <- simpson

  evenness <- simpson/log(species_number)
  evenness[is.infinite(evenness)] <- NA
  data_num$Evenness_simpson_a <- evenness
  return(data_num)
}

# input datasets for loop
datasets <- list(
  birds_biodiv = birds_crossed_num,
  bats_biodiv = bats_crossed_num
)

# loop for calculating biodiversity indices workflow
# and adding them to the final dataset
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  dataset <- biodiv_indices(dataset)
  assign(dataset_name, dataset)
}
##############################

# add column "PlotID" (sampling plots of the study site) to dataset
birds_biodiv <- birds_biodiv %>%
  dplyr::mutate(PlotID = birds_crossed$PlotID)
bats_biodiv <- bats_biodiv %>%
  dplyr::mutate(PlotID = bats_crossed$PlotID)
birds_biodiv$Species_group <- "birds"
bats_biodiv$Species_group <- "bats"
# extract important columns (last five columns)
birds_indices <- birds_biodiv[, tail(names(birds_biodiv), 6)]
bats_indices <- bats_biodiv[, tail(names(bats_biodiv), 6)]
# rename columns
bats_indices <- bats_indices %>%
  rename(
    n_Species = Species_count,
    n_Individuals = Individuals_count
  )
birds_indices <- birds_indices %>%
  rename(
    n_Species = Species_count,
    n_Individuals = Individuals_count
  )
##############################
# save dataframes as dataset
write.csv(birds_indices, "output/predictors/birds_biodiv_indices.csv")
write.csv(bats_indices, "output/predictors/bats_biodiv_indices.csv")
birds_indices <- readr::read_csv("output/predictors/birds_biodiv_indices.csv")
bats_indices <- readr::read_csv("output/predictors/bats_biodiv_indices.csv")

# combine data to one dataset and aggregate values per sampling plot
biodiv_indices<- rbind(birds_indices, bats_indices)


biodiv_indices <- biodiv_indices %>%
  group_by(PlotID) %>%
  summarise(
    n_Species = sum(n_Species, na.rm = TRUE),
    n_Individuals = sum(n_Individuals, na.rm = TRUE),
    Simpson_diversity_a = mean(Simpson_diversity_a, na.rm = TRUE),
    Evenness_simpson_a = mean(Evenness_simpson_a, na.rm = TRUE),
    )

write.csv(biodiv_indices, "output/predictors/biodiv_indices.csv")

##############################
##### End of Workflow 01 #####
##############################









