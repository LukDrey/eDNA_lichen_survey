# load required packages
library(here)
library(dplyr)
library(janitor)  

#####
# read in Data 
#####

# lichen mapping
lichen_forest <- read.table(here('Data', 'lichen_diversity_forests', '4460.txt'), header = T)

# read in MatrixData
lichen_matrix  <- read.table(here('Data', 'lichen_diversity_forests', 'MatrixData.txt'), header = T, sep = '|')

#transform matrix data into usable format 
trans_mat <- t(lichen_matrix)
trans_dat <- as.data.frame(trans_mat)

trans_dat <- tibble::rownames_to_column(trans_dat)

lichen_species_pres <- trans_dat %>% 
  janitor::clean_names() %>%
  janitor::row_to_names(.,1, remove_row = T)

lichen_species_pres <- lichen_species_pres %>% 
  dplyr::rename(., PlotID = Plot_ID)

# basic plot info
plot_info <- read.table(here('Data', '20826.txt'), header = T, sep = '\t')

# create big dataset containing all the info we need 

full_lichen_div <- dplyr::left_join(plot_info, lichen_species_pres, by = 'PlotID')


#####
# data clean up
#####

#only forest plots 

forest_lichen_div <- dplyr::filter(full_lichen_div,
                                        grepl('HEW|AEW|SEW', full_lichen_div$EP_PlotID))

# make a list of the column names
species_list <- gsub("\\^E.*", "", colnames(forest_lichen_div[,28:812])) %>% 
  gsub("\\^G.*", "",.) %>% 
  gsub("\\^T.*", "", .) %>% 
  gsub("\\^B.*", "", .) %>% 
  gsub("\\^M.*", "", .) %>%
  unique() 
  
species_df <- data.frame(species = species_list)

write.csv(species_df, "boch_species_list.csv")

