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
plot_info <- read.table(here('Data', 'plot_info', '20826.txt'), header = T, sep = '\t')

# forest management intensity index ForMI

formi <- read.table(here('Data', 'ForMI', '16466.txt'), header = T, sep = '\t')

formi <- formi %>% rename(., EP_PlotID = EP_Plotid)

# create big dataset containing all the info we need 

plots_landuse <- dplyr::left_join(plot_info, formi, by = 'EP_PlotID')

full_lichen_div <- dplyr::left_join(plots_landuse, lichen_species_pres, by = 'PlotID')


#####
# data clean up
#####

#select all columns of epiphytic species growing on Fagus sylvatica

clean_full_lichen_div <- dplyr::select(full_lichen_div, EP_PlotID, PlotID,
                                Exploratory, ActivePlot, Comments,
                                VIP, MIP, RW, HW, Elevation, Slope, 
                                Aspect, MainTreeSpecies, ForMI, Iharv,
                                Inonat, Idwcut, contains("^E"))

#only forest plots 

clean_full_lichen_div <- dplyr::filter(clean_full_lichen_div,
                                        grepl('HEW|AEW|SEW', clean_full_lichen_div$EP_PlotID))

# calculate which species occurs on how many plots 
# first make the pres abs data numeric 

for (i in 18:613) {
  clean_full_lichen_div[,i] <- as.numeric(clean_full_lichen_div[,i])
  
}

# calculate column sum

clean_sum <- as.data.frame(colSums(clean_full_lichen_div[,18:613], na.rm = T))

# rename column 
clean_sum <- clean_sum %>%
  rename(., num_plot = `colSums(clean_full_lichen_div[, 18:613], na.rm = T)`)

#filter anything that is 0 
clean_sum <- dplyr::filter(clean_sum, num_plot > 0 )

clean_sum <- tibble::rownames_to_column(clean_sum)

#order according to abundance
clean_sum <- clean_sum[order(clean_sum$num_plot, decreasing = T ), ]
clean_sum$rowname

alb_sum <- 

#look for promising species
dplyr::filter(alb_sum, grepl("Evernia", rowname))
dplyr::filter(clean_sum, grepl("Pseudevernia", rowname))
dplyr::filter(alb_sum, grepl("Melanelia", rowname))
dplyr::filter(clean_sum, grepl("Hypogymnia", rowname))

col_names_full <- c('EP_PlotID', 
                   'ActivePlot', 'ForMI') 


# on which plots does Evernia prunastri occur

evernia_abund <- select(clean_full_lichen_div, all_of(col_names_full),
                           contains("Evernia") )

which(evernia_abund$`Evernia prunastri^E(Fag_syl)` > 0) #row 146

evernia_abund[146,]

# Lecanora

leconora_abund <- select(clean_full_lichen_div, all_of(col_names_full),
                        contains("Lecanora") )

rowSums(leconora_abund[4:73])
colSums(leconora_abund[4:73], na.rm = T)



# on which plots does Pseudevernia furfuracea occur 

pseudevernia_abund <- select(clean_full_lichen_div, all_of(col_names_full),
                             contains("Pseudevernia") )

which(pseudevernia_abund$`Pseudevernia furfuracea^E(Pic_abi)` > 0)

pseudevernia_abund[4,]

# on which plots does Melanelia glabratula occur 

melanelia_abund <- select(clean_full_lichen_div, all_of(col_names_full),
                          contains("Melanelia") )

which(melanelia_abund$`Melanelia glabratula^E(Fag_syl)` > 0)

melanelia_abund[c(1,10,12,19,118,120,126,128,130,135,139,140,146,149,150),1]

which(melanelia_abund$`Melanelia subaurifera^E(Fag_syl)` > 0)

melanelia_abund[146, 1]

# on which plots does hypogymnia occur 

hypogymnia_abund <- select(clean_full_lichen_div, all_of(col_names_full), 
                           contains("Hypogymnia"))

which(hypogymnia_abund$`Hypogymnia physodes^E(Fag_syl)` > 0)

hypogymnia_abund[c(118,128,146), 1]

# on which plots does parmelia sulcata occur 

parmelia_abund <- select(clean_full_lichen_div, all_of(col_names_full), 
                           contains("Parmelia"))

which(parmelia_abund$`Parmelia sulcata^E(Fag_syl)` > 0)

parmelia_abund[c(1, 10,12,118,119,120,128,130,133,139,140,146,149,150), 1]

# on which plots does Phaephyscia orbicularis occur

phaephyscia_abund <- select(clean_full_lichen_div, all_of(col_names_full), 
                         contains("Phaeophysica"))

which(phaephyscia_abund$`Phaeophysica orbicularis^E(Fag_syl)` > 0)

parmelia_abund[c(8, 9,10,118,119,122,123,129,139,140,148,149,151), 1]

# all species occuring on fagus 

clean_full_lichen_div_fag <- dplyr::select(clean_full_lichen_div, EP_PlotID, PlotID,
                                        Exploratory, ActivePlot, Comments,
                                        VIP, MIP, RW, HW, Elevation, Slope, 
                                        Aspect, MainTreeSpecies, ForMI, Iharv,
                                        Inonat, Idwcut, contains("Fag_syl"))



# calculate which species occurs on how many plots 
# first make the pres abs data numeric 

for (i in 18:112) {
  clean_full_lichen_div_fag[,i] <- as.numeric(clean_full_lichen_div_fag[,i])
  
}

# calculate column sum

clean_sum_fag <- as.data.frame(colSums(clean_full_lichen_div_fag[,18:112], na.rm = T))

# rename column 
clean_sum_fag <- clean_sum_fag %>%
  rename(., num_plot = `colSums(clean_full_lichen_div_fag[, 18:112], na.rm = T)`)

#filter anything that is 0 
clean_sum_fag <- dplyr::filter(clean_sum_fag, num_plot > 0 )

clean_sum_fag <- tibble::rownames_to_column(clean_sum_fag)

#order according to abundance
clean_sum_fag <- clean_sum_fag[order(clean_sum_fag$num_plot, decreasing = T ), ]
clean_sum_fag$rowname

#####
# Hainich
#####

# subset to Hainich forests Ep's only

hainich_lichen_div <- dplyr::filter(full_lichen_div, grepl('HEW', full_lichen_div$EP_PlotID) )

# only select columns that contain data on epiphytic lichen but keep additional information

tokeep_hai <- colnames(hainich_lichen_div[,c(1:3,5,28)])

hainich_epi_lichen_div <- dplyr::select(hainich_lichen_div, all_of(tokeep_hai),  contains('^E'))

# calculate which species occurs on how many plots 
# first make the pres abs data numeric 

for (i in 6:601) {
  hainich_epi_lichen_div[,i] <- as.numeric(hainich_epi_lichen_div[,i])
  
}

# calculate column sum
hai_sum <- as.data.frame(colSums(hainich_epi_lichen_div[,6:601], na.rm = T))

# rename column 
hai_sum <- hai_sum %>%
  rename(., num_plot = `colSums(hainich_epi_lichen_div[, 6:601], na.rm = T)`)

#filter anything that is 0 
hai_sum <- dplyr::filter(hai_sum, num_plot > 0 )

hai_sum <- tibble::rownames_to_column(hai_sum)

#order according to abundance
hai_sum <- hai_sum[order(hai_sum$num_plot, decreasing = T ), ]
hai_sum

#####
#where do the 5 most common species occur
#####

#which columns do we need
col_names_hai <- c('EP_PlotID', 
                   'ActivePlot', 'ForMI', 
                   'Porina aenea^E(Fag_syl)', 
                   'Dimerella pineti^E(Fag_syl)',
                   'Lepraria incana^E(Fag_syl)',
                   'Arthonia spadicea^E(Fag_syl)', 
                   'Porina aenea^E(Fra_exc)') 

# select columns of our most abundant species 
hainich_epi_lichen_div_abund <- select(hainich_epi_lichen_div, all_of(col_names_hai))



#####
# Schwaebische Alb
#####

# subset to Schwäbische Alb forests Ep's only

alb_lichen_div <- dplyr::filter(full_lichen_div, grepl('AEW', full_lichen_div$EP_PlotID) )

# only select columns that contain data on epiphytic lichen but keep additional information

tokeep_alb <- colnames(alb_lichen_div[,c(1:3,5,28)])

alb_epi_lichen_div <- dplyr::select(alb_lichen_div, all_of(tokeep_alb),  contains('^E'))

# Lecanora

leconora_abund <- select(alb_epi_lichen_div, all_of(col_names_full),
                         contains("Lecanora") )

rowSums(leconora_abund)
colSums(leconora_abund[6:601,], na.rm = T)

# calculate which species occurs on how many plots 
# first make the pres abs data numeric 

for (i in 6:601) {
  alb_epi_lichen_div[,i] <- as.numeric(alb_epi_lichen_div[,i])
  
}

# calculate column sum
alb_sum <- as.data.frame(colSums(alb_epi_lichen_div[,6:601], na.rm = T))

# rename column 
alb_sum <- alb_sum %>%
  rename(., num_plot = `colSums(alb_epi_lichen_div[, 6:601], na.rm = T)`)

#filter anything that is 0 
alb_sum <- dplyr::filter(alb_sum, num_plot > 0 )

alb_sum <- tibble::rownames_to_column(alb_sum)

#order according to abundance
alb_sum <- alb_sum[order(alb_sum$num_plot, decreasing = T ), ]
alb_sum

#####
#where do the most common species occur
#####

#which columns do we need
col_names_alb <- c('EP_PlotID',  
                   'ActivePlot', 'ForMI', 
                   'Lecanora chlarotera^E(Fag_syl)', 
                   'Phlyctis argena^E(Fag_syl)',
                   'Graphis scripta^E(Fag_syl)',
                   'Pertusaria leioplaca^E(Fag_syl)',
                   'Dimerella pineti^E(Fag_syl)',
                   'Porina aenea^E(Fag_syl)') 

# select columns of our most abundant species 
alb_epi_lichen_div_abund <- select(alb_epi_lichen_div, all_of(col_names_alb))




#####
# Schorfheide-Chorin
#####


# subset to Schwäbische Alb forests Ep's only

sch_lichen_div <- dplyr::filter(full_lichen_div, grepl('SEW', full_lichen_div$EP_PlotID) )

# only select columns that contain data on epiphytic lichen but keep additional information

tokeep_sch <- colnames(sch_lichen_div[,c(1:3,5,28)])

sch_epi_lichen_div <- dplyr::select(sch_lichen_div, all_of(tokeep_sch),  contains('^E'))

# calculate which species occurs on how many plots 
# first make the pres abs data numeric 

for (i in 6:601) {
  sch_epi_lichen_div[,i] <- as.numeric(sch_epi_lichen_div[,i])
  
}

# calculate column sum
sch_sum <- as.data.frame(colSums(sch_epi_lichen_div[,6:601], na.rm = T))

# rename column 
sch_sum <- sch_sum %>%
  rename(., num_plot = `colSums(sch_epi_lichen_div[, 6:601], na.rm = T)`)

#filter anything that is 0 
sch_sum <- dplyr::filter(sch_sum, num_plot > 0 )

sch_sum <- tibble::rownames_to_column(sch_sum)

#order according to abundance
sch_sum <- sch_sum[order(sch_sum$num_plot, decreasing = T ), ]
sch_sum

#####
#where do the most common species occur
#####

#which columns do we need
col_names_sch <- c('EP_PlotID',  
                   'ActivePlot', 'ForMI', 
                   'Lepraria incana^E(Fag_syl)',
                   'Dimerella pineti^E(Fag_syl)',
                   'Porina aenea^E(Fag_syl)',
                   'Lecanora conizaeoides^E(Pin_syl)',
                   'Bacidina arnoldiana agg.^E(Fag_syl)') 

# select columns of our most abundant species 
sch_epi_lichen_div_abund <- select(sch_epi_lichen_div, all_of(col_names_sch))

test <- dplyr::semi_join(hai_sum, alb_sum, by = "rowname")
test2 <- dplyr::semi_join(test, sch_sum, by = "rowname")

test3 <- dplyr::setdiff(hai_sum,alb_sum, "rowname")

dplyr::filter(sch_sum, grepl("Melanelia", rowname))
