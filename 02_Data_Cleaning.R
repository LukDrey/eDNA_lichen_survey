#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")

library(dplyr); packageVersion("dplyr")

library(phyloseq); packageVersion("phyloseq")

library(stringr); packageVersion("stringr")

library(GGally); packageVersion("GGally")

##################################################################
##                          Section 2                           ##
##         Data Loading, Data Cleaning and Initialization       ##
##################################################################

##----------------------------------------------------------------
##                           Metadata                            -
##----------------------------------------------------------------

# Load the metadata table.
own_metadata <- utils::read.csv(here::here("Data","sample_data_2.csv"), header = T, sep = ",")

# Get the EP_Plot_ID from the Sample names
own_metadata <- own_metadata %>% tidyr::separate(sample, c(NA, 'Plot_ID', NA) , sep = '_', remove = F)

# Insert EW after the frist character of the string in a given cell. 
own_metadata$Plot_ID <- base::gsub("^(.{1})(.*)$",
                             "\\1EW\\2",
                             own_metadata$Plot_ID)

# Load the Climate data. 
climate <- utils::read.csv(here::here("Data","plots_climate_weekly.csv"), header = T, sep = ",")

# Subset to the sampling weeks and 2 weeks before for each Explo. Alb = weeks 16-18,
# Hai = weeks 17-19, Sch = weeks 18-20. 
sampling_weeks_clim <- climate %>% 
  dplyr::filter(base::grepl("AEW", plotID) & datetime %in% c("2021W16", "2021W17", "2021W18") |
                base::grepl("HEW", plotID) & datetime %in% c("2021W17", "2021W18", "2021W19") |
                base::grepl("SEW", plotID) & datetime %in% c("2021W18", "2021W19", "2021W20"))

# Remove the datetime column since it is not needed anymore.
sampling_weeks_clim$datetime <- NULL

# Calculate the average of the climate variable for each plot. 
sampling_weeks_clim_avg <- sampling_weeks_clim %>% 
  dplyr::group_by(plotID) %>%  
  dplyr::summarise_each(dplyr::funs(base::round(base::mean(., na.rm = T), 2))) %>% 
  dplyr::rename(Plot_ID = plotID)

# One of the plots (HEW33) had a broken thermometer so we replace the temperature 
# value with that of the nearest station (HEW26).
sampling_weeks_clim_avg$Ta_200[sampling_weeks_clim_avg$Plot_ID == "HEW33"] <-
  sampling_weeks_clim_avg$Ta_200[sampling_weeks_clim_avg$Plot_ID == "HEW26"]

# Load the stand structure dataset. 
stand_structure_full <-  utils::read.csv(here::here("Data","22766_3_data.csv"), header = T, sep = ";")

# Subset to the variables we are interested in and give them more meaningful names. 
stand_structure_variables <- stand_structure_full %>% 
  dplyr::select(EP_Plotid, ssm_N, ssm_BA, sp_BA_1D, d_m, d_SD, d_gini) %>% 
  dplyr::rename(Plot_ID = EP_Plotid, 
                stand_density_abundance = ssm_N,
                stand_density_basal_area = ssm_BA,
                stand_evenness_basal_area = sp_BA_1D,
                DBH_avg = d_m)

# Load the dataset of the effective number of forest layers = measure of vertical heterogeneity. 
vertical_heterogeneity <- utils::read.csv(here::here("Data","27826_3_data.csv"), header = T, sep = ";")

# We are only interested in the newest data. 
vertical_heterogeneity_new <- vertical_heterogeneity %>% 
  dplyr::select(EP, enl_2019) %>% 
  dplyr::rename(Plot_ID = EP)

# Load the dataset describing the proportion of forest in a 2 km radius around the plot to see 
# if the plot is situated within a forest or is a smaller patch of trees. 
surroundings <- utils::read.csv(here::here("Data","15929_2_data.csv"), header = T, sep = ";")

# Subset to only keep info on surrounding forest areas. 
surroundings_forest <- surroundings %>% 
  dplyr::select(Plot, RA_forest) %>% 
  dplyr::rename(Plot_ID = Plot)

# Load the Canopy Openness dataset. 
canopy_openness <- utils::read.csv(here::here("Data","canopy_openness_2019.csv"), header = T, sep = ",") %>% 
  dplyr::rename(Plot_ID = EP)

# Load in the inventory data of single trees.
single_trees <- utils::read.csv(here::here('Data', '21426_3_data.csv'), header = T, sep = ';')

# Calculate plot wise number of present tree species. 
plot_struc <- single_trees %>% 
  dplyr::group_by(EP_Plotid) %>% 
  dplyr::count(species) %>% 
  dplyr::rename(number = n, Plot_ID = EP_Plotid)

# Find the most abundant tree species on each plot.
plot_dom <- plot_struc %>% 
  dplyr::group_by(Plot_ID) %>%
  dplyr::mutate(num_dom = max(number)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(number == num_dom) %>%
  dplyr::select(-number) %>%
  dplyr::rename(dom_species = species)

# Get the total number of trees on each plot. 
plot_trees <- plot_struc %>% 
  dplyr::group_by(Plot_ID) %>%
  dplyr::summarise(num_tot = sum(number))

# Combine both into one table. 
plot_comp <- dplyr::left_join(plot_dom, plot_trees, by = 'Plot_ID')

# Calculate most abundant tree proportion. 
plot_comp <- plot_comp %>%
  dplyr::mutate(dom_tot_ratio = base::round(num_dom/num_tot , 2))

# Remove intermediate columns we are not interested in. 
plot_comp_final <- plot_comp %>% 
  dplyr::select("Plot_ID", "dom_tot_ratio")

# Load the geographical distances between plots. 
################################################################################################
####################################################################################
###############################################################################
#################################################################
geo <- fread(here::here("Data", "basic_plot_info.csv"))
geo <- geo[Landuse == "Forest", .(EP_Plot_ID, Latitude, Longitude)]
geo <- geo[EP_Plot_ID != "na"]
setnames(geo, old = "EP_Plot_ID", new = "Plot")

# Combine all of our metadata into one big metadata table. 
metadata_full <- plyr::join_all(base::list(own_metadata,
                                     sampling_weeks_clim_avg,
                                     stand_structure_variables,
                                     vertical_heterogeneity_new, 
                                     surroundings_forest, 
                                     canopy_openness, 
                                     plot_comp_final),
                                by = "Plot_ID",
                                type = "inner")

# Set the sample IDs as the rownames.
base::rownames(metadata_full) <- metadata_full$sample
metadata_full$sample <- NULL

# Load in the algal ASV IDs. 
algae_asv_IDs <- base::readRDS(here("Data", "algae_asv_IDs.rds"))

# Get the sample names from the algal ASV table. 
sample_names <- base::colnames(algae_asv_IDs) %>% 
  stringr::str_subset(pattern = ("sequence|MPC|PCN|BL"), negate = T)

# Some Plots had been clear cut and need to be removed. 
# An additional one has very low read counts for algae and fungi and also needs to be removed.
metadata_full <- metadata_full[!(metadata_full$Plot_ID == "HEW3" |
                                   metadata_full$Plot_ID == "HEW13" |
                                   metadata_full$Plot_ID == "HEW12"),]

##----------------------------------------------------------------
##                          ASV tables                           -
##----------------------------------------------------------------

##---------
##  Fungi  
##---------
ASV_table_fungi_cur <- base::readRDS(here("Data", "ASV_table_fungi_cur.rds"))

# Keep only samples that do represent real tree swabs. Cut Controls. 
asv_fungi <- ASV_table_fungi_cur$curated_table %>%
  dplyr::select(all_of(sample_names))

##---------------------------------------------------------------
##                        Taxonomy Tables                       -
##---------------------------------------------------------------

##---------
##  Fungi  
##---------

# Load the fungal taxonomy table.
# (Available as supplementary data)
tax_fungi <- base::readRDS(here::here("Data", 'tax_table_fungi.rds'))
tax_fungi <- base::as.data.frame(tax_fungi) %>%
  tibble::rownames_to_column('sequence')
tax_fungi <- tax_fungi %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the taxonomy table and the representative sequences
tax_clean_fungi <- dplyr::left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi')

# Split the taxonomy into different columns of taxonomic levels.

fungi_tax_fin <- tidyr::separate(tax_clean_fungi, Kingdom, c(NA, 'Kingdom') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Phylum, c(NA, 'Phylum') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Class, c(NA, 'Class') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Order, c(NA, 'Order') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Family, c(NA, 'Family') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Genus, c(NA, 'Genus') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Species, c(NA, 'Species') , sep = '__')

# Rename the ASV_ID column. 
fungi_tax_fin <- dplyr::rename(fungi_tax_fin, ASV_ID = seq_name_fungi)

# Set rownames.
base::row.names(fungi_tax_fin) <- fungi_tax_fin$ASV_ID

fungi_tax_fin$sequence_fungi <- NULL
fungi_tax_fin$ASV_ID <- NULL

##---------------------------------------------------------------
##               Create the Phyloseq Objects                    -
##---------------------------------------------------------------

##---------
##  Fungi  
##---------

# Transform dataframe to matrix.
asvmat_fungi <- base::as.matrix(asv_fungi)

# Transform dataframe to matrix.
taxmat_fungi <- base::as.matrix(fungi_tax_fin) 

# Create ASV table for phyloseq.
ASV_FUN <- phyloseq::otu_table(asvmat_fungi, taxa_are_rows = T) 

# Create taxonomy table for phyloseq.
TAX_FUN <- phyloseq::tax_table(taxmat_fungi) 

# Metadata for phyloseq. 
sampledata <- phyloseq::sample_data(metadata_full) 

# Combine in phyloseq object. 
phy_fungi <- phyloseq::phyloseq(ASV_FUN, TAX_FUN, sampledata) 
phy_fungi

# Remove taxa that have zero reads remaining. 
phy_fungi <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy_fungi) != 0, phy_fungi)

########
# Make a table that is usable with FUNGUILD. 
########

# Pull out the taxonomy table of the curated fungal ASVs. 
tax_fun_prep <- base::data.frame(phyloseq::tax_table(phy_fungi)) %>% 
  tibble::rownames_to_column("ASV_ID")

# Create an empty data frame to fill with the format FUNGUILD requires. 
tax_funguild <- base::data.frame(ASV_ID = tax_fun_prep$ASV_ID)

tax_funguild$taxonomy <- base::paste(tax_fun_prep$Kingdom, tax_fun_prep$Phylum, tax_fun_prep$Class, tax_fun_prep$Order,
                               tax_fun_prep$Family, tax_fun_prep$Genus, tax_fun_prep$Species, sep = ";")

utils::write.table(tax_funguild, here::here("Data", "tax_funguild.tsv"), sep = "\t")

# Funguild will be excuted on the server again since it is a Python Script. 

# These are the results and we can append them back to our phyloseq object. 
# For the fungi we can append the data obtained from the functional annotation with FUNGUILD. 
funguild <- utils::read.csv(here::here( "Data", "tax_funguild.guilds.txt"), header = T, sep = "\t")

# Keep only the trophic mode since that is what we are most interested in. 
guilds <- funguild %>% 
  dplyr::select(ASV_ID, Guild) 

# Join with the taxonomy file. 
fungi_tax_fin_guilds <- dplyr::left_join(fungi_tax_fin %>% 
                                            tibble::rownames_to_column( var = "ASV_ID"),
                                         guilds)

# Set the rownames.
base::rownames(fungi_tax_fin_guilds) <- fungi_tax_fin_guilds$ASV_ID 
fungi_tax_fin_guilds$ASV_ID <- NULL 

# Get the guild info into the taxonomy table. 

phyloseq::tax_table(phy_fungi) <- phyloseq::tax_table(base::as.matrix(fungi_tax_fin_guilds))


# Now we can extract a phyloseq object with just lichenized fungi and
# another containing the rest. 
# Look for rows in the column Guild that contain info on Lichenized OTUs. 
lichen_otu_names <- funguild %>%
  filter(stringr::str_detect(Guild, 'Lichenized')) %>%
  .$ASV_ID 

# Extract these from the raw dataset.
sub_lichen_otu <- subset(otu_table(phy_fungi), 
                         rownames(otu_table(phy_fungi)) %in% lichen_otu_names)

#merge into a phyloseq object
phy_lichen <- merge_phyloseq(sub_lichen_otu, 
                                tax_table(phy_fungi), 
                                sample_data(phy_fungi))
# Subset the phyloseq object to split the dataset into the bark and soil samples. 

##---------
##  Lichen  
##---------

phy_lichen_bark <- phyloseq::subset_samples(phy_lichen, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

phy_lichen_soil <- phyloseq::subset_samples(phy_lichen, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)


#################################################################
##                          Section 5                          ##
##                   Miscallenous Data Export                  ##
#################################################################

###
# Species Lists
###

lichen_species_list <- metagMisc::phyloseq_to_df(phy_lichen_bark, addtax = T, sorting = "taxonomy") %>% 
  dplyr::select("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Guild") %>% 
  dplyr::arrange(Class)

write.csv(lichen_species_list, "lichen_species_list.csv")
