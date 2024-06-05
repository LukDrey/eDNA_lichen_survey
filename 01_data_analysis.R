#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")
# ‘1.0.1’

library(tidyverse); packageVersion("tidyverse")
# ‘2.0.0’

library(phyloseq); packageVersion("phyloseq")
# ‘1.44.0’

library(stringr); packageVersion("stringr")
# ‘1.5.0’

library(GGally); packageVersion("GGally")
# 2.1.2

library(sf); packageVersion("sf")
#1.0.14

library(ggspatial); packageVersion("ggspatial")
#1.1.9

##################################################################
##                          Section 2                           ##
##         Data Loading, Data Cleaning and Initialization       ##
##################################################################

##----------------------------------------------------------------
##                           Metadata                            -
##----------------------------------------------------------------

# Load the metadata table that includes the dominant tree species, substrate it was
# collected from and the Sample ID.
metadata <- utils::read.csv(here::here("Data","sample_data_2.csv"),
                            header = T, sep = ",") %>% 
  # Get the Plot ID from the Sample ID.
  tidyr::separate(sample, c(NA, 'Plot_ID', NA),
                                                 sep = '_', remove = F)

# Insert EW after the first character of the string in a given cell.
# With this we can later match them to the floristic study. 
metadata$Plot_ID <- base::gsub("^(.{1})(.*)$",
                             "\\1EW\\2",
                             metadata$Plot_ID)

# Set the sample IDs as the rownames and remove the now unnecessary column. 
base::rownames(metadata) <- metadata$sample
metadata$sample <- NULL

# Some Plots had been clear cut and need to be removed. 
# An additional one has very low read counts for fungi
# and thus also needs to be removed.
metadata <- metadata[!(metadata$Plot_ID == "HEW3" |
                       metadata$Plot_ID == "HEW13" |
                       metadata$Plot_ID == "HEW12"),]

##----------------------------------------------------------------
##                          ASV tables                           -
##----------------------------------------------------------------

# Load the fungal ASV table that was previously curated with the LULU algorithm. 
ASV_table_fungi_cur <- base::readRDS(here("Data", "ASV_table_fungi_cur.rds"))$curated_table

# Get the names of the real samples without the negative controls and blanks
# from the fungal ASV table. 
sample_names <- base::colnames(ASV_table_fungi_cur) %>% 
  stringr::str_subset(pattern = ("sequence|MPC|PCN|BL"), negate = T)

# Keep only samples that do represent real tree swabs. Cut Controls. 
asv_fungi <- ASV_table_fungi_cur %>%
  dplyr::select(all_of(sample_names))

##---------------------------------------------------------------
##                        Taxonomy Tables                       -
##---------------------------------------------------------------

# Load the fungal taxonomy table.
# (Available as supplementary data in the Github repository).
tax_fungi <- base::readRDS(here::here("Data", 'tax_table_fungi.rds')) %>% 
  base::as.data.frame() %>%
  tibble::rownames_to_column('sequence') %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads from a FASTA file.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
fungi_rep_seqs <- base::data.frame(seq_name_fungi = base::names(fungi_seqs_fasta),
                                   sequence_fungi = base::paste(fungi_seqs_fasta))

# Join the taxonomy table and the representative sequences
tax_fungi_names <- dplyr::left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi')

# Split the taxonomy into different columns of taxonomic levels.
fungi_tax_fin <- tidyr::separate(tax_fungi_names, Kingdom, c(NA, 'Kingdom') , sep = '__') %>% 
  tidyr::separate(., Phylum, c(NA, 'Phylum') , sep = '__') %>% 
  tidyr::separate(., Class, c(NA, 'Class') , sep = '__') %>% 
  tidyr::separate(., Order, c(NA, 'Order') , sep = '__') %>% 
  tidyr::separate(., Family, c(NA, 'Family') , sep = '__') %>% 
  tidyr::separate(., Genus, c(NA, 'Genus') , sep = '__') %>% 
  tidyr::separate(., Species, c(NA, 'Species') , sep = '__')

# We have some assignments to plants as well, since they also carry ITS copies.
# We want to keep only Fungi. 
fungi_tax_fin <- fungi_tax_fin %>% 
  dplyr::filter(Kingdom == "Fungi") %>% 
  dplyr::rename(, ASV_ID = seq_name_fungi)

# Set the ASV IDs as the rownames as neeeded for phyloseq.
base::row.names(fungi_tax_fin) <- fungi_tax_fin$ASV_ID

# Remove information that we do not need anymore going further. 
fungi_tax_fin$sequence_fungi <- NULL
fungi_tax_fin$ASV_ID <- NULL

##---------------------------------------------------------------
##               Create the Phyloseq Object                     -
##---------------------------------------------------------------

# Create ASV table for phyloseq. 
# Note that the function requires a matrix as input so we need to transform our dataframe.
ASV_FUN <- phyloseq::otu_table(base::as.matrix(asv_fungi),
                               taxa_are_rows = T) 

# Create taxonomy table for phyloseq.
TAX_FUN <- phyloseq::tax_table(base::as.matrix(fungi_tax_fin)) 

# Metadata for phyloseq. 
sampledata <- phyloseq::sample_data(metadata) 

# Combine in phyloseq object. 
phy_fungi <- phyloseq::phyloseq(ASV_FUN,
                                TAX_FUN,
                                sampledata) 

# Some taxa have zero because they were only found in cut samples. 
# Remove those taxa.
phy_fungi <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy_fungi) != 0,
                                  phy_fungi)
phy_fungi

# Subset the phyloseq object so it contains only bark samples. 
phy_fungi_bark <- phyloseq::subset_samples(phy_fungi, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

##---------------------------------------------------------------
##            Prepare fungal sequences for Martin7              -
##---------------------------------------------------------------

# Instead of the already existing fungal assignments with the UNITE database,
# we want to build our taxonomy on a specialized lichen database 
# which is called Martin7 (https://doi.org/10.23855/preslia.2023.311).
# We use a simple nucleotide BLAST with a 97% identity cutoff
# against the fasta file. 

# First we need to prepare a FASTA file that we can blast.
pure_fun_seqs <- data.frame(tax_table(phy_fungi_bark)) %>% 
  tibble::rownames_to_column(var = "ASV_ID") %>% 
  dplyr::select(ASV_ID) %>% 
  dplyr::left_join(., dplyr::rename(fungi_rep_seqs, ASV_ID = seq_name_fungi)) %>% 
  dplyr::rename(Sequence = sequence_fungi)

# Write a FASTA file for the ASVs so we can blast them against 
# the Martin7 database.
# First we need to define a function that lets us do that. 
writeFasta <- function(data, filename, sequence_col, name_col){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum, name_col], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum, sequence_col]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Then write the FASTA file.
# writeFasta(data = pure_fun_seqs,
#            filename = here("Data", "all_seqs.fa"),
#            name_col = "ASV_ID",
#            sequence_col = "Sequence")

# BLAST is a standalone software that we run through a Linux server. 
# Here is the commands we used. 
# makeblastdb -in m7_ITS.fas -parse_seqids -dbtype nucl
# 
# blastn \ 
# -db m7_ITS.fas \
# -outfmt '6 qseqid sseqid pident' \
# -out all_M7.txt \
# -qcov_hsp_perc 80 \
# -perc_identity 97 \
# -num_threads 20 \
# -query all_seqs.fa

# Load in the blast results. 
# Filter the eDNA dataset to only include the ASVs that were assigned to lichens by MARTIN7
full_martin7_blast <- read.table(here("Data", "all_M7.txt"),
                               header = F,
                               sep = "") %>% 
  dplyr::rename(ASV_ID = V1,
                Species = V2,
                perc_ident = V3) 

# Filter to keep only the top results in term of identity. 
filt_martin7_blast <- full_martin7_blast %>% 
  dplyr::group_by(ASV_ID) %>% 
  dplyr::slice(which.max(perc_ident)) %>% 
  dplyr::arrange(Species)

# Clean the Species column and split it into Genus and Species. 
filt_martin7_blast_clean <- cbind(filt_martin7_blast, 
                                  clean = stringr::str_extract(filt_martin7_blast$Species, "[^_]*_[^_]*")) %>% 
  dplyr::select("ASV_ID", "clean", "perc_ident") %>% 
  tidyr::separate_wider_delim(., clean, delim = "_", names = c("Genus", "Species"))

# Join the martin7 data frame with the taxonomy table from UNITE. 
unite_tax <- data.frame(tax_table(phy_fungi_bark)) %>%  
  tibble::rownames_to_column(var = "ASV_ID") %>% 
  dplyr::select(ASV_ID, Genus, Species) %>% 
  dplyr::rename(Genus_unite = Genus,
                Species_unite = Species)

##---------------------------------------------------------------
##                     FUNGUILD preparation                     -
##---------------------------------------------------------------

# FUNGUILD is a python script that matches our ASVs against a
# database and assigns information to which guild the fungus belongs along
# with some functional information. 
# We are using it to find which ASVs are assigned as lichenized based
# on the taxonomic assignment with UNITE. 

# Get the taxonomy table of the curated fungal ASVs. 
tax_fun_prep <- base::data.frame(phyloseq::tax_table(phy_fungi)) %>% 
  tibble::rownames_to_column("ASV_ID")

# Create an empty data frame to fill with the format FUNGUILD requires. 
tax_funguild <- base::data.frame(ASV_ID = tax_fun_prep$ASV_ID)

tax_funguild$taxonomy <- base::paste(tax_fun_prep$Kingdom,
                                     tax_fun_prep$Phylum,
                                     tax_fun_prep$Class,
                                     tax_fun_prep$Order,
                                     tax_fun_prep$Family,
                                     tax_fun_prep$Genus,
                                     tax_fun_prep$Species,
                                     sep = ";")

# utils::write.table(tax_funguild,
#                   here::here("Data", "tax_funguild.tsv"),
#                   sep = "\t",
#                   row.names = F, 
#                   quote = F)

# Funguild will be executed on the server. 
# Here are the commands. 
# python Guilds_v1.0.py -otu tax_funguild.tsv  -db fungi

# Append funguild information so we can exclude Unite assignments that funguild
# does not classify as lichenized. 
# For the fungi we can append the data obtained from the functional annotation with FUNGUILD. 
funguild <- utils::read.csv(here::here( "Data", "tax_funguild.guilds.txt"),
                            header = T, sep = ",")

# Keep only the trophic mode since that is what we are most interested in. 
guilds <- funguild %>% 
  dplyr::select(ASV_ID, Guild) 

# Join with the UNITE taxonomy file. 
unite_tax_guilds <- dplyr::left_join(unite_tax, guilds)

# Remove taxa that are not assigned as lichens by FUNGUILD. 
unite_tax_lichens <- unite_tax_guilds %>% 
  filter(stringr::str_detect(Guild, 'Lichenized'))

# Remove ASVs from the UNITE file that are already assigned by Martin7. 
unite_lichen_not_martin7 <- unite_tax_lichens %>% 
  filter(!(ASV_ID %in% filt_martin7_blast_clean$ASV_ID)) %>% 
  # Rename the columns so we can join the df with the martin7 df. 
  rename(Genus = Genus_unite,
         Species = Species_unite) %>% 
  mutate(unite = "yes")

# Join with blast results against Martin7. 
combined_tax <- dplyr::full_join(filt_martin7_blast_clean,
                                unite_lichen_not_martin7,
                                by = c("ASV_ID", "Genus", "Species")) %>% 
  mutate(across(unite,  ~replace(., is.na(.), "no")))

##---------------------------------------------------------------
##                 Data Cleaning of eDNA results                -
##---------------------------------------------------------------

# Calculate the prevalence of each ASV and put it in a table with
# the taxonomic information as well. 
# NOTE: The function prevalence from the metagMisc is not exported to the namespace, 
# hence we need the ::: to get it. 
eDNA_abund_prev <- metagMisc:::prevalence(phy_fungi_bark) %>% 
  tibble::rownames_to_column(var = "ASV_ID") %>% 
  dplyr::select("ASV_ID", "Prevalence", "TotalAbundance") %>% 
  # Add a column stating if the ASV is found in the eDNA.
  dplyr::mutate(in_eDNA = "yes") %>% 
  dplyr::rename(Count = TotalAbundance,
                num_of_plots_edna = Prevalence) %>% 
  dplyr::left_join(combined_tax, .) %>% 
  dplyr::select(-dplyr::one_of("Guild")) %>% 
  dplyr::arrange(Genus)

# write.csv(eDNA_abund_prev, "lichen_edna_martin7.csv", row.names = F)

##---------------------------------------------------------------
##          Data Cleaning of floristic mapping results          -
##---------------------------------------------------------------

# Read in the data from Steffen Boch which can be found at https://www.bexis.uni-jena.de/ under Dataset ID 4460.  
# This is a species by sample matrix.
lichen_matrix  <- utils::read.table(here::here('Data',
                                               'lichen_diversity_forests',
                                               'MatrixData.txt'),
                                    header = T, sep = '|')

# Transform the matrix data into a tidy data format.  
lichen_df <- data.frame(t(lichen_matrix)) %>% 
  tibble::rownames_to_column() %>% 
  # Make the first row the header.
  janitor::row_to_names(.,1, remove_row = T) %>% 
  dplyr::rename(., PlotID = Plot_ID)

# Read in the basic plot info that can be found here https://www.bexis.uni-jena.de/ under Dataset ID 20826.
# And directly remove information that is not interesting for our study. 
plot_info <- read.table(here('Data', '20826.txt'), header = T, sep = '\t') %>% 
  dplyr::select("EP_PlotID", "PlotID", "Exploratory")

# Join to get a data frame containing the lichen presence and plot information. 
full_lichen_div <- dplyr::left_join(plot_info, lichen_df, by = 'PlotID')

# Remove the information for plots that are not within the 
# 150 experimental plots. 
forest_lichen_div <- dplyr::filter(full_lichen_div,
                                   grepl('HEW|AEW|SEW',
                                         full_lichen_div$EP_PlotID)) %>% 
  # Turn the presence/absence columns to numeric.
  dplyr::mutate(dplyr::across(!EP_PlotID & !PlotID & !Exploratory, as.numeric)) %>% 
  # Replace NAs with 0 for absence.
  mutate(across(everything(),  ~replace(., is.na(.), 0)))

# Pivot the table to longer format.
forest_lichen_div_long <- tidyr::pivot_longer(forest_lichen_div,
                                              cols = colnames(forest_lichen_div 
                                                              %>% select(where(is.numeric))),
                             names_to = "Species", values_to = "Occurence")

# Find the names of the lichens that were found epiphytic within the 150 EPs.
epiphyte_records <- forest_lichen_div_long %>% 
  filter(stringr::str_detect(Species, '\\^E'),
         Occurence > 0)

epiphyte_records$Species <- gsub("\\^E.*", "", epiphyte_records$Species) 

epiphyte_record_species <- unique(epiphyte_records$Species)

# Remove the information on which substrate the species was found.
forest_lichen_div_long$Species <- gsub("\\^E.*", "", forest_lichen_div_long$Species) %>% 
  gsub("\\^G.*", "",.) %>% 
  gsub("\\^T.*", "", .) %>% 
  gsub("\\^B.*", "", .) %>% 
  gsub("\\^M.*", "", .) 

# Calculate how often the species was observed in each plot 
# (number of tree species the species was observed on).
lichen_substrate_occur <- forest_lichen_div_long %>%
  group_by(Species, EP_PlotID) %>%
  summarise(Occurence = sum(Occurence))

# If a species was observed on more than one tree species we set that number
# to 1 to signify the presence in the plot. We are not interested in which tree.
lichen_substrate_presence <- lichen_substrate_occur %>%
  mutate(across(Occurence, ~1 * (. != 0)))

# Keep only species that were found epiphytically in the study.
epiphyte_substrate_presence <- lichen_substrate_presence %>% 
  dplyr::filter(Species %in% epiphyte_record_species)

# Calculate the number of plots a species was observed in. 
lichen_plot_presence <- lichen_substrate_presence %>%
  group_by(Species) %>%
  summarise(num_of_plots_flor = sum(Occurence))

# Keep only species that were found epiphytically in the study
# and were observed in at least one of the 150 experimental plots. 
epiphyte_plot_presence <- lichen_plot_presence %>% 
  dplyr::filter(Species %in% epiphyte_record_species,
                num_of_plots_flor > 0)

# Split the species column into two columns Genus and Species
# based on the space in the middle. 
epiphyte_plot_presence_clean <- epiphyte_plot_presence %>% 
  tidyr::separate_wider_regex(cols = Species,
                              patterns = c(Genus = "(?:^Cf. .*)*(?:.*)*(?: )*(?!:agg.)", " ",
                                           Species = "(?: )*(?:.*)")) %>% 
  tidyr::separate_wider_regex(cols = Genus, 
                              patterns = c(Genus = "(?:^Cf. \\w+)*(?:\\w+)",
                                           Helper = ".*"), 
                              too_few = "align_start")  %>%
  dplyr::mutate(Helper = stringr::str_remove(Helper, "^[\\s]+")) %>%
  dplyr::mutate(Species = stringr::str_remove(Species, "^[\\s]+")) %>% 
  tidyr::unite("Species", Helper:Species, sep = " ") %>%
  dplyr::mutate(Species = stringr::str_remove(Species, "^[\\s]+")) %>% 
  dplyr::mutate(Species = stringr::str_remove(Species, "^[.]")) %>% 
  dplyr::mutate(Genus = stringr::str_replace_all(Genus, "Cf$", "Cf.")) %>%
  dplyr::mutate(Species = stringr::str_remove(Species, "^[\\s]+")) %>% 
  dplyr::mutate(in_floristic = "yes")

##---------------------------------------------------------------
##                      Combining both studies                  -
##---------------------------------------------------------------

# Combine the floristic mapping results with the eDNA dataframe.
combined_species_df <- dplyr::full_join(eDNA_abund_prev, epiphyte_plot_presence_clean,
                                        by = c("Genus", "Species")) %>% 
  dplyr::filter((Genus != "Umbilicaria") %>% tidyr::replace_na(TRUE)) %>% 
  dplyr::mutate(in_eDNA = tidyr::replace_na(in_eDNA, "no"),
                in_floristic = tidyr::replace_na(in_floristic, "no"),
                num_of_plots_flor = tidyr::replace_na(num_of_plots_flor, 0)) 

# Save the table.
# write.csv(combined_species_df, "combined_species_df_martin7.csv", row.names = F)

##---------------------------------------------------------------
##             BLASTing ASVs that have no assignment            -
##---------------------------------------------------------------

# Some ASVs have no taxonomic assignment from the UNITE database. 
# For additional assignments we query the NCBI nt database using BLASTn.
# For that we need a FASTA file of the unknown ASVs.

# Add the sequences so we can blast taxa that are NA for Genus or Species. 
combined_species_df_seqs <- combined_species_df %>% 
  dplyr::left_join(., dplyr::rename(fungi_rep_seqs, ASV_ID = seq_name_fungi)) %>% 
  dplyr::rename(Sequence = sequence_fungi)

unknown_asv_seqs <- combined_species_df_seqs %>%
  filter(is.na(Species) | Species %in% c("sp.", "sp.1", "cf.")) %>%
  dplyr::select(ASV_ID, Sequence)

# Write a FASTA file for the unknown species so we can blast them.
# writeFasta(data = unknown_asv_seqs,
#            filename = here("Data", "unknown_asv_seqs.fa"),
#            name_col = "ASV_ID",
#            sequence_col = "Sequence")

##---------------------------------------------------------------
##                   Taxonomic Assessment                       -
##---------------------------------------------------------------

# We cleaned the combined species table and updated the taxonomy to the
# names accepted as "current" in the MycoBank database. 
# We also added any additional assignment that we could make via a BLAST search. 
# This was the case for 2 ASVs. 
# Furthermore we evaluated if the species assignments make sense for our study
# region. This unfortunately could only be done by hand outside of R. 
# All following steps are based on this cleaned and curated dataset. 

##---------------------------------------------------------------
##                   Last data cleaning steps                   -
##---------------------------------------------------------------

# Read in the cleaned and curated combined species dataframe. 
cleaned_df <- read.csv("combined_species_df_martin7.csv", header = T, sep = ",") 

# Some of the records do not have a species assignment that we can be 
# 100% certain of. We remove these records to be conservative in our analysis. 
# This is true for both floristic records and eDNA assignments. 
cleaned_df_final <-  cleaned_df %>%
  filter(!is.na(Species)) %>%  
  filter(str_detect(Species, "sp\\.|sp.1|cf.|spec.", negate = T)) %>% 
  filter(str_detect(Genus, "Cf.", negate = T))

#################################################################
##                          Section 3                          ##
##                        Data Analysis                        ##
#################################################################

##---------------------------------------------------------------
##               Find species with multiple ASVs                -
##---------------------------------------------------------------

eDNA_duplicate_species <- cleaned_df_final %>% 
  dplyr::filter(!is.na(ASV_ID)) %>%
  group_by(Genus, Species) %>%
  arrange(.by_group = T) %>% 
  filter(n()>1) %>%
  ungroup() %>% 
  tidyr::unite(Gen_Spec, c( "Genus", "Species"), sep = " ")

# Find the number of species with several ASVs assigned to it. 
num_dup_species <- length(unique(eDNA_duplicate_species$Gen_Spec))
num_dup_species

##---------------------------------------------------------------
##                        Distribution Maps                     -
##---------------------------------------------------------------

##----------------------
##  Data Preparation  
##----------------------

# We have identified 5 species of interest that have varying degrees of overlap between 
# eDNA and floristic survey. We want to plot them. 

# Read in the shapefiles provided by the Biodiversity Exploratories.
# Available at https://www.bexis.uni-jena.de/ dataset ID 31234.
SCH_map <- sf::read_sf(here("Data", "31234_8_Exploratorium_sch.gpkg"))
ALB_map <- sf::read_sf(here("Data", "31234_8_Exploratorium_alb.gpkg"))
HAI_map <- sf::read_sf(here("Data", "31234_8_Exploratorium_hai.gpkg"))

# Load the coordinates for the plots. 
# Available at https://www.bexis.uni-jena.de/ dataset ID 1000.
geo <- utils::read.csv(here::here("Data", "basic_plot_info.csv"), sep = ";") %>% 
  dplyr::filter(Landuse == "Forest") %>% 
  dplyr::filter(EP_Plot_ID != "na") %>% 
  dplyr::rename(Plot = EP_Plot_ID) %>% 
  dplyr::select(Plot, Latitude, Longitude)

# Some of the floristic records have updated names in MycoBank
# We kept a record of the old names so we can match the records
# of the old names to their new names. 

# Find the species that had a name change. 
cleaned_df_names <- cleaned_df_final %>% 
  dplyr::select(ASV_ID, Genus, Species, old_genus, old_species) %>% 
  dplyr::filter(is.na(ASV_ID)) %>% 
  dplyr::filter(old_genus != "")  %>% 
  tidyr::unite(Gen_Spec, c( "Genus", "Species"), sep = " ") %>% 
  tidyr::unite(old_Gen_Spec, c( "old_genus", "old_species"), sep = " ") %>% 
  dplyr::select(-ASV_ID)

# Change the names of the species in the dataframe that contains the records 
# of presence/absence of the floristic records. 
epiphyte_substrate_presence_new <- epiphyte_substrate_presence %>%
  dplyr::rename(old_Gen_Spec = Species) %>% 
  dplyr::full_join(., cleaned_df_names) %>% 
  dplyr::mutate(new_Gen_Spec = dplyr::case_when(is.na(Gen_Spec) ~ old_Gen_Spec,
                                            TRUE ~ Gen_Spec)) %>%
  dplyr::ungroup() %>%
  dplyr::select(new_Gen_Spec, EP_PlotID, Occurence) %>% 
  dplyr::filter(str_detect(new_Gen_Spec, "sp\\.|sp.1|cf.|spec.|Cf.", negate = T))

# Separate the combined Genus and Species column. 
epiphyte_substrate_presence_clean <- epiphyte_substrate_presence_new %>% 
  tidyr::separate(col = new_Gen_Spec,
                  into = c("Genus", "Species"),
                  sep = " ")

# Get the original ASV table from the phyloseq object and 
# match it to the ASV IDs from the final cleaned taxonomy dataset. 
eDNA_ASV_samples <- subset(otu_table(phy_fungi_bark),
                         rownames(otu_table(phy_fungi_bark)) %in%
                           c(cleaned_df_final$ASV_ID)) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "ASV_ID") %>% 
  dplyr::left_join(cleaned_df_final %>% 
                     dplyr::select(ASV_ID, Genus, Species)) %>% 
  dplyr::relocate(ASV_ID, Genus, Species) %>% 
  tidyr::pivot_longer(cols = dplyr::starts_with("Sample"),
                      names_to = "Sample_ID",
                      values_to = "Read_Count") 

# Now get the plotID into the table so we can compare eDNA and
# floristic in the distribution maps. 
sample_plot_match <- data.frame(phyloseq::sample_data(phy_fungi_bark)) %>% 
  tibble::rownames_to_column(var = "Sample_ID") %>% 
  dplyr::select(Sample_ID, Plot_ID, exploratory)

# In a last step turn the count data into presence/absence. 
edna_plot_occur <- eDNA_ASV_samples %>% 
  # and subset to the 5 species of interest. 
  dplyr::filter(Genus %in% c("Buellia", "Graphis", "Lepraria",
                              "Phlyctis", "Physcia")) %>% 
  dplyr::filter(Species %in% c("incana", "griseovirens", "scripta", 
                               "argena", "adscendens")) %>%   
  dplyr::group_by(Genus, Species, Sample_ID) %>% 
  dplyr::summarise(sum = sum(Read_Count)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate_if(is.numeric, ~1 * (. != 0)) %>% 
  dplyr::left_join(sample_plot_match) %>% 
  dplyr::rename(eDNA_occur = sum)

# Combine floristic with the eDNA subset. Add the geographic information.
combined_plot_occur <- edna_plot_occur %>% 
  dplyr::left_join(epiphyte_substrate_presence_clean %>%
                     dplyr::rename(Plot_ID = EP_PlotID)) %>% 
  dplyr::left_join(geo %>% dplyr::rename(Plot_ID = Plot)) %>% 
  data.frame() %>% 
  dplyr::rename(floristic_occur = Occurence) %>% 
  dplyr::filter(!is.na(.$floristic_occur)) 

combined_plot_occur$eDNA_occur <- if_else(combined_plot_occur$eDNA_occur == 1,
                                          "yes", "no") %>% as.factor()
combined_plot_occur$floristic_occur <- if_else(combined_plot_occur$floristic_occur == 1,
                                               "yes", "no") %>% as.factor()

# Add a column for easier plotting of the different possibilities
# with which method a species was found. 
combined_plot_occur$combination <- paste(combined_plot_occur$eDNA_occur,
                                         sep="-",
                                         combined_plot_occur$floristic_occur) %>%
  as.factor()

# Subset to the regions.
combined_plot_occur_SCH <- combined_plot_occur %>% 
  dplyr::filter(exploratory == "Schorfheide")  
combined_plot_occur_ALB <- combined_plot_occur %>%
  dplyr::filter(exploratory == "Alb")
combined_plot_occur_HAI <- combined_plot_occur %>%
  dplyr::filter(exploratory == "Hainich")

# Subset to each of the different species of interest.
genera <- c("Physcia", "Buellia",  "Graphis",  "Lepraria", "Phlyctis")

# Schorfheide-Chorin
for (i in genera) {
  df <- subset(combined_plot_occur_SCH, Genus == i)
  assign(paste("combined_plot_occur_SCH", i, sep = "_") , df)
}

############
# Hainich-Dün
for (i in genera) {
  df <- subset(combined_plot_occur_HAI, Genus == i)
  assign(paste("combined_plot_occur_HAI", i, sep = "_") , df)
}

############
# Swabian Alb
for (i in genera) {
  df <- subset(combined_plot_occur_ALB, Genus == i)
  assign(paste("combined_plot_occur_ALB", i, sep = "_") , df)
}

##----------------------
##      Plotting  
##----------------------

# All our plots for one region will follow the same structure and will be plotted the same
# so we can define a function of the plotting for each region and make the script more readable. 

plot_SCH <- function(map, point_data){
  ggplot2::ggplot() +
    geom_sf(data = map, fill = "#91be1c", alpha = 0.75) +
    geom_point(aes(x=Longitude, y=Latitude,
                   color = combination, fill = combination),
               data = point_data, size = 1, shape = 21, stroke = 0.2,
               position = position_jitter(h=0.01, w=0.01, seed = 42)) +
    scale_color_manual(values = c("no-no" = "black",
                                  "no-yes" = "#251605",
                                  "yes-no" = "#BFFFF1",
                                  "yes-yes" = "#508484"),
                       labels = c("in none",
                                  "only floristic",
                                  "only eDNA",
                                  "in both"),
                       name = "Species occurence", 
                       drop = FALSE) +
    scale_fill_manual(values = c("no-no" = 0,
                                 "no-yes" = "#251605",
                                 "yes-no" = "#BFFFF1",
                                 "yes-yes" = "#508484"),
                      labels = c("in none",
                                 "only floristic",
                                 "only eDNA",
                                 "in both"),
                      name = "Species occurence", 
                      drop = FALSE) +
    coord_sf(default_crs = sf::st_crs(4326)) +
    theme_void() + 
    theme(title = element_blank())
}

plot_ALB <- function(map, point_data, title){
  ggplot2::ggplot() +
    geom_sf(data = map, fill = "#91be1c", alpha = 0.75)+
    geom_point(aes(x=Longitude, y=Latitude,
                   color = combination, fill = combination),
               data = point_data, size = 1, shape = 21, stroke = 0.2,
               position = position_jitter(h=0.0015, w=0.0015, seed = 42)) +
    scale_color_manual(values = c("no-no" = "black",
                                  "no-yes" = "#251605",
                                  "yes-no" = "#BFFFF1",
                                  "yes-yes" = "#508484"),
                       labels = c("in none",
                                  "only floristic",
                                  "only eDNA",
                                  "in both"),
                       name = "Species occurence", 
                       drop = FALSE) +
    scale_fill_manual(values = c("no-no" = 0,
                                 "no-yes" = "#251605",
                                 "yes-no" = "#BFFFF1",
                                 "yes-yes" = "#508484"),
                      labels = c("in none",
                                 "only floristic",
                                 "only eDNA",
                                 "in both"),
                      name = "Species occurence", 
                      drop = FALSE) +
    coord_sf(default_crs = sf::st_crs(4326)) +
    theme_void() +
    labs(title = title) +
    theme(title = element_text(size = 7, face = "italic"),
          plot.title = element_text(hjust = 0.5))  +
    guides(color = guide_legend(override.aes = list(size = 3, stroke = 1) ) )
}

plot_HAI <- function(map, point_data){
  ggplot2::ggplot() +
    geom_sf(data = map, fill = "#91be1c", alpha = 0.75) +
    geom_point(aes(x=Longitude, y=Latitude,
                   color = combination, fill = combination),
               data = point_data, size = 1, shape = 21, stroke = 0.2,
               position = position_jitter(h=0.015, w=0.015, seed = 42)) +
    scale_color_manual(values = c("no-no" = "black",
                                  "no-yes" = "#251605",
                                  "yes-no" = "#BFFFF1",
                                  "yes-yes" = "#508484"),
                       labels = c("in none",
                                  "only floristic",
                                  "only eDNA",
                                  "in both"),
                       name = "Species occurence", 
                       drop = FALSE) +
    scale_fill_manual(values = c("no-no" = 0,
                                 "no-yes" = "#251605",
                                 "yes-no" = "#BFFFF1",
                                 "yes-yes" = "#508484"),
                      labels = c("in none",
                                 "only floristic",
                                 "only eDNA",
                                 "in both"),
                      name = "Species occurence", 
                      drop = FALSE) +
    coord_sf(default_crs = sf::st_crs(4326)) +
    theme_void()  + 
    theme(title = element_blank())
}
##----------------------
##  Schorfheide-Chorin  
##----------------------

###################### 
# Physcia adscendens
physcia_plot_SCH <- plot_SCH(map = SCH_map, point_data = combined_plot_occur_SCH_Physcia)
physcia_plot_SCH

###################### 
# Buellia griseovirens
buellia_plot_SCH <- plot_SCH(map = SCH_map, point_data = combined_plot_occur_SCH_Buellia) +
  ggspatial::annotation_scale(width_hint = 0.1, location = "tr", text_cex = .6,
                              height = unit(1, "mm"))
buellia_plot_SCH

###################### 
# Graphis scripta
graphis_plot_SCH <- plot_SCH(map = SCH_map, point_data = combined_plot_occur_SCH_Graphis)
graphis_plot_SCH

###################### 
# Lepraria incana
lepraria_plot_SCH <- plot_SCH(map = SCH_map, point_data = combined_plot_occur_SCH_Lepraria)
lepraria_plot_SCH

###################### 
# Phlyctis argena
phlyctis_plot_SCH <- plot_SCH(map = SCH_map, point_data = combined_plot_occur_SCH_Phlyctis)
phlyctis_plot_SCH


##----------------------
##  Hainich-Dün  
##----------------------

###################### 
# Physcia adscendens
physcia_plot_HAI <- plot_HAI(map = HAI_map, point_data = combined_plot_occur_HAI_Physcia)
physcia_plot_HAI

###################### 
# Buellia griseovirens
buellia_plot_HAI <- plot_HAI(map = HAI_map, 
                             point_data = combined_plot_occur_HAI_Buellia) +
  ggspatial::annotation_scale(width_hint = 0.12, location = "tr", text_cex = .6,
                              height = unit(1, "mm"))
buellia_plot_HAI

###################### 
# Graphis scripta
graphis_plot_HAI <- plot_HAI(map = HAI_map, point_data = combined_plot_occur_HAI_Graphis)
graphis_plot_HAI

###################### 
# Lepraria incana
lepraria_plot_HAI <- plot_HAI(map = HAI_map, point_data = combined_plot_occur_HAI_Lepraria)
lepraria_plot_HAI

###################### 
# Phlyctis argena
phlyctis_plot_HAI <- plot_HAI(map = HAI_map, point_data = combined_plot_occur_HAI_Phlyctis)
phlyctis_plot_HAI

##----------------------
##  Swabian Alb  
##----------------------

###################### 
# Physcia adscendens
physcia_plot_ALB <- plot_ALB(map = ALB_map,
                             point_data = combined_plot_occur_ALB_Physcia,
                             title = "Physcia adscendens")
physcia_plot_ALB 

###################### 
# Buellia griseovirens
buellia_plot_ALB <- plot_ALB(map = ALB_map,
                             point_data = combined_plot_occur_ALB_Buellia,
                             title = "Buellia griseovirens") +
  ggspatial::annotation_scale(width_hint = 0.15, location = "tr", text_cex = .6,
                              height = unit(1, "mm"))
buellia_plot_ALB 

###################### 
# Graphis scripta
graphis_plot_ALB <- plot_ALB(map = ALB_map,
                             point_data = combined_plot_occur_ALB_Graphis,
                             title = "Graphis scripta") 
graphis_plot_ALB 

###################### 
# Lepraria incana
lepraria_plot_ALB <- plot_ALB(map = ALB_map,
                              point_data = combined_plot_occur_ALB_Lepraria,
                              title = "Lepraria incana") 
lepraria_plot_ALB 

###################### 
# Phlyctis argena
phlyctis_plot_ALB <- plot_ALB(map = ALB_map,
                              point_data = combined_plot_occur_ALB_Phlyctis,
                              title = "Phlyctis argena") 
phlyctis_plot_ALB 

##----------------------
##  Combine the plots   
##----------------------
# First extract the overall legend for all plots
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  physcia_plot_ALB + theme(legend.position = "bottom")
)

# Combine the maps.
# Add  the legend to the plots and control the relative size and also position of the legend 
# relative to the plot, i.e. decrease the plots margin to move the legend closer. 
legend_maps <- cowplot::plot_grid(buellia_plot_ALB + theme(legend.position="none"),
                                  graphis_plot_ALB + theme(legend.position="none"),
                                  lepraria_plot_ALB + theme(legend.position="none"),
                                  phlyctis_plot_ALB + theme(legend.position="none"),
                                  physcia_plot_ALB + theme(legend.position="none"),
                                  buellia_plot_HAI + theme(legend.position="none"),
                                  graphis_plot_HAI + theme(legend.position="none"),
                                  lepraria_plot_HAI + theme(legend.position="none"),
                                  phlyctis_plot_HAI + theme(legend.position="none"),
                                  physcia_plot_HAI + theme(legend.position="none"),
                                  buellia_plot_SCH + theme(legend.position="none"),
                                  graphis_plot_SCH + theme(legend.position="none"),
                                  lepraria_plot_SCH + theme(legend.position="none"),
                                  phlyctis_plot_SCH + theme(legend.position="none"),
                                  physcia_plot_SCH + theme(legend.position="none"),
                                  NULL, NULL, legend, 
                                  nrow = 4, ncol = 5, 
                                  rel_heights = c(1,1,1,0.1))
legend_maps
 
# Save the figure as a pdf.
# cowplot::save_plot("final_maps_new.pdf", legend_maps,
#                   base_width = 6.88976, base_height = NULL, bg = "white")
# 
# 
# cowplot::save_plot("final_maps_new.png", legend_maps,
#                   base_width = 11.29, base_height = NULL, bg = "white")

##----------------------------------------------------------------
##                        Venn Diagrams                          -
##          Species overlap between the two methods              -
##----------------------------------------------------------------

# First we need a dataframe that only contains each species once and the info
# on if it was found in the eDNA and/or the floristic study.
# We also add a helper column that allows for easier plotting of the venn diagram.  
venn_df <- cleaned_df_final %>% 
  group_by(Genus, Species) %>% 
  summarise(across(c("num_of_plots_edna", "num_of_plots_flor"), ~sum(., na.rm = TRUE))) %>%  
  dplyr::mutate(in_eDNA = if_else(num_of_plots_edna > 0, "yes", "no"),
                in_floristic = if_else(num_of_plots_flor > 0, "yes", "no")) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(helper_col = paste("helper", 1:nrow(.), sep = ""))

in_eDNA <- venn_df %>% 
  dplyr::filter(in_eDNA == "yes") 

in_floristic <- venn_df %>% 
  dplyr::filter(in_floristic == "yes") 

survey_list <- list(
  in_eDNA = in_eDNA$helper_col,
  in_floristic = in_floristic$helper_col)

names(survey_list) <- c("eDNA metabarcoding", "floristic survey")

overlap_venn <- ggvenn::ggvenn(survey_list, 
               stroke_size = 0.5,
               set_name_size = 4,
               stroke_color = "black",
               fill_color = c("#BFFFF1", "#EEA243"),
               fill_alpha = 1, 
               auto_scale = TRUE)
overlap_venn

# ggplot2::ggsave(here::here("Figures", "overlap_venn_new.png"), plot = overlap_venn, device = png,
#                 width = 175, height = 150, units = "mm", bg = "white")

##----------------------------------------------------------------
##                Bar Charts of species prevalence               -
##----------------------------------------------------------------

cleaned_df_final_eDNA <- cleaned_df_final %>% 
  dplyr::filter(!is.na(.$ASV_ID)) %>% 
  tidyr::unite(Species,  Genus:Species, sep = " ")

eDNA_ASV_presence_table <- subset(otu_table(phy_fungi_bark),
                                  rownames(otu_table(phy_fungi_bark)) %in%
                                    c(cleaned_df_final$ASV_ID)) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "ASV_ID") %>% 
  dplyr::left_join(cleaned_df_final %>% 
                     dplyr::select(ASV_ID, Genus, Species)) %>% 
  dplyr::relocate(ASV_ID, Genus, Species) %>% 
  tidyr::unite(Species,  Genus:Species, sep = " ")

# Replace the read count with a 1 for to signify presence. 
eDNA_ASV_presence_table[3:154] <- replace(eDNA_ASV_presence_table[,3:154],
                                          eDNA_ASV_presence_table[,3:154] > 0, 1)

# Split into a list of dataframes with one dataframe for each species. 
species_ASV_presence <- split(eDNA_ASV_presence_table,
                              eDNA_ASV_presence_table$Species)

# Use a for loop to calculate the number of plots we found each species 
# in with the eDNA methods. 
eDNA_species_prevalence <- data.frame()
for (i in names(species_ASV_presence)) {
  col_sums <- colSums(species_ASV_presence[[i]][3:154])
  sum <- sum(replace(col_sums, col_sums > 0, 1))
  species <- i
  df <- data.frame(species = species, num_of_plots_eDNA = sum)
  eDNA_species_prevalence <- rbind(eDNA_species_prevalence , df)
}

floristic_species_prevalence <- epiphyte_plot_presence %>%
  dplyr::rename(old_Gen_Spec = Species) %>% 
  dplyr::full_join(., cleaned_df_names) %>%  
  dplyr::mutate(new_Gen_Spec = dplyr::case_when(is.na(Gen_Spec) ~ old_Gen_Spec,
                                                TRUE ~ Gen_Spec)) %>% 
  dplyr::ungroup() %>%
  dplyr::select(new_Gen_Spec, num_of_plots_flor) %>% 
  dplyr::filter(str_detect(new_Gen_Spec, "sp\\.|sp.1|cf.|spec.|Cf.", negate = T)) %>% 
  dplyr::rename(species = new_Gen_Spec)

# Combine the two datasets. 
combined_species_prevalence <- full_join(floristic_species_prevalence,
                                         eDNA_species_prevalence) %>% 
  tidyr::replace_na(list(num_of_plots_flor = 0,
                         num_of_plots_eDNA = 0)) 

# Subset it to the species that occur in more than 25 plots in one of the two methods.
combined_species_prevalence_over_25 <- combined_species_prevalence %>% 
  dplyr::filter(num_of_plots_eDNA >= 25 | num_of_plots_flor >= 25)


combined_species_prevalence_over_25_plot <- reshape2::melt(combined_species_prevalence_over_25[,c('species','num_of_plots_eDNA','num_of_plots_flor')],id.vars = 1)

combined_species_prevalence_over_25_plot$species <- as.factor(combined_species_prevalence_over_25_plot$species)

ordered <- combined_species_prevalence_over_25[order(combined_species_prevalence_over_25$num_of_plots_eDNA,
                                                          decreasing = T),]

names_axis <- c(ordered$species)

prevalence_plot <- ggpubr::ggbarplot(data = combined_species_prevalence_over_25_plot,
                  x = "species", y = "value",
                  fill = "variable",
                  position = position_dodge(0.7)) +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(face = "italic", size = 7),
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 7)) +
  scale_fill_manual(values = c("#BFFFF1", "#EEA243"),
                    labels = c("eDNA metabarcoding", "floristic survey")) +
  scale_x_discrete(limits = rev(names_axis),
                   expand = c(0.006, 0)) +
  labs(y = "Number of occurences") +
  coord_flip()
prevalence_plot

# ggplot2::ggsave(here::here("Figures", "prevalence_plot_new.png"), plot = prevalence_plot, device = png,
#                 width = 175, height = 150, units = "mm", bg = "white")


##----------------------------------------------------------------
##                    Calculate lichen percentage                -
##----------------------------------------------------------------

# Subset to the regions. 

phy_fungi_bark_SCH <- phyloseq::subset_samples(phy_fungi_bark,
                                                       exploratory == "Schorfheide") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

phy_fungi_bark_ALB <- phyloseq::subset_samples(phy_fungi_bark,
                                               exploratory == "Alb") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

phy_fungi_bark_HAI <- phyloseq::subset_samples(phy_fungi_bark,
                                               exploratory == "Hainich") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Subset to the lichen ASVs.

lichen_ASV_count_table_SCH <- phy_fungi_bark_SCH %>% 
  phyloseq::otu_table() %>% 
  data.frame() %>% 
  subset(.,
         rownames(.) %in%
           c(cleaned_df_final$ASV_ID)) %>% 
  as.data.frame()

lichen_ASV_count_table_ALB <- phy_fungi_bark_ALB %>% 
  phyloseq::otu_table() %>% 
  data.frame() %>% 
  subset(.,
         rownames(.) %in%
           c(cleaned_df_final$ASV_ID)) %>% 
  as.data.frame()

lichen_ASV_count_table_HAI <- phy_fungi_bark_HAI %>% 
  phyloseq::otu_table() %>% 
  data.frame() %>% 
  subset(.,
         rownames(.) %in%
           c(cleaned_df_final$ASV_ID)) %>% 
  as.data.frame()

# Calculate the percentage of lichens in the total reads.

eDNA_total_reads_SCH <- sum(colSums(otu_table(phy_fungi_bark_SCH)))
eDNA_total_reads_ALB <- sum(colSums(otu_table(phy_fungi_bark_ALB)))
eDNA_total_reads_HAI <- sum(colSums(otu_table(phy_fungi_bark_HAI)))

lichen_reads_SCH <- sum(colSums(lichen_ASV_count_table_SCH))
lichen_reads_ALB <- sum(colSums(lichen_ASV_count_table_ALB))
lichen_reads_HAI <- sum(colSums(lichen_ASV_count_table_HAI))

(lichen_reads_SCH / eDNA_total_reads_SCH) * 100
(lichen_reads_ALB / eDNA_total_reads_ALB) * 100
(lichen_reads_HAI / eDNA_total_reads_HAI) * 100

##----------------------------------------------------------------
##             Number of floristic species per region            -
##----------------------------------------------------------------

epiphyte_records_SCH <- subset(epiphyte_records, Exploratory == "SCH")  %>%   
  filter(str_detect(Species, "sp\\.|sp.1|cf.|spec.|Cf.", negate = T)) %>% 
  dplyr::filter(Species %in% epiphyte_record_species,
                Occurence > 0)

epiphyte_records_ALB <- subset(epiphyte_records, Exploratory == "ALB")  %>%   
  filter(str_detect(Species, "sp\\.|sp.1|cf.|spec.|Cf.", negate = T)) %>% 
  dplyr::filter(Species %in% epiphyte_record_species,
                Occurence > 0)

epiphyte_records_HAI <- subset(epiphyte_records, Exploratory == "HAI")  %>%   
  filter(str_detect(Species, "sp\\.|sp.1|cf.|spec.|Cf.", negate = T)) %>% 
  dplyr::filter(Species %in% epiphyte_record_species,
                Occurence > 0)

# Get the number of species.
length(unique(epiphyte_records_SCH$Species))

length(unique(epiphyte_records_ALB$Species))

length(unique(epiphyte_records_HAI$Species))

##----------------------------------------------------------------
##                      Species List export                      -
##----------------------------------------------------------------

species_list <- combined_species_prevalence %>% 
  tidyr::separate_wider_delim(., species,
                              delim = " ",
                              names = c("Genus", "Species"),
                              too_many = "merge") %>% 
  dplyr::select(Genus, Species,
                num_of_plots_eDNA,
                num_of_plots_flor) %>% 
  dplyr::mutate(in_eDNA = if_else(num_of_plots_eDNA > 0, "yes", "no"),
                in_floristic = if_else(num_of_plots_flor > 0, "yes", "no"))

# write.csv(species_list,
#           "species_list_martin7.csv", row.names = F)

