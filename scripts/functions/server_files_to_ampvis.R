####################
# Load the results from 16S from the server to R 
##################################

server_files_to_ampvis <- function(
    DataPath = "C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_RETHINK_sample_collection/data/amplicon_data/",
    metadata_file_name = "JUN2023_samples-metadata.xlsx",
    folder = "OnlineKontrol_2023_JUNE/",  
    append_barcode = "DK_2306",
    ID_col = "SampleID", 
    total_reads_file_name = "totalreads",
    metadata_file_path = paste0(DataPath, "metadata/", metadata_file_name) 
  ){


# Adapted form Kaspers pipeline but more simpel
 #library("ampvis2")
 #library("data.table")

# Download data from server at set directories 
metadata_file <- metadata_file_path #paste0(DataPath, "metadata/", metadata_file_name)    # <--- adjust this

###################################### 

#nanopore data
if (total_reads_file_name == "totalreads" ){
    reads <- fread(
    paste0(DataPath, folder, "totalreads.csv"),
    header = FALSE,
    col.names = c(paste0(ID_col), "Reads")
    )
  } else if(total_reads_file_name == "summary"){
    reads <- vroom::vroom(
      paste0(DataPath, folder, "summary.txt"), 
      #header = FALSE,
      #col.names = c(paste0(ID_col), "Reads")
    ) %>% 
      select(barcode, total_reads#, mapped_reads, nfilt
             ) %>% 
      data.table::as.data.table() %>% 
      mutate(
        #mapped_reads = mapped_reads - nfilt
        #total_reads = total_reads + nfilt
        ) #%>% select(-nfilt)
  colnames(reads) <- c(paste0(ID_col), "total_reads"#, "mapped_reads"
                       )
  }else{
  print("Total reads / Summary file is missing")
    }

metadata_ext <- tools::file_ext(metadata_file)
if (metadata_ext %in% c("xls", "xlsx")) {
  metadata <- openxlsx::read.xlsx(
    metadata_file,
    detectDates = TRUE
  )
  data.table::setDT(metadata)
} else if (metadata_ext %in% c("txt", "csv", "tsv")) {
  metadata <- fread(
    metadata_file,
    encoding = "UTF-8",
    nrows = 24
  )
}

# Append run specific ID 
metadata[[ID_col]] <- with(metadata, paste0(append_barcode, "_", Barcode))
reads[[ID_col]] <- with(reads, paste0(append_barcode, "_", reads[[ID_col]]))

# Merge the reads info from summery file to the metadata 
metadata <- reads[metadata, on = ID_col] 
# Filter samles not in the summery file
metadata <- metadata %>% filter(!is.na(total_reads))

# Impart the OTU table WITH  where unclassifed 
mappedreads <- vroom::vroom(paste0(DataPath, folder,"otutable_mappedreads.tsv"), 
                     show_col_types = FALSE) 

# Get the current column names
current_colnames <- colnames(mappedreads)
# Do not rename the OTU and Taxonomy columns 
columns_to_rename <- c(2:(ncol(mappedreads) - 7))
# Rename the selected columns by adding the prefix
new_colnames <- paste(append_barcode, current_colnames[columns_to_rename], sep = "_")
# Update the column names of the data frame
colnames(mappedreads)[columns_to_rename] <- new_colnames

# Load ampvis with the raw reads 
d_raw <- amp_load(
  mappedreads,
  metadata = metadata
)

# this assumes order of samples are the same, but only one sample this time
  # unclassified <- d_raw$metadata$Read - colSums(d_raw$abund)
  # d_norm <- d_raw
  # d_norm$abund["unclassified", ] <- unclassified
  # d_norm$tax["unclassified", ] <- "unclassified"


# import of nomalised amp_object (cannot be rarefied)
norm_abun <- vroom::vroom(paste0(DataPath, folder,"otutable_normalised.tsv"), 
                     show_col_types = FALSE) 
# Get the current column names
current_colnames <- colnames(norm_abun)
# Do not rename the OTU and Taxonomy columns 
columns_to_rename <- c(2:(ncol(norm_abun) - 7))
# Rename the selected columns by adding the prefix
new_colnames <- paste(append_barcode, current_colnames[columns_to_rename], sep = "_")
# Update the column names of the data frame
colnames(norm_abun)[columns_to_rename] <- new_colnames


d_norm <- amp_load(
  norm_abun,
  metadata = metadata
)

# d_noctrl <- amp_subset_samples(
#   d_g,
#   !grepl("ctrl", tolower(PlantID)),
#   !grepl("pcrpos|pcrneg", tolower(LibID)),
#   !grepl("ctrl", tolower(SampleSite)),
#   !grepl("ctrl", tolower(SampleName))
# )

# unclassified <- 
#   unclassified %>% enframe(., name = {{ID_col}}, value = "unclassified_reads")

# tot_reads_with_unclassified <- 
#   amp_alpha_diversity(d_raw) %>% 
#   left_join(unclassified, by = join_by({{ID_col}})) %>% 
#   #select(ID_col, LibID,SampleName,ProjectName, seq_ID, Reads, unclassified_reads) %>% 
#   rename("classified_reads" = "Reads")




# Assign names
list <- list(d_norm, d_raw#, tot_reads_with_unclassified
             )

}



















