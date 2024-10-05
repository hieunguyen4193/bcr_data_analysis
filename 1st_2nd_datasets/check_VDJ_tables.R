gc()
rm(list = ls())

source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/00_helper_functions.R")
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

library("Biostrings")

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "1st_2nd_BSimons_Datasets"

#####----------------------------------------------------------------------#####
##### read file from VDJ output, after pre-processing. 
#####----------------------------------------------------------------------#####
path.to.VDJ.ouptut <- file.path(outdir, PROJECT, "VDJ_output")
all.VDJ.files <- Sys.glob(file.path(path.to.VDJ.ouptut, "annotated_contigs*.csv"))
names(all.VDJ.files) <- to_vec(
  for (item in all.VDJ.files) str_replace(str_replace(basename(item), "annotated_contigs_clonaltype_", ""), ".csv", "")
)
path.to.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.save.output <- file.path(path.to.output, "VDJ_data_output")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### read file "filtered_contig_annotations.csv" from the raw CellRanger pipeline
#####----------------------------------------------------------------------#####
path.to.main.input <- "/media/hieunguyen/GSHD_HN01/storage/1st_2nd_BSimons_Datasets/VDJ"
all.raw.VDJ.files <- Sys.glob(file.path(path.to.main.input, "*", "filtered_contig_annotations.csv"))
names(all.raw.VDJ.files) <- to_vec(
  for (item in all.raw.VDJ.files) basename(dirname(item))
)
all.raw.VDJ.files <- all.raw.VDJ.files[names(all.VDJ.files)]
all.raw.VDJ.fasta <- Sys.glob(file.path(path.to.main.input, "*", "filtered_contig.fasta"))
names(all.raw.VDJ.fasta) <- to_vec(
  for (item in all.raw.VDJ.fasta) basename(dirname(item))
)
all.raw.VDJ.fasta <- all.raw.VDJ.fasta[names(all.VDJ.files)]

#####----------------------------------------------------------------------#####
##### merge pre-processed data with the raw VDJ data
#####----------------------------------------------------------------------#####
maindf <- list()
for (sampleid in names(all.VDJ.files)){
  rawdf <- read.csv(all.raw.VDJ.files[[sampleid]])
  prepdf <- read.csv(all.VDJ.files[[sampleid]]) 
  prepdf$prep_barcode <- prepdf$barcode
  prepdf <- prepdf %>% rowwise() %>%
    mutate(barcode = str_replace(barcode, sprintf("%s_%s_", sampleid, sampleid), ""))
  prepdf <- subset(prepdf, select = c(barcode, CTgene, CTnt, CTaa, CTstrict))
  tmpdf1 <- merge(rawdf, prepdf, by.x = "barcode", by.y = "barcode")
  tmp.fasta <- readDNAStringSet(all.raw.VDJ.fasta[[sampleid]])
  tmpdf2 <- as.data.frame(tmp.fasta) %>% rownames_to_column("contig_id")
  colnames(tmpdf2) <- c("contig_id", "fasta_seq")
  final.tmpdf <- merge(tmpdf1, tmpdf2, by.x = "contig_id", by.y = "contig_id")
  maindf[[sampleid]] <- final.tmpdf
  writexl::write_xlsx(maindf[[sampleid]], file.path(path.to.save.output, sprintf("%s.xlsx", sampleid)))
}

##### Check fasta files
tmpdf <- maindf$`17_MM9_Ecoli`

to_vec(for (item in colnames(tmpdf))
  if (grepl("_nt", item)) item
  )





