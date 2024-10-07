gc()
rm(list = ls())

source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/00_helper_functions.R")
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

library("Biostrings")
new.pkgs <- c("msa")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() ==  FALSE){
    BiocManager::install(pkg, update = FALSE)
  }
}
library("msa")

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "1st_2nd_BSimons_Datasets"
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/get_germline_sequences.R")
names(s.V.genes) <- to_vec(
  for (item in names(s.V.genes)) str_split(item, "[*]")[[1]][[1]]
)
names(s.J.genes) <- to_vec(
  for (item in names(s.J.genes)) str_split(item, "[*]")[[1]][[1]]
)

names(cellranger.ref) <- to_vec(
  for (item in names(cellranger.ref)) str_split(str_split(item, " ")[[1]][[1]], "[|]")[[1]][[2]]
)
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
##### get metadata sheets and single cell data from the 1st and 2nd datasets
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "1st_2nd_BSimons_Datasets"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.06.output <- file.path(path.to.main.output, "06_output")
dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
s.obj.1st <- readRDS(file.path(path.to.main.output, 
                               "01_output",
                               "1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1", 
                               "1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1.rds"))
metadata.1st <- s.obj.1st@meta.data %>% rownames_to_column("barcode")

s.obj.2nd <- readRDS(file.path(path.to.main.output, 
                               "01_output", 
                               "2nd_dataset_removed_5_6.without_reInt.res1", 
                               "2nd_dataset_removed_5_6.without_reInt.res1.rds"))
metadata.2nd <- s.obj.2nd@meta.data %>% rownames_to_column("barcode")

all.samples <- c(unique(s.obj.1st$name), unique(s.obj.2nd$name))

#####----------------------------------------------------------------------#####
##### clone dataframes
#####----------------------------------------------------------------------#####
path.to.clone.dfs <- file.path(path.to.main.output, "VDJ_data_output")
clonedf <- data.frame()
for (sample.id in all.samples){
  tmpdf <- readxl::read_excel(file.path(path.to.clone.dfs, sprintf("%s.xlsx", sample.id)))
  tmpdf$SampleID <- sample.id
  clonedf <- rbind(clonedf, tmpdf)
}

nt_cols <- to_vec(
  for (item in colnames(clonedf)) if (grepl("_nt", item) == TRUE) item
)
clonedf.raw <- clonedf 

##### keep IGH clonedf only
clonedf <- subset(clonedf, clonedf$chain == "IGH")

i <- 1
V.gene <- clonedf[i, ]$v_gene
J.gene <- clonedf[i, ]$j_gene
CDR3.seq <- clonedf[i, ]$cdr3_nt
CDR3.length <- nchar(CDR3.seq)

# GL.V.gene <- s.V.genes[[V.gene]] %>% as.character()
# GL.J.gene <- s.J.genes[[J.gene]] %>% as.character()

GL.V.gene <- cellranger.ref[[V.gene]] %>% as.character()
GL.J.gene <- cellranger.ref[[J.gene]] %>% as.character()

repN.seq <- paste(replicate(n = CDR3.length, expr = "N"), collapse = "")
GL.seq <- sprintf("%s%s%s", GL.V.gene, repN.seq, GL.J.gene)
clone.seq <- clonedf[i, ]$prep.full.seq

input.seqs <- c(clone.seq, GL.seq) %>% DNAStringSet()
msa.output <- msa(input.seqs, method = "Muscle")
aligned.clone.seq <- toString(unmasked(msa.output)[[1]])
aligned.GL.seq <- toString(unmasked(msa.output)[[2]])



         