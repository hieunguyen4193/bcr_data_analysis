# Introduction

## Input data
Bulk BCR sequencing data analysis repo. 

In this analysis, we use the input `FASTQ` files from `/mnt/storage/data/local/mol_med/bcr/220701_etc_biopsies/samples` in the `molmed` server. A copy of this `FASTQ` dataset is stored in `Coscine` *220701_etc_biopsies_Simons.tar.gz*.

## Metadata

See the file `mid_labels.csv` in the folder `220701_etc_biopsies_data_analysis` in github repo. 

# Run `mixcr` pipeline from `FASTQ` files

- Run the following command to run `mixcr` pipeline: 

```bash 
./mixcr_pipeline/run_mixcr_pipeline.sh -i /path/to/input/FASTQ/files -o /path/to/save/output -e ".fastq"
```

A successful output run is stored at: `molmed-server: /home/hieu/outdir/bcrtree_outdir/20240903/mid_based_output` or `hieu-server: /media/hieunguyen/HNSD01/outdir/220701_etc_biopsies`.

# Run downstream analysis

## Generate trees by the built-in BCR tree algorithm in `mixcr`
### Preprocessing steps

### Main run

# Run `gctree` pipeline 

## Preprocessing steps
### Generate mouse-based trees

### Generate 


