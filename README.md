# General code-repo for the BCR-bulk and single-cell data analysis.

## Description

- `first_two_scRNA_datasets`: analysis scripts for the first two batches of scRNA + BCR datasets, including `220504_Simons_Pabst_MolMedizin_scVDJseq`, `211216_BSimons`, `211221_BSimons`, `220104_BSimons`. The sub-folder `run_CellRanger_Bonn_data` contains preprocessing scripts and `CellRanger` pipeline for the external datasets (Bonn) `BSimons_Bonn_data`. Raw data and `CellRanger` outputs are stored in NAS at `/volume1/CRC1382_NAS01/CRC1382/OFFICIAL/raw_data`. Note that, for the dataset `220504_Simons_Pabst_MolMedizin_scVDJseq`, we re-run the `CellRanger` pipeline and add **YFP** gene sequence to the mm10 reference genome during alignment steps. 

- `220701_etc_biopsies_data_analysis`: data analysis scripts including `mixcr` pipeline and `gctree` pipeline for the old bulk-BCR dataset. 

- `240805_data_analysis`: data analysis scripts for the new bulk-BCR dataset (240805_BSimons). 

- `240826_data_analysis`: data analysis scripts for the new single-cell BCR VDJ dataset (240826_BSimons).

All raw data are stored in NAS and Coscine RWTH. 





