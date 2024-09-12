##### batch 1
path_to_input="/mnt/storage/data/local/mol_med/bcr/220701_etc_biopsies/samples";
##### batch 2:
# path_to_input="/home/hieu/storage/240826_Simons_pabst_MolMedcine_AmpliSeq/Fastq";

output_version=20240903
outputdir=/home/hieu/outdir/bcrtree_outdir/${output_version};
mkdir -p $outputdir;

files=$(ls ${path_to_input}/*R1*.fastq*);
for file in $files;do \
filename=$(echo $file | xargs -n 1 basename);
sampleid=${filename%_S*};
fastq1=${file}
fastq2=${file%R1*}R2_001.fastq; # add .gz depending on the input fastq format ext
echo -e "------------------------------------------"
echo -e "WORKING ON SAMPLE " $sampleid "\n"
echo -e "------------------------------------------"
echo $sampleid;
echo -e $fastq1;
echo -e $fastq2;

bash mixcr_pipeline.sh $sampleid $fastq1 $fastq2 $outputdir;
done