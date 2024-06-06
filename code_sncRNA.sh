# Quality control

## RAW reads - fastqc

sample_dir="Samples"
fqc_raw_dir="QC/RAW_Reads"
ext_v=".fastq.gz"
threads_v=60

mkdir -p $fqc_raw_dir

#fastqc -t 60 $sample_dir/*$ext_v -o $fqc_raw_dir
multiqc ./$fqc_raw_dir -o ./$fqc_raw_dir -n "Raw_Reads_QC"

## RAW reads - adapter remotion and quality filter
cutadapt_dir="Samples_AfterCutadapt"
mkdir -p $cutadapt_dir

for name_r1 in `ls $sample_dir/*_R1*$ext_v`; do 
	outname=$(echo $name_r1 | awk -F"\/" '{print$NF}'); 
        logfile=$(echo $outname | sed "s/$ext_v/.log/g"); 
	echo -e "\nR1=$name_r1\nOutname1=$outname\nlogfile=$logfile";

	#NEB sequence Ref: https://www.neb.com/en-us/faqs/2017/07/17/how-should-my-nebnext-small-rna-library-be-trimmed
	cutadapt -j $threads_v -m 15 -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $name_r1 -o $cutadapt_dir"/"$outname 2>&1 | tee $cutadapt_dir"/"$logfile ;

done;

## Clean reads - fastqc

fqc_clean_dir="QC/Clean_Reads"
mkdir -p $fqc_clean_dir

fastqc -t 60 $cutadapt_dir/*$ext_v -o $fqc_clean_dir
multiqc ./$fqc_clean_dir -o ./$fqc_clean_dir -n "Clean_Reads_QC"

## Alignment Bowtie2/AASRA

ref_file="Rattus_novergicus_SmallncRNA_Ref_IDs.fa"
ref_dir="rat_ncRNA_ref"

echo "Genome dir="$ref_dir
echo "Genome file="$ref_file

init_dir=$(pwd)
# Indexing

cd $ref_dir
saf_f=$(echo $ref_file | sed "s/\.fa//g" | sed "s/\.fasta//g" | sed "s/\.fna//g")".saf"
echo "index_file="$saf_f
AASRA-index -i $ref_file -l CCCCCCCCCC -r GGGGGGGGGG -s $saf_f
cd -

# Alignment
dirout_align="Alignment_AASRA"
mkdir -p $dirout_align

cd $cutadapt_dir
anchor_file_name=$ref_dir"/anchored_"$ref_file
for i in *$ext_v; do
        echo -e "\n\n>>>>>$i<<<<<\n\n" 
	prefix_v=$(echo $i | sed "s/$ext_v//g")
	echo -e "PREFIX="$prefix_v
	echo -e "File_temp="$prefix_v"_temp"
        zcat $i > $prefix_v"_temp"
	echo "genome_path="$ref_dir
	echo "anchor_file_name="$anchor_file_name
        AASRA -p $threads_v -i $prefix_v"_temp" -l CCCCC -r GGGGG -b $init_dir"/"$anchor_file_name  2>&1 | tee ${prefix_v}_aasra.log ;
	#remove fastq files descompacted 
	rm $prefix_v"_temp"
	#remove fastq files descompacted with anchor
	rm "anchored_"$prefix_v"_temp";
	samtools view -b "anchored_${prefix_v}_temp.sam" > ${prefix_v}.bam
	rm "anchored_${prefix_v}_temp.sam"
	mv *.bam $dir_ref/$dirout_align
	mv *aasra.log $dir_ref/$dirout_align
done

cd -

multiqc ./$dirout_align -o ./$dirout_align -n "Alignment_Hisat2_QC"

## FeatureCounts

featurecounts_outdir="Read_Counts"
mkdir -p $featurecounts_outdir
featurecounts_outfile=$featurecounts_outdir"/featureCounts_AASRA.counts"
echo "Outfile="$featurecounts_outfile

featureCounts -T $threads_v -a $ref_dir"/"$saf_f -F SAF -o $featurecounts_outfile $dirout_align/*.bam 2>&1 | tee $featurecounts_outdir'/featureCounts.log'

multiqc ./$featurecounts_outdir -o ./$featurecounts_outdir -n "Assignment_FeatureCounts_QC"





