# Quality control

## RAW reads QC - fastqc

sample_dir="Samples"
fqc_raw_dir="QC/RAW_Reads"
ext_v=".fastq.gz"
threads_v=60

mkdir -p $fqc_raw_dir

fastqc -t 60 $sample_dir/*$ext_v -o $fqc_raw_dir
multiqc ./$fqc_raw_dir -o ./$fqc_raw_dir -n "Raw_Reads_QC"

## RAW reads - adapter remotion and quality filter
cutadapt_dir="Samples_AfterCutadapt"
mkdir -p $cutadapt_dir

for name_r1 in `ls $sample_dir/*_R1*$ext_v`; do 
	name_r2=$(echo $name_r1 | sed "s/_R1/_R2/g"); 
	outname=$(echo $name_r1 | awk -F"\/" '{print$NF}'); 
	outname2=$(echo $outname | sed "s/_R1/_R2/g"); 
        logfile=$(echo $outname | sed "s/$ext_v/.log/g"); 
	echo -e "\nR1=$name_r1\nR2=$name_r2\nOutname1=$outname\nOutname2=$outname2\nlogfile=$logfile";

    cutadapt -j $threads_v -m 20 -q 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $cutadapt_dir"/"$outname -p $cutadapt_dir"/"$outname2 $name_r1 $name_r2 2>&1 | tee $cutadapt_dir"/"$logfile ;
done;

## Clean reads QC - fastqc
fqc_clean_dir="QC/Clean_Reads"
mkdir -p $fqc_clean_dir

fastqc -t 60 $cutadapt_dir/*$ext_v -o $fqc_clean_dir
multiqc ./$fqc_clean_dir -o ./$fqc_clean_dir -n "Clean_Reads_QC"

## Alignment using Hisat2
genome_file="genome/Rattus_norvegicus_Rnor_6p0.fa"
echo "Genome file="$genome_file

# Indexing
index_f=$(echo $genome_file | sed "s/\.fa//g" | sed "s/\.fasta//g" | sed "s/\.fna//g")"_hisat2index"
echo "index_file="$index_f
hisat2-build -p $threads_v $genome_file $index_f

dirout_align="Alignment"
mkdir -p $dirout_align

for map_r1 in `ls $cutadapt_dir/*_R1*.gz`; do
	echo "Read1="$map_r1;
	map_r2=$(echo $map_r1 | sed "s/_R1/_R2/g");
	echo "Read2="$map_r2;
	
	outname_prefix=$dirout_align"/"$(basename $map_r1 | sed "s/_R1//g" | sed "s/$ext_v//g");
	echo "prefix ="$outname_prefix

	{ hisat2 -p $threads_v --dta -x $index_f -1 $map_r1 -2 $map_r2 2> >(tee $outname_prefix".hisat2.log" >&2); } | tee >(samtools flagstat - > $outname_prefix".flagstat") | samtools sort -@ 4 -O BAM | tee $outname_prefix".bam" | samtools index - $outname_prefix".bam.bai";

done;

multiqc ./$dirout_align -o ./$dirout_align -n "Alignment_Hisat2_QC"

## FeatureCounts
gtf_f="genome/Rattus_norvegicus.Rnor_6.0.104.gtf" 
featurecounts_outdir="Read_Counts"
mkdir -p $featurecounts_outdir
featurecounts_outfile=$featurecounts_outdir"/featureCounts.counts"
echo "Outfile="$featurecounts_outfile

featureCounts -O -M -s 2 -T $threads_v -p -t gene -a $gtf_f -g gene_id -o $featurecounts_outfile $dirout_align/*.bam 2>&1 | tee $featurecounts_outdir'/featureCounts.log'

multiqc ./$featurecounts_outdir -o ./$featurecounts_outdir -n "Assignment_FeatureCounts_QC"


