#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_finalBam
#SBATCH --output=ATACpipelineDev_fullpipeline.%a.out
#SBATCH --error=ATACpipelineDev_fullpipeline.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-1
#SBATCH --time=72:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

##THINGS TO DO BEFORE RUNNING THIS SCRIPT:
	#1) Make directory whose name is the GEO Project Number of the dataset you are trying to analyze
	#2) Download the SRR_Acc_List.txt file from GEO and ensure it is in the proper directory


module load sratoolkit/3.0.0
vdb-config --interactive
module load samtools
module load R
module load java
module load gatk/4.0.3.0

GEO_Project_Num=test

##FASTER-Q##
echo "---FASTERQ-DUMP---"
prefetch --option-file /Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/SRR_Acc_List.txt

for file in /Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/SSR*/
do
fasterq-dump $file
done

#fastqs1=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*1.fastq)
#fastqs2=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*2.fastq)

#this_fastq1=${fastqs1[${SLURM_ARRAY_TASK_ID}]}
#this_fastq2=${fastqs2[${SLURM_ARRAY_TASK_ID}]}

#mkdir /Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/QC_output

#fastqc $this_fastq1 -o /Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/QC_output/
#fastqc $this_fastq2 -o /Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/QC_output/
#wait 

##CUTADAPT##
#echo "---CUT ADAPT---"
#fastq_files=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*_1.fastq)

#this_fastq1=${fastq_files[${SLURM_ARRAY_TASK_ID}]}

#this_fastq2=`echo ${this_fastq1} | sed 's/_1/_2/g'`

#this_fastq1_out=`echo ${this_fastq1} | sed 's/fastq/trimmed.fastq/g' | sed 's@.*/@@'`

#this_fastq2_out=`echo ${this_fastq2} | sed 's/fastq/trimmed.fastq/g' | sed 's@.*/@@'`

#module load python/3.6.4
#pip3 install --user --upgrade cutadapt
	#~/.local/bin/cutadapt -m 5 -e 0.2 -a CTGTCTCTTATA -A CTGTCTCTTATA -o $this_fastq1_out -p $this_fastq2_out $this_fastq1 $this_fastq2
	#~/.local/bin/cutadapt -m 5 -e 0.2 -a CTTATACACATCT -A CTTATACACATCT -o $this_fastq1_out -p $this_fastq2_out $this_fastq1 $this_fastq2
#cutadapt -m 5 -e 0.2 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o $this_fastq1_out -p $this_fastq2_out $this_fastq1 $this_fastq2  #Illumina R1 + R2 Adapter trimming
#wait

##BOWTIE2##

#echo "---BOWTIE2---"

#trimmed_fastqs=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*_1.trimmed.fastq)

#this_fastq=${trimmed_fastqs[${SLURM_ARRAY_TASK_ID}]}

#fastq1=`echo ${this_fastq}`

#fastq2=`echo ${this_fastq} | sed 's/_1/_2/g'`

#prefix2=`echo ${this_fastq} | sed 's/_1.trimmed.fastq//g'`

#export BOWTIE2_INDEXES="/Genomics/ayroleslab2/emma/ATAC-Seq/hg19"

#/Genomics/ayroleslab2/emma/bin/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 \
#-x hg19 \
#-1 ${fastq1} \
#-2 ${fastq2} \
#-k 4 \
#-X 2000 \
#--local \
#-S ${prefix2}.sam \
#2>$prefix2.align.log
#wait

##SAMtoBAM##
#echo "---GET RAW BAM---"
#sams=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*sam)

#this_sam=${sams[${SLURM_ARRAY_TASK_ID}]}

#prefix3=`echo ${this_sam} | sed 's/.sam//g'`

#samtools view -h -b ${this_sam} | samtools sort -o $prefix3.raw.bam -
#wait

#samtools index ${prefix3}.raw.bam

#get stats for mtDNA and other alignment info
#samtools idxstats ${prefix3}.raw.bam > ${prefix3}.raw.bam.stats
#wait

##FILTER BAM###
#echo "---FILTER BAM---"
#bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*raw.bam)

#this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

#prefix4=`echo ${this_bam} | sed 's/.raw.bam$//g'`

# keep proper pair
# filter out mate unmapped
# filter out poor quality

#samtools view -bh -f 2 -F 524 -q 30 ${this_bam} > ${prefix4}.filt.bam

#sort by read name, fixmate for marking duplicates
#samtools sort -n ${prefix4}.filt.bam > ${prefix4}.filt.qnameSort.bam

#samtools fixmate ${prefix4}.filt.qnameSort.bam ${prefix4}.fixmate.bam

#samtools view -bh -f 2 -F 1804 ${prefix4}.fixmate.bam > ${prefix4}.fixmate.filt.bam

#samtools sort ${prefix4}.fixmate.filt.bam > ${prefix4}.toMark.bam

#rm ${prefix4}.fixmate.filt.bam

#rm ${prefix4}.fixmate.bam

#rm ${prefix4}.filt.qnameSort.bam

#rm ${prefix4}.filt.bam
#wait

##MARK DUPLICATES##
#echo "---MARK DUPLICATES---"
#filtbams=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*toMark.bam)

#this_filtbam=${filtbams[${SLURM_ARRAY_TASK_ID}]}

#prefix5=`echo ${this_filtbam} | sed 's/.toMark.bam$//g'`

#in_picard=/Genomics/grid/users/alea/programs/picard-tools-1.141/picard.jar

#in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0


#gatk  MarkDuplicates \
#-I ${this_filtbam} \
#-O ${prefix5}.markedDup.bam \
#-M ${prefix5}.markedDup.metrics.txt \
#--ASSUME_SORTED true \
#--VALIDATION_STRINGENCY LENIENT
#wait

##REMOVE DUPLICATES##
#"---REMOVE DUPLICATES---"
#markedbams=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*markedDup.bam)

#this_markedbam=${markedbams[${SLURM_ARRAY_TASK_ID}]}

#prefix6=`echo ${this_markedbam} | sed 's/.markedDup.bam$//g'`


# keep proper pair
# filter out duplicates
# filter out mate unmapped
# filter out non-primary alignments

#samtools view -bh -f 2 -F 1804 ${this_markedbam} > ${prefix6}.filt.dedup.bam

#samtools index ${prefix6}.filt.dedup.bam

# final clean alignment stats

#samtools idxstats ${prefix6}.filt.dedup.bam > ${prefix6}.stats
#wait

##FINAL BAM##
#echo "---GET FINAL BAM---"
#uniqbams=(/Genomics/ayroleslab2/emma/ATAC-Seq/$GEO_Project_Num/*dedup.bam)

#this_uniqbam=${uniqbams[${SLURM_ARRAY_TASK_ID}]}

#prefix7=`echo ${this_uniqbam} | sed 's/.dedup.bam$//g'`

#CHROMOSOMES=$(samtools view -H ${this_uniqbam} | grep '^@SQ' | cut -f2 | grep -v -e Y -e X -e M -e gl -e _random | sed 's/SN://' | xargs echo)

#CHROM=$(samtools view -H ${this_bam} | grep '^@SQ' | cut -f2 | grep -v -e Y -e X -e M -e gl -e _random | sed 's/SN:chr//' | xargs echo)

#prefix8=`echo ${this_uniqbam} | sed 's/trimmed.fastq.filt.dedup.bam$//'`

## filter out all chrX/Y, mtDNA, and scaffold alignments

#samtools view -bh -o ${prefix8}.final.bam ${this_uniqbam} $CHROMOSOMES

#samtools index ${prefix8}.final.bam

#samtools idxstats ${prefix8}.final.bam > ${prefix8}.final.bam.idxstats
