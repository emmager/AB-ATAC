#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_MACS2
#SBATCH --output=call_concensus_peaks
##SBATCH --error=ATACpipelineDev_MACS2.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=12:00:00
##SBATCH --array=207
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu
​
module load MACS2/2.2.7.1 

#Create Empty Files to bin samples by age
#touch /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/youngsamples.txt
#touch /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/oldersamples.txt


#SRR_file=$(sed -n ${SLURM_ARRAY_TASK_ID}p /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/PRJNA744248_metadata | cut -f 1 | sed '/^#/ ! s/$/.no-blacklist.bam/')
#Age=$(sed -n ${SLURM_ARRAY_TASK_ID}p /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/PRJNA744248_metadata | cut -f 2)
#Cell_Type=$(sed -n ${SLURM_ARRAY_TASK_ID}p /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/PRJNA744248_metadata | cut -f 3)


##Call Peaks for each individual sample##
#bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/ATACseqQC/*_shifted.bam)
#this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}
#prefix=`echo ${this_bam} | sed 's/_shifted.bam$//g'|sed 's/^\/Genomics\/ayroleslab2\/emma\/ATAC-Seq\/PRJNA744248\/ATACseqQC\///'`
#macs2 callpeak -t ${this_bam} -n /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/MACS2/${prefix} -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100


##Call Concensus Peaks for all samples pooled##
#echo "Call Concensus Peaks for all samples pooled"
macs2 callpeak -t /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/ATACseqQC/highqualsamples/*_shifted.bam -n /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/MACS2/highqualsamples/all_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100  --tempdir ./

##Call Concensus Peaks for all Young Samples##
#echo "Call Concensus Peaks for all Young Samples"
#youngbams=$(awk '$2 == "Young" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${youngbams} -n ./MACS2/young_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all Old Samples##
#echo "Call Concensus Peaks for all Old Samples"
#olderbams=$(awk '$2 == "Older" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${olderbams} -n ./MACS2/older_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD4 Bulk Naive Samples##
#echo "Call Concensus Peaks for all CD4 Bulk Naive Samples"
#CD4_BulkNaive=$(awk '$2 == "CD4_BulkNaive" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD4_BulkNaive} -n ./MACS2/CD4_BulkNaive_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD4 Bulk NonNaive Samples##
#echo "Call Concensus Peaks for all CD4 Bulk Naive Samples"
#CD4_BulkNonNaive=$(awk '$2 == "CD4_BulkNonNaive" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD4_BulkNonNaive} -n ./MACS2/CD4_BulkNonNaive_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD4 Tfh Samples##
#echo "Call Concensus Peaks for all CD4 Tfh Samples"
#CD4_Tfh=$(awk '$2 == "CD4_Tfh" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD4_Tfh} -n ./MACS2/CD4_Tfh_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD4 Treg Samples##
#echo "Call Concensus Peaks for all CD4 Treg Samples"
#CD4_Treg=$(awk '$2 == "CD4_Treg" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD4_Treg} -n ./MACS2/CD4_Treg_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 Bulk Naive Samples##
#echo "Call Concensus Peaks for all CD8 Bulk Naive Samples"
#CD8_BulkNaive=$(awk '$2 == "CD8_BulkNaive" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_BulkNaive} -n ./MACS2/CD8_BulkNaive_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 Bulk NonNaive Samples##
#echo "Call Concensus Peaks for all CD8 Bulk NonNaive Samples"
#CD8_BulkNonNaive=$(awk '$2 == "CD8_BulkNonNaive" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_BulkNonNaive} -n ./MACS2/CD8_BulkNonNaive_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 CM Samples##
#echo "Call Concensus Peaks for all CD8 CM Samples"
#CD8_CM=$(awk '$2 == "CD8_CM" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_CM} -n ./MACS2/CD8_CM_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 EM1 Samples##
#echo "Call Concensus Peaks for all CD8 EM1 Samples"
#CD8_EM1=$(awk '$2 == "CD8_EM1" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_EM1} -n ./MACS2/CD8_EM1_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 EM2 Samples##
#echo "Call Concensus Peaks for all CD8 EM2 Sample"
#CD8_EM2=$(awk '$2 == "CD8_EM2" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_EM2} -n ./MACS2/CD8_EM2_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 EMRA Samples##
#echo "Call Concensus Peaks for all CD8 EMRA Samples"
#CD8_EMRA=$(awk '$2 == "CD8_EMRA" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_EMRA} -n ./MACS2/CD8_EMRA_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 Naive Samples##
#echo "Call Concensus Peaks for all CD8 Naive Samples"
#CD8_Naive=$(awk '$2 == "CD8_Naive" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_Naive} -n ./MACS2/CD8_Naive_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 PD1+CD39+ Samples##
#echo "Call Concensus Peaks for all CD8 PD1+CD39+ Samples"
#CD8_PD1+CD39+=$(awk '$2 == "CD8_PD1+CD39+" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_PD1+CD39+} -n ./MACS2/CD8_PD1+CD39+_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 SCM-R3- Samples##
#echo "Call Concensus Peaks for all CD8 SCM-R3- Sample"
#CD8_SCM-R3-=$(awk '$2 == "CD8_SCM-R3-" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_SCM-R3-} -n ./MACS2/CD8_SCM-R3-_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

##Call Concensus Peaks for all CD8 SCM-R3+ Samples##
#echo "Call Concensus Peaks for all CD8 SCM-R3+ Samples"
#CD8_SCM-R3+=$(awk '$2 == "CD8_SCM-R3+" {print $1}' PRJNA744248_metadata | sed '/^#/ ! s/$/.no-blacklist.bam/')
#macs2 callpeak -t ${CD8_SCM-R3+} -n ./MACS2/CD8_SCM-R3+_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

