#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_MACS2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=15:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load MACS2/2.2.7.1 

#Create Empty Files to bin samples by age
#touch /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/youngsamples.txt
#touch /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/oldersamples.txt


#SRR_file=$(sed -n ${SLURM_ARRAY_TASK_ID}p /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/PRJNA744248_metadata | cut -f 1 | sed '/^#/ ! s/$/.no-blacklist.bam/')
#Age=$(sed -n ${SLURM_ARRAY_TASK_ID}p /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/PRJNA744248_metadata | cut -f 2)
#Cell_Type=$(sed -n ${SLURM_ARRAY_TASK_ID}p /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/PRJNA744248_metadata | cut -f 3)


#if [$Age = "Young"]; then
#       cp $SRR_file  ./YoungSamples
#else
#       cp $SRR_file ./OlderSamples
#fi

##Call Peaks for each individual sample##
#bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*.no-blacklist.bam)

#this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

#prefix=`echo ${this_bam} | sed 's/.no-blacklist.bam$//g'|sed 's/^\/Genomics\/ayroleslab2\/emma\/ATAC-Seq\/PRJNA744248\/testMACS\///'`

#macs2 callpeak -t /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/testMACS/* -n ./testMACS2/grouped -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100


##Call Concensus Peaks for all samples pooled##
echo "Call Concensus Peaks for all samples pooled"
macs2 callpeak -t /Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*.no-blacklist.bam -n ./MACS2/all_samples_concensus_peaks -f BAMPE -g 2.7e9 -q 0.5 --keep-dup all -B  --nomodel --extsize 200 --shift -100

