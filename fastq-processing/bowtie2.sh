#!/bin/bash

#SBATCH --job-name=ATACpipelineDev_bowtie2
#SBATCH --output=ATACpipelineDev_bowtie2.%a.out
#SBATCH --error=ATACpipelineDev_bowtie2.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --array=1-72%40
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

fastq_files=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA734635/*.trimmed.fastq)

this_fastq=${fastq_files[${SLURM_ARRAY_TASK_ID}]}

fastq1=`echo ${this_fastq}`

#fastq2=`echo ${this_fastq} | sed 's/-1/-4/g'`

prefix=`echo ${this_fastq} | sed 's/.trimmed.fastq//g'`

export BOWTIE2_INDEXES="/Genomics/ayroleslab2/emma/ATAC-Seq/mm9/"

/Genomics/argo/users/emmarg/lab/bin/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 -U ${fastq1} -x mm9 -k 4 --local -S ${prefix}.sam 2>$prefix.align.log

# -U option is for single end reads. Replace w/ -1 and -2 for paired end reads
#-2 ${fastq2} \
#-X 2000 \
