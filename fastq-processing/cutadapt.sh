#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_fastqc
#SBATCH --output=ATACpipelineDev_cutadapt-%A.%a.out
#SBATCH --error=ATACpipelineDev_cutadapt-%A.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00:00
#SBATCH --array=100-122
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

fastq_files=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA734635/*.fastq)

this_fastq1=${fastq_files[${SLURM_ARRAY_TASK_ID}]}

#this_fastq1=SRR15054980_1.fastq

#this_fastq2=`echo ${this_fastq1} | sed 's/_1/_2/g'`

this_fastq1_out=`echo ${this_fastq1} | sed 's/fastq/trimmed.fastq/g' | sed 's@.*/@@'`

#this_fastq2_out=`echo ${this_fastq2} | sed 's/fastq/trimmed.fastq/g' | sed 's@.*/@@'`

module load python/3.6.4
pip3 install --user --upgrade cutadapt
#~/.local/bin/cutadapt -m 5 -e 0.2 -a CTGTCTCTTATA -A CTGTCTCTTATA -o $this_fastq1_out -p $this_fastq2_out $this_fastq1 $this_fastq2
#~/.local/bin/cutadapt -m 5 -e 0.2 -a CTTATACACATCT -A CTTATACACATCT -o $this_fastq1_out -p $this_fastq2_out $this_fastq1 $this_fastq2
#cutadapt -m 5 -e 0.2 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o $this_fastq1_out -p $this_fastq2_out $this_fastq1 $this_fastq2  #Illumina R1 + R2 Adapter trimming
cutadapt -m 5 -e 0.2 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o $this_fastq1_out $this_fastq1
#pip3 uninstall cutadapt


