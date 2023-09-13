#!/bin/bash
#SBATCH --job-name=R_test
#SBATCH --output=R_test.%a.out
#SBATCH --error=R_test.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=06:00:00
#SBATCH --array=207
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load R/4.1.1

Rscript /Genomics/argo/users/emmarg/lab/ATAC-Seq/scripts/ATACseqQC.R ${SLURM_ARRAY_TASK_ID}

#sed -i 's/","/\t/g' /Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/*_bamQC.txt
#sed -i 's/,/\t/g' /Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/*_bamQC.txt
#sed -i 's/^"	//' /Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/*_bamQC.txt 
#sed -i '2s/^"1"//' /Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/*_bamQC.txt
#sed -i 's/"$//' /Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/*_bamQC.txt
#sed -i '2s/^	//' /Genomics/argo/users/emmarg/lab/ATAC-Seq/PRJNA744248/ATACseqQC/*_bamQC.txt


