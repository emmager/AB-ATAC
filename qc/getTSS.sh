#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_getTSS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

tail -n +2 $1 | awk '
BEGIN{OFS="\t"}{
if($4=="+") {print $3,$5,$5+1,$2 "_" $13,".",$4}
 else {print $3,$6-1,$6,$2 "_" $13,".",$4}
 }' > "$1".TSS.bed
sort -k1,1 -k2,2n -k3,3n -k6,6 -u "$1".TSS.bed > "$1".TSS.unique.bed
awk -v "OFS=\t" '(($2-1000) >=0){print $1,$2-1000,$3+1000,$4,$5,$6}' "$1".TSS.unique.bed | grep -v "chrM" > "$1".TSS.unique.2K.bed 
#/Genomics/ayroleslab2/emma/ATAC-Seq/ucsc_refseq_hg38_ncbi.txt
