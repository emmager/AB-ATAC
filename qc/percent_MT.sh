#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_finalreads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

touch percent_MT.txt

SEEDS='/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*.stats'

for SEED in ${SEEDS}
do
prefix=`echo ${SEED} | sed 's/.stats//'`

echo ${prefix} >> percent_MT.txt
done


touch temp.txt

for SEED in ${SEEDS}
do
cut -f3 ${SEED} | paste -s -d+ - | bc > mapped
awk -v reads=`cat mapped` '$1 == "chrM" {print $3/reads}' ${SEED} >> temp.txt
done

paste percent_MT.txt temp.txt > tmp
#awk '{printf "%s%s",$0,NR%3?"\t":RS}' tmp > temp.txt
mv tmp percent_MT.txt
rm tmp

