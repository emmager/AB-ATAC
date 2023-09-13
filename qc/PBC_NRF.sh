#!/bin/bash
#SBATCH --job-name=ATACpipelineDev_PBC_NRF
##SBATCH --output=ATACpipelineDev_PBC_NRF.%a.out
##SBATCH --error=ATACpipelineDev_PBC_NRF.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=6:00:00
##SBATCH --array=50
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load samtools
module load bedtools/2.29.1


#bams=(/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*markedDup.bam)

#this_bam=${bams[${SLURM_ARRAY_TASK_ID}]}

#prefix=`echo ${this_bam} | sed 's/.markedDup.bam$//g'`

#bamToBed -i ${this_bam} > ${prefix}.bed
​
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' SRR15054929.bed | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1}END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > SRR15054929_PBC.txt
​
echo -e 'TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair' | cat - SRR15054929_PBC.txt > SRR15054929.er1 && mv SRR15054929.er1 SRR15054929_PBC.txt


#Make Summary Files for Each of the Statistics

touch NRF.txt
SEEDS='/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*_PBC.txt'

for SEED in ${SEEDS}
do
awk 'NR%2==0' ${SEED} | awk '{print $5}' >> NRF.txt
done
paste SRR_Acc_List.txt NRF.txt > all_NRFs.txt

touch PBC1.txt
SEEDS='/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*_PBC.txt'

for SEED in ${SEEDS}
do
prefix=`echo ${SEED} | sed 's/_PBC.txt$//g'`
awk 'NR%2==0' ${prefix}_PBC.txt | awk '{print $6}' >> PBC1.txt
done
paste SRR_Acc_List.txt PBC1.txt > all_PBC1s.txt


touch PBC2.txt
SEEDS='/Genomics/ayroleslab2/emma/ATAC-Seq/PRJNA744248/*_PBC.txt'

for SEED in ${SEEDS}
do
prefix=`echo ${SEED} | sed 's/_PBC.txt$//g'`
awk 'NR%2==0' ${prefix}_PBC.txt | awk '{print $7}' >> PBC2.txt
done
paste SRR_Acc_List.txt PBC2.txt > all_PBC2s.txt
