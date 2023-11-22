#!/bin/sh
#$ -S /bin/sh
#$ -t 1-25
#$ -tc 25
#$ -l s_vmem=16G -l mem_req=16G
#$ -cwd
#$ -o out
#$ -e error

chr_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

length_chr_num=${#chr_num[*]}

i_length_chr_num=$((SGE_TASK_ID%length_chr_num))

i_chr=${chr_num[$i_length_chr_num]}

$HOME/local/bin/R --vanilla --args chr$i_chr $HOME/tRNA/tss/GSE132660_SupplementaryData1.txt $HOME/hg38/chr${i_chr}.fa sequence <tRNA_sequence.R> error.log
