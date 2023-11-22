#!/bin/sh
#$ -S /bin/sh
#$ -t 1-250
#$ -tc 250
#$ -l s_vmem=16G -l mem_req=16G
#$ -cwd
#$ -o out
#$ -e error

chr_num=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
p_name=(DNA_bending_stiffness DNA_denaturation Duplex_disrupt_energy Duplex_free_energy Protein_induced_deformability_Bp Stabilizing_energy_of_Z_DNA_AS Stabilizing_energy_of_Z_DNA_SA Stacking_energy Bruckner Packer)

length_chr_num=${#chr_num[*]}
length_p_name=${#p_name[*]}

i_length_chr_num=$((SGE_TASK_ID%length_chr_num))
i_length_p_name=$((SGE_TASK_ID/length_chr_num%length_p_name))

i_chr=${chr_num[$i_length_chr_num]}
i_p_name=${p_name[$i_length_p_name]}

$HOME/local/bin/R --vanilla --args $i_p_name chr$i_chr $HOME/tRNA/tss/GSE132660_SupplementaryData1.txt $HOME/DPP_hg38/${i_p_name}/chr${i_chr}.txt avgDPP <tRNA_avgDPP.R> error.log
