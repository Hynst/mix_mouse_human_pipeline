#!/bin/bash

### 1) create reference assembly where only archer sequences occur for GRCh37 and GRm38 ###

REF_HUMN="/mnt/ssd/ssd_3/references/homsap/GRCh37-p13/seq/GRCh37-p13.fa"
REF_MOUS="/mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/ref/GRCm38.p6-93_MOUS_contigs.fa"

# use bedtools getfasta (ref vs. BED)
#BED_ARCHER_HUMN="/mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/bed/jana_archer_unique_plus2nt.bed"
#BED_ARCHER_MOUS="/mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/bed/jana_archer_unique_plus2nt_MYS.bed"

# human
#echo "Human"
#bedtools getfasta -fi $REF_HUMN -bed $BED_ARCHER_HUMN -fo /mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/ref/GRCh37-p13_Archer.fa

# mouse
#echo "Mouse"
#bedtools getfasta -fi $REF_MOUS -bed $BED_ARCHER_MOUS -fo /mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/ref/GRCm38.p6-93_Archer.fa

# mix
cat $REF_HUMN \
> /mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/ref/MIX_GRCh37-p13_GRCm38.p6-93.fa
cat $REF_MOUS \
>> /mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/ref/MIX_GRCh37-p13_GRCm38.p6-93.fa

# create BWA index for MIX reference
/mnt/ssd/ssd_2/install/dir/anaconda/envs/lymphopanel/bin/bwa index -a bwtsw \
	/mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/ref/MIX_GRCh37-p13_GRCm38.p6-93.fa \
	-p /mnt/ssd/ssd_2/bioda_temp/kuba/MartinCul_mysi/ref/BWA/MIX_GRCh37-p13_GRCm38.p6-93
