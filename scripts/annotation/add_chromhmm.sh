#!/bin/bash

cd ~/Desktop/SSP/scripts/annotation/bedtools2
# cd ~/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-01/data/comprehensive/ATAC_comprehensive.bytype

module load ./

# PBMC
outdir=~/Desktop/SSP/DA_7_12_17
data=~/Desktop/SSP/DA_7_12_17/monot_merged_consensus_raw_atac_whitelisted.txt
hmmanno_segments=~/Desktop/SSP/Data/support/E062_18_core_K27ac_segments.bed

echo "now processing $data"
# touch ${outdir}/$(basename $data .txt).bed
bed=${outdir}/$(basename $data .txt).bed
perl -ne 'print if $.>1' $data | cut -f1-3 > $bed
fullanno=${outdir}/$(basename $data .txt)_fullanno.txt
anno=$(basename $data .txt)_annotated

head -1 $data | cut -f1-3 | sed 's/$/	chromHMMstate/g' > h.txt
source ~/.bash_profile
bedtools intersect -wo -a $bed -b $hmmanno_segments | cut -f1-3,7 | sortBed | uniq > $(basename $data .txt)_chromHMM.tmp
sed s/E18/Quies/g $(basename $data .txt)_chromHMM.tmp \
	| sed s/E17/ReprPCWk/g | sed s/E16/ReprPC/g | sed s/E15/EnhBiv/g | sed s/E14/TssBiv/g | sed s/E13/Het/g | sed s/E12/ZNF_Rpts/g \
	| sed s/E11/EnhWk/g | sed s/E10/EnhA2/g | sed s/E9/EnhA1/g | sed s/E8/EnhG2/g | sed s/E7/EnhG1/g | sed s/E6/TxWk/g | sed s/E5/Tx/g \
	| sed s/E4/TssFlnkD/g | sed s/E3/TssFlnkU/g | sed s/E2/TssFlnk/g | sed s/E1/TssA/g | cat h.txt - \
	> ${outdir}/$(basename $data .txt)_chromHMM.txt


annotatePeaks.pl $bed hg19 > $fullanno
cut -f2-4,8,10,16,19 $fullanno | sed 's/ //g' | sed s/\(.*\)//g | sed 's/\.[1-9]	/	/g' | sed 's/\.[1-9][0-9]	/	/g' > ${anno}.tmp
head -1 ${anno}.tmp | sed 's/Chr/chr/g' | sed 's/Start/start/g' | sed 's/End/end/g' > h.txt
perl -ne 'print if $.>1' ${anno}.tmp | sortBed | perl -lane 'print $F[0],"\t",$F[1]-1,"\t",$F[2],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[6]' | cat h.txt - > ${outdir}/${anno}.txt
rm *tmp
rm h.txt