0#!/bin/bash

#PBS -w /home/s-sengua/scripts/
#PBS -l nodes=1:ppn=16,walltime=48:00:00 
#PBS -m a
#PBS -t 12-86
#PBS -M Anooshree.Sengupta@jax.org 
#PBS -e /home/s-sengua/scripts/fetchperr
#PBS -o /home/s-sengua/scripts/fetchpout

module load bowtie2/2.2.3
module load bwa/0.7.12
module load samtools/0.1.19
module load fastx/0.14
module load Trimmomatic/0.33
module load Python/2.7.3
module load R
module load perl/5.10.1
module load bedtools/2.26.0

sra="SRR51988${PBS_ARRAYID}"

# ~/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump -Z -X 5 /sdata/dbGaP/13995/SRR5198832.sra

# sra="SRR390728"

## 0) Pre-fetch SRA files from NCBI server, move to work folder
## NOT NECESSARY SINCE SRA FILES ALREADY IN HELIX
# prefetch $sra
# mv ~/ncbi/public/sra/${sra}.sra /sdata/dbGaP/13995

## 1) run sra-tools to convert sra files into fastq files
## Q = 33 - stick to the sanger offset - latest version of illumina uses this too!
sraName="/sdata/dbGaP/13995/${sra}.sra"
fileName=$(basename ${sraName}.sra)
filePath="/sdata/dbGaP/13995/FASTQ"
[[ ! -d "${filePath}" ]] && mkdir "${filePath}";
cd ~/ncbi/dbGaP-13995/
~/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --outdir $filePath --split-files $sraName

## 2) fastx quality filter to remove sequences that have phred <30 on 70% seq
ffilePath="/sdata/dbGaP/13995/FASTQ_Filter"
[[ ! -d "${ffilePath}" ]] && mkdir "${ffilePath}";
outname="${ffilePath}/${fileName}_1_filtered.fastq"
fastq_quality_filter -v -q 30 -p 70 -i ${filePath}/${sra}_1.fastq -o $outname > ${outname}.out

## 3) fastQC quality metrics
qfilePath="/sdata/dbGaP/13995/FASTQC"
[[ ! -d "${qfilePath}" ]] && mkdir "${qfilePath}";
/opt/compsci/FastQC/0.11.3/fastqc -t 2 --noextract ${outname} -o ${qfilePath}

## 4) trimmomatic read trimming + adapter trimming, re-run fastQC
trfilePath="/sdata/dbGaP/13995/trimmomatic"
[[ ! -d "${trfilePath}" ]] && mkdir "${trfilePath}";
troutname=${trfilePath}/$(basename "${outname}" .fastq)_trimmed.fastq
java -jar /opt/compsci/Trimmomatic/0.33/trimmomatic-0.33.jar SE -threads 2 ${outname} ${troutname} ILLUMINACLIP:/sdata/dbGaP/13995/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
q2filePath="/sdata/dbGaP/13995/FASTQC_trimmed"
[[ ! -d "${q2filePath}" ]] && mkdir "${q2filePath}";
/opt/compsci/FastQC/0.11.3/fastqc -t 2 --noextract ${troutname} -o ${q2filePath}

## 5) bwa alignment (hg19)
mfilePath="/sdata/dbGaP/13995/bwa"
[[ ! -d "${mfilePath}" ]] && mkdir "${mfilePath}";
fileNameSAM=${mfilePath}/$(basename "${troutname}" .fastq).sam
bwa mem -M /data/mby/INDEXES/HUMAN/BWA/hg19.fa ${troutname} > ${fileNameSAM}

## 6) shifting, cleaning dupes, and converting to bam
sfilePath="/sdata/dbGaP/13995/SAM"
[[ ! -d "${sfilePath}" ]] && mkdir "${sfilePath}";
fileSorted=${sfilePath}/$(basename "${fileNameSAM}" .sam)_sorted.sam
java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/SortSam.jar INPUT=${fileNameSAM} OUTPUT=${fileSorted} SO=coordinate
        
fileRmdup=${sfilePath}/$(basename "${fileSorted}" .sam)_rmdup.sam
fileMetrics=${sfilePath}/$(basename "${fileSorted}" .sam)_rmdup_metrics.txt
# java -Xms1g -Xmx4g -jar /opt/compsci/picard/1.95/MarkDuplicates.jar INPUT=${fileSorted} OUTPUT=${fileRmdup} METRICS_FILE=${fileMetrics} REMOVE_DUPLICATES=true READ_NAME_REGEX=null
        
## shift reads and convert to BAM
bfilePath="/sdata/dbGaP/13995/BAM"
[[ ! -d "${bfilePath}" ]] && mkdir "${bfilePath}";
fileBam=${bfilePath}/$(basename "${fileRmdup}" .sam).bam
samtools view -bS ${fileRmdup} > ${fileBam}

fileShift=${bfilePath}/$(basename "${fileRmdup}" .bam)_shifted
bedtools bamtobed -i ${fileBam} | awk -v OFS="\t" '{if($6=="+"){print $1,$2+4,$3,$4,$5,$6} else if($6=="-"){print $1,$2-5,$3,$4,$5,$6}}' | awk -v OFS="\t" '{if($2<0){print $1,0,$3,$4,$5,$6} else if($2>=0){print $1,$2,$3,$4,$5,$6}}' > ${fileShift}.bed
bedtools bedtobam -i ${fileShift}.bed -g /opt/compsci/BEDtools/2.26.0/genomes/human.hg19.genome > ${fileShift}.bam

# sort & index BAM, create final output BAM file
outputPath="/sdata/dbGaP/13995/output"
[[ ! -d "${outputPath}" ]] && mkdir "${outputPath}";
fileSort=${outputPath}/$(basename "${fileShift}" .bam)_sorted
samtools sort ${fileShift}.bam ${fileSort}
samtools index ${fileSort}.bam