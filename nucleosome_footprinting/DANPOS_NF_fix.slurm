#!/bin/bash
#SBATCH --job-name=novogene-job
#SBATCH --output=/shared/MingmaData/Colon_case_control_574samples/nucleosome_footprinting_80_200/job_message/%x_%j.out
#SBATCH --ntasks-per-node=12


### JOB PREP ###

### check $HOSTNAME of host running job
#echo "hostname:  $HOSTNAME"

### check $HOME of user submitting job
#echo "pbs_o_home:  $PBS_O_HOME"

### check local scratch on compute node running job
#echo "tmpdir:  $TMPDIR"

echo "host:  $PBS_O_HOST"
### check $CWD at job submission
echo "pbs_o_workdir:  $PBS_O_WORKDIR"
date

source activate py37

#name=LNT_8_5B
dir="/shared/MingmaData/Colon_case_control_574samples/nucleosome_footprinting_80_200/"${name}
regionFile="/shared/colon_bed/colon_PAS_promoter_hg19.bed6"
hg19ChrSize="/shared/colon_bed/hg19.chrom.sizes"
f_rscript="/shared/software/NF_scripts_yumei/scripts/columnMeanByFactor.R"

cd ${dir}
size_range="80_200"
f_bam=${name}'_psorted_filtered_'${size_range}'.bam'
f_cleaned_bam=${name}'_filtered_'${size_range}'_psorted_cleaned.bam'
f_sorted_bam=${name}'_filtered_'${size_range}'_psorted.bam'


samtools sort -O bam -o ${f_sorted_bam} ${f_bam}
samtools index ${f_sorted_bam}
samtools view -F 268 -b ${f_sorted_bam} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 > ${f_cleaned_bam}
samtools index ${f_cleaned_bam}

inFile=${f_cleaned_bam}
python /shared/software/DANPOS3-master/danpos.py dpos ${f_cleaned_bam} --paired 1 -u 0 -c 1000000 -o Rst_danpos2
outPrefix=$(echo $inFile|sed 's/.bam//');
prefix=${name}
wigToBigWig -clip Rst_danpos2/pooled/${outPrefix}*.wig $hg19ChrSize Rst_danpos2/pooled/${outPrefix}.bw
## occupancy
bigWigAverageOverBed Rst_danpos2/pooled/${outPrefix}.bw $regionFile ${prefix}.out.occupancy.tsv
## fuzziness
mkdir intersectRst meanByMarker
awk -v OFS="\t" 'NR>1{print $1,$4-1,$4,$5,$6}' Rst_danpos2/pooled/${outprefix}*.positions.xls|grep -vE "Lambda_NEB|pUC19"|bedtools intersect -a $regionFile -b stdin -wo >intersectRst/${prefix}.intersect.tsv
cut -f4,11  intersectRst/${prefix}.intersect.tsv | Rscript ${f_rscript} -c=2 -f=1 -o=meanByMarker/${prefix}.outMean.fuzziness.tsv;


date
