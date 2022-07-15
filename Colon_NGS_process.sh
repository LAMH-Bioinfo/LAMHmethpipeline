#!/bin/bash
#SBATCH --job-name=helio_methy
#SBATCH --output=job_message/%x_%j.out
#SBATCH --ntasks-per-node=9


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

function help_msg(){
	echo "first parameter The path of the project";
	echo "the second parameter is the name of sample";
	return 0;
}


if [ $# != 2 ]; then 
    help_msg;
    exit;
fi

## define the reference file
fa="/shared/reference/hg19_lambda/hg19_Lambda_pUC19.fa"
f_panel='/shared/colon_bed/colon_marker_V2_hg19_sorted_merged_lambda_hg19.bed'

## input parameter
dir=$1
name=$2
cd ${dir}/${name}
size_range='80_200'

fastqc *fastq.gz

mkdir trim_galore_out
sample_list_1=($(ls -1 | grep '_R1.fastq.gz'))
sample_list_2=($(ls -1 | grep '_R2.fastq.gz'))

for((i=0;i<${#sample_list_1[@]};i++))
do
trim_galore --paired  -q 20  --clip_R1 5 --clip_R2 10 --three_prime_clip_R1 30 --three_prime_clip_R2 30 -o trim_galore_out/ ${sample_list_1[i]}  ${sample_list_2[i]}
done


cd trim_galore_out
mv *val_1.fq.gz ../
mv *val_2.fq.gz ../


cd ../
### bsmap to mapping reads

sample_list_1=($(ls -1 | grep '_R1_val_1.fq.gz'))
sample_list_2=($(ls -1 | grep '_R2_val_2.fq.gz'))

if [ ${#sample_list_2[@]} -gt 0 ]
then
for((i=0;i<${#sample_list_1[@]};i++))
do
bsmap -a ${sample_list_1[i]} -b ${sample_list_2[i]} -d ${fa} -S 100 -R -p 1 -o $(echo ${sample_list_1[i]} | cut -d'.' -f1)'.bam' &> ${name}'_bsmap_log.txt'
done
else
for((i=0;i<${#sample_list_1[@]};i++))
do
bsmap -a ${sample_list_1[i]} -d ${fa}  -R -p 16 -o $(echo ${sample_list_1[i]} | cut -d'.' -f1)'.bam' &> ${name}'_bsmap_log.txt'
done
fi

### samtools to merge the output bamfiles    
mkdir bsmap_out
if [ ${#sample_list_1[@]} -gt 1 ]
then
samtools merge -ru 'bsmap_out/'$name'.bam' *.bam
else
mv $(echo ${sample_list_1[0]} | cut -d'.' -f1)'.bam' 'bsmap_out/'$name'.bam'
fi

cd bsmap_out
## change -F 260 to 268 
samtools view  -F 268 -b  ${name}'.bam' > ${name}'_mapped.bam'
samtools view -b -L ${f_panel} ${name}'_mapped.bam'  > ${name}'_mapped_onTarget.bam'
samtools sort ${name}'_mapped_onTarget.bam' ${name}'_psorted'
samtools index ${name}'_psorted.bam'
/shared/software/anaconda/bin/bamPEFragmentSize  --histogram ${name}'_fragmentSize.png' -T 'Fragment size' --maxFragmentLength 500 -b ${name}'_psorted.bam' --samplesLabel ${name} --outRawFragmentLengths ${name}'_fragment_length.tsv'
/shared/software/Miniconda/bin/alignmentSieve -b ${name}'_psorted.bam' --samFlagExclude 268 --ignoreDuplicates --maxFragmentLength 185 --minFragmentLength 65 -o ${name}'_psorted_filtered_'${size_range}'.bam' --filterMetrics metrics.txt
/shared/software/Miniconda/bin/alignmentSieve -b ${name}'_psorted.bam' --samFlagExclude 268 --maxFragmentLength 185 --minFragmentLength 65 -o ${name}'_psorted_filtered_'${size_range}'_with_dup.bam' --filterMetrics metrics.txt
## for CHH calculation
/shared/software/Miniconda/bin/alignmentSieve -b ${name}'_psorted.bam' --samFlagExclude 268 --ignoreDuplicates -o ${name}'_psorted_filtered_CHH.bam' --filterMetrics metrics.txt

f_out=${dir}/${name}/"/bsmap_out/"${name}"_meth_combined_length_"${size_range}".txt"
f_bam=${dir}/${name}"/bsmap_out/"${name}'_psorted_filtered_'${size_range}'.bam'

python /shared/software/bsmap-2.90/bin/methratio.py -d ${fa} -x CG -i no-action -u -g -p \
-o ${f_out} ${f_bam}

f_out=${dir}/${name}/"/bsmap_out/"${name}"_meth_combined_length_"${size_range}"_with_dup.txt"
f_bam=${dir}/${name}/"/bsmap_out/"${name}'_psorted_filtered_'${size_range}'_with_dup.bam'
python /shared/software/bsmap-2.90/bin/methratio.py -d ${fa} -x CG -i no-action -u -g -p \
-o ${f_out} ${f_bam}

##
f_out=${dir}/${name}"/bsmap_out/"${name}"_CHALM_combined_length_"${size_range}".txt"
f_bam=${dir}/${name}"/bsmap_out/"${name}'_psorted_filtered_'${size_range}'.bam'
python /shared/python_script/CHALM/src/CHALM.py trad -d ${fa} \
-x CG -i no-action -u -g -p -l 1 -o ${f_out} ${f_bam}

##
f_out=${dir}/${name}/"/bsmap_out/"${name}"_meth_combined_dedup_CHH.txt"
f_bam=${dir}/${name}/"/bsmap_out/"${name}"_psorted_filtered_CHH.bam"
python /shared/software/bsmap-2.90/bin/methratio.py -d ${fa} -x CHH -i no-action -u -g -p \
-o ${f_out} ${f_bam}

## select CpG sites in the target regions
awk '{if(NR >= 2) print $1"\t"$2"\t"($2+1)"\t"$5"\t"$6"\t"$3}' $name'_meth_combined_length_'${size_range}'.txt' >  $name'_meth_combined_length_'${size_range}'.bed'
bedtools intersect -wa -wb -a ${f_panel} -b $name'_meth_combined_length_'${size_range}'.bed' > $name'_meth_combined_length_'${size_range}'_selected.bed'

awk '{if(NR >= 2) print $1"\t"$2"\t"($2+1)"\t"$5"\t"$6"\t"$3}' $name'_meth_combined_length_'${size_range}'_with_dup.txt' >  $name'_meth_combined_length_'${size_range}'_with_dup.bed'
bedtools intersect -wa -wb -a ${f_panel} -b $name'_meth_combined_length_'${size_range}'_with_dup.bed' > $name'_meth_combined_length_'${size_range}'_selected_with_dup.bed'

awk '{if(NR >= 2) print $1"\t"$2"\t"($2+1)"\t"$5"\t"$6"\t"$3}' $name'_CHALM_combined_length_'${size_range}'.txt' >  $name'_CHALM_combined_length_'${size_range}'.bed'
bedtools intersect -wa -wb -a ${f_panel} -b $name'_CHALM_combined_length_'${size_range}'.bed' > $name'_CHALM_combined_length_'${size_range}'_selected.bed'


cp $name'_meth_combined_length_'${size_range}'_selected.bed' ${dir}'/all_meth_deeptools_'${size_range}
cp $name'_meth_combined_length_'${size_range}'_selected_with_dup.bed' ${dir}'/all_meth_deeptools_'${size_range}'_with_dup'
cp $name'_CHALM_combined_length_'${size_range}'_selected.bed' ${dir}'/all_CHALM_deeptools_'${size_range}
cp ${name}"_meth_combined_dedup_CHH.txt" ${dir}'/all_CHH'

cd ../
cp *_fastqc.zip ${dir}/all_QC
cd trim_galore_out 
cp *_trimming_report.txt ${dir}/all_QC



date
