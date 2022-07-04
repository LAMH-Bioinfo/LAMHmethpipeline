#!/bin/bash
#SBATCH --job-name=novogene-job
#SBATCH --output=/shared/job_message/%x_%j.out
#SBATCH --ntasks-per-node=15


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
	echo "first parameter is the derectory of project";
	echo "second parameter is the name of samples"
	return 0;
}

if [ $# != 2 ];then 
	help_msg;
	exit;
fi

#name=LNT_8_5B
dir=$1
f_panel="/shared/colon_bed/colon_marker_V2_hg19_sorted_merged.bed"
name=$2

cd $1/$2"/bsmap_out/"
size_range="80_150"
alignmentSieve -b ${name}'_psorted_filtered_80_200.bam' --samFlagExclude 260 --ignoreDuplicates --maxFragmentLength 135 --minFragmentLength 65 -o ${name}'_psorted_filtered_'${size_range}'.bam' --filterMetrics metrics.txt
f_out=${dir}/${name}/bsmap_out/${name}"_CHALM_region_MBS_combined_length_"${size_range}".txt"
f_bam=${dir}/${name}/bsmap_out/${name}'_psorted_filtered_'${size_range}'.bam'
python /shared/python_script/CHALM_MBS/CHALM_methylation_block.py CHALM -d /shared/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
-x CG -R ${f_panel}   -L 120 -u -p -o ${f_out} ${f_bam}
awk '{if(NR >= 2) print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4}' ${name}'_CHALM_region_MBS_combined_length_'${size_range}'.txt' >  $name'_CHALM_region_MBS_combined_length_'${size_range}'.bed'
bedtools intersect -wa -wb -a ${f_panel}  -b $name'_CHALM_region_MBS_combined_length_'${size_range}'.bed' > $name'_CHALM_region_MBS_combined_length_'${size_range}'_selected.bed'
cp $name'_CHALM_region_MBS_combined_length_'${size_range}'_selected.bed'   ${dir}/all_MBS_deeptools_${size_range}

## 150_200bp 
size_range2="150_200"
alignmentSieve -b ${name}'_psorted_filtered_80_200.bam' --samFlagExclude 260 --ignoreDuplicates --maxFragmentLength 185 --minFragmentLength 135 -o ${name}'_psorted_filtered_'${size_range2}'.bam' --filterMetrics metrics.txt

f_out=${dir}/${name}/bsmap_out/${name}"_CHALM_region_MBS_combined_length_"${size_range2}".txt"
f_bam=${dir}/${name}/bsmap_out/${name}'_psorted_filtered_'${size_range2}'.bam'
python /shared/python_script/CHALM_MBS/CHALM_methylation_block.py CHALM -d /shared/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
-x CG -R ${f_panel}   -L 120 -u -p -o ${f_out} ${f_bam}

awk '{if(NR >= 2) print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4}' ${name}'_CHALM_region_MBS_combined_length_'${size_range2}'.txt' >  $name'_CHALM_region_MBS_combined_length_'${size_range2}'.bed'
bedtools intersect -wa -wb -a ${f_panel}  -b $name'_CHALM_region_MBS_combined_length_'${size_range2}'.bed' > $name'_CHALM_region_MBS_combined_length_'${size_range2}'_selected.bed'
cp $name'_CHALM_region_MBS_combined_length_'${size_range2}'_selected.bed'   ${dir}/all_MBS_deeptools_${size_range2}

date 
