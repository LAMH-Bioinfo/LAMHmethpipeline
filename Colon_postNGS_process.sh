#!/bin/bash
#SBATCH --job-name=novogene-job
#SBATCH --output=%x_%j.out
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

#name=LNT_8_5B
dir=$1
cd ${dir}
#multiqc all_QC



f_panel='/shared/colon_bed/colon_marker_V2_hg19_sorted_merged_lambda_hg19.bed'
f_panel_no_control='/shared/colon_bed/colon_marker_V2_hg19_sorted_merged.bed'

#v13='/shared/novogene-data/X202SC20060567-Z01-F042/raw_data/model_v1.3_maker_sorted.bed'
#cd all_CHH
array=($(cat ${dir}/sample_id.txt))

cd all_CHH
for i in ${array[@]}; do awk '{total+=$5*$6;all+=$6}END{print total/all}' $i'_meth_combined_dedup_CHH.txt'; done > ct_conversion.txt
mv ct_conversion.txt ../
cd ../
for i in ${array[@]}; do cd $i; grep 'aligned pairs' $i'_bsmap_log.txt' | awk -F'(' '{print $2}' | awk -F')' '{print $1}'; cd ../; done > mapping_rate.txt
for i in ${array[@]}; do cd $i; grep 'unique pairs' $i'_bsmap_log.txt' | awk -F'unique pairs:' '{print $2}' | awk -F'(' '{print $2}' | awk -F')' '{print $1}'; cd ../; done > unique_mapping_rate.txt
for i in ${array[@]}; do cd $i'/bsmap_out'; samtools view $i'_mapped_onTarget.bam' | wc -l; cd ../../; done > onTarget_reads.txt
for i in ${array[@]}; do cd $i'/bsmap_out'; samtools view $i'_mapped.bam' | wc -l; cd ../../; done > total_reads.txt

for i in ${array[@]};
	do 
#	cd $i'/bsmap_out/';
	samtools index $dir/$i'/bsmap_out/'$i'_psorted_filtered_80_200.bam';
       	samtools view -L ${f_panel_no_control}  $dir/$i'/bsmap_out'/$i'_psorted_filtered_80_200.bam'|wc -l  >> target_reads.txt;  
        samtools view $dir/$i'/bsmap_out/'$i'_psorted_filtered_80_200.bam' 'pUC19'|wc -l >> pUC19_reads.txt;
	samtools view $dir/$i'/bsmap_out/'$i'_psorted_filtered_80_200.bam' 'Lambda_NEB'|wc -l >> Lambda_reads.txt;


# with dup 
	samtools index $dir/$i'/bsmap_out/'$i'_psorted_filtered_80_200_with_dup.bam'; 
	samtools view -L ${f_panel_no_control}  $dir/$i'/bsmap_out/'$i'_psorted_filtered_80_200_with_dup.bam'|wc -l >> target_reads_with_dup.txt;
	samtools view $dir/$i'/bsmap_out/'$i'_psorted_filtered_80_200_with_dup.bam' 'pUC19'|wc -l >> pUC19_reads_with_dup.txt;
	samtools view $dir/$i'/bsmap_out/'$i'_psorted_filtered_80_200_with_dup.bam' 'Lambda_NEB'|wc -l  >> Lambda_reads_with_dup.txt;

done;
cd all_meth_deeptools_80_200
for i in ${array[@]}; do grep -v pUC19 $i'_meth_combined_length_80_200_selected.bed' | grep -v Lambda_NEB | awk '{total+=$8}END{print total/(NR+0.1)}'; done > coverage.txt
for i in ${array[@]}; do grep -v pUC19 $i'_meth_combined_length_80_200_selected.bed' | grep -v Lambda_NEB | awk '{total+=$7}END{print total/(NR+0.1)}'; done > mean_meth.txt
for i in ${array[@]}; do grep pUC19 $i'_meth_combined_length_80_200_selected.bed' | awk -F'\t' 'NR>=20 && NR <=153{a+=$7*$8; b+=$8}END{print a/(b+0.1)}'; done > pUC19_meth.txt
for i in ${array[@]}; do grep pUC19 $i'_meth_combined_length_80_200_selected.bed' | awk -F'\t' '{total+=$8}END{print total/(NR+0.1)}'; done > pUC19_coverage.txt
for i in ${array[@]}; do grep Lambda_NEB $i'_meth_combined_length_80_200_selected.bed' | awk -F'\t' 'NR>=30 && NR <=3073{a+=$7*$8; b+=$8}END{print a/(b+0.1)}'; done > lambda_meth.txt
for i in ${array[@]}; do grep Lambda_NEB $i'_meth_combined_length_80_200_selected.bed' | awk -F'\t' '{total+=$8}END{print total/(NR+0.1)}'; done > lambda_coverage.txt


mv coverage.txt ../
mv mean_meth.txt ../
mv pUC19_meth.txt ../
mv lambda_meth.txt ../
mv pUC19_coverage.txt ../
mv lambda_coverage.txt ../

cd ../
cd all_meth_deeptools_80_200_with_dup
for i in ${array[@]}; do grep -v pUC19 $i'_meth_combined_length_80_200_selected_with_dup.bed' | grep -v Lambda_NEB | awk '{total+=$8}END{print total/NR}'; done > coverage_with_dup.txt
mv coverage_with_dup.txt ../
cd ../
paste -d "\t" sample_id.txt  mapping_rate.txt  unique_mapping_rate.txt  onTarget_reads.txt  total_reads.txt ct_conversion.txt coverage.txt  coverage_with_dup.txt mean_meth.txt target_reads.txt target_reads_with_dup.txt  Lambda_reads.txt Lambda_reads_with_dup.txt lambda_meth.txt lambda_coverage.txt  pUC19_reads.txt  pUC19_reads_with_dup.txt pUC19_meth.txt pUC19_coverage.txt > sample_QC.txt

sed -i '1s/^/sample_id\tmapping_rate\tunique_mapping_rate\ton_target_reads\ttotal_reads\tchh_ratio\tcoverage\tcoverage_with_dup\tmean_meth\ttarget_reads\ttarget_reads_with_dup\tlambda_reads\tlambda_reads_with_dup\tlambda_meth\tlambda_coverage\tpuc19_reads\tpuc19_reads_with_dup\tpuc19_meth\tpuc19_coverage\n/' sample_QC.txt
#awk -F'·\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"}' sample_QC.txt > sample_QC_cleaned.txt

