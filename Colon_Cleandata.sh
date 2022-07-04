#!/bin/bash


function help_message()
    {
     echo "the first parameter is the path of the project";
     echo "the second parameter is the prefix of the samples";
     return 0;
  }
##
if [ $# == 2 ] ; then
    dir=$1;
    sample=$2;
    for i in $(cat ${sample});
    do cd ${dir}/$i;
    rm  *_val_*;
    cd bsmap_out;
    rm $i".bam" $i"_mapped.bam" $i"_mapped_onTarget.bam" $i"_psorted.bam" $i"_psorted_filtered_80_200_with_dup.bam"  $i"_psorted_filtered_CHH.bam"; 
    done
else 
    help_message;
fi

echo "done"
