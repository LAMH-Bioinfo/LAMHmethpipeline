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
    prefix=$2;
    cd ${dir};
    mkdir all_meth_deeptools_80_200 all_meth_deeptools_80_200_with_dup all_CHALM_deeptools_80_200  all_QC all_CHH job_message;
    for i in $(ls|grep ${prefix}); do echo $i; done > sample_id.txt;
    echo "done";
else
    help_message;
fi
