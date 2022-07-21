# Liver NGS Pipeline  
## prepare for NGS process   
```bash
bash /shared/script_github/Liver_LAMH/HCC_preNGS_process.sh  [abs path of derectory]  [prefix of samples]
```
经过第一步，可以产生流程所需要的文件夹与样本列表文件。
如果样本prefix比较多，可以手动生成sample_id.txt文件。  
## run analysis pipeline  
可以使用以下命令进行所有样本的分析  
```
for i in $(cat sample_id.txt);do sbatch /shared/script_github/Liver_LAMH/HCC_NGS_process.sh [abs path of derectory] $i;done
```
运行上述命令后，所有样本的处理结果会投递到Amazon集群上去。可以等待结果生成。  
## run summary script   
sbatch /shared/script_github/Liver_LAMH/HCC_postNGS_process.sh [abs path of derectory]   
既可以完成样本结果的统计


test sample
