#/usr/bin/env python3 
#-*-conding:utf-8-*-

import os
import sys
import click
import argparse
import re
import json

class SampleInfo(object):
    '''
    collect the sample information of the sample map information
    '''

    def __init__(self,samplename ,path):
        '''
        :param samplename: name of the sample
        :param path: result file path of the sample
        '''
        self.name = samplename
        self.path = path
        self.dict = {}


    def getmapinfo(self):
        mapfile = self.path+"/"+self.name+"_R1_bismark_bt2_PE_report.txt"
        with open(mapfile ,'r') as f1:
            for line in f1:
                if re.search("Sequence pairs analysed in total",line):
                    self.dict['Total Reads'] = int(line.strip().split(":\t")[1])*2
                elif re.search("Number of paired-end",line):
                    self.dict['TOTAL_MAPPED_READS'] = int(line.strip().split(":\t")[1])*2
                elif re.search("Mapping efficiency",line):
                    self.dict['Mapping Rate'] = line.strip().split(":\t")[1]
                elif re.search("Total number of C's analysed",line):
                    self.dict["Total number of C's analysed"] = line.strip().split(":\t")[1]
                elif re.search("Total methylated C's in",line) or re.search("Total unmethylated C's in",line):
                    keys ,value = line.strip().split(":\t")[0],line.strip().split(":\t")[1]
                    self.dict[keys] = value
                elif re.search("C methylated in" ,line):
                    keys, value = line.strip().split(":\t")[0], line.strip().split(":\t")[1]
                    self.dict[keys] = value
                else:
                    continue
        mapfile2 = self.path+"/"+self.name+"_mapstat.txt"
        with open(mapfile2,'r') as f2:
            for line2 in f2:
                if re.search("insert size average",line2):
                    self.dict['Average_insert'] = line2.strip().split("\t")[2]


    def getdupinfo(self):
        dupfile =self.path+"/"+self.name+"_sort.deduplication_report.txt"
        with open(dupfile,'r') as dupf:
            for line in dupf:
                if re.search("Total number duplicated alignments",line):
                    info = line.strip().split(":\t")
                    remove_read = int(info[1].split("(")[0])*2
                    duprate     = info[1].split("(")[1].strip(")")
                    self.dict["Overall_duplicated_rate"] = duprate
                    self.dict["Dedup_Reads"] = remove_read
                else:
                    continue
        return self.dict



    def getmethyinfo(self):
        methyfile = self.path+"/"+self.name+"_sort_stat.txt"
        info = os.popen("cat %s |head -8 |tail -2"%methyfile).read()
        templist = [x.split("\t") for x in info.split("\n")]
        picard_dict = dict(zip(templist[0],templist[1]))
#        print(picard_dict)
        extractlist = ["ON_TARGET_BASES","PCT_SELECTED_BASES","MEAN_TARGET_COVERAGE","MEDIAN_TARGET_COVERAGE","PCT_TARGET_BASES_1X",
                       "PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X",
                       "PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X","FOLD_80_BASE_PENALTY","FOLD_ENRICHMENT"]
        for key,value in picard_dict.items():
            if key in extractlist:
                self.dict[key] = value
            else:
                continue

    def getcytotate(self):
        lamda = self.path +"/" + self.name+"_lambda_rate.txt"
        chrm  = self.path +"/" + self.name+"_chrM_rate.txt"
        with open(lamda ,'r') as f2:
            for line in f2:
                if re.search("Number of cytosines",line):
#                    print(line.strip().split(": "))
                    self.dict['LAMBDA_Cs'] = line.strip().split(": ")[1]
                elif re.search("Bisulfite conversion rate",line):
                    self.dict["LAMBDA_CONVERSION_EFFICIENCY"] = line.strip().split(": ")[1]
                else:
                    continue 
        
        with open(chrm ,'r') as f3:
            for line in f3:
                if re.search("Number of cytosines",line):
                    self.dict['CHRM_Cs'] = line.strip().split(": ")[1]
                elif re.search("Bisulfite conversion rate",line):
                    self.dict["CHRM_CONVERSION_EFFICIENCY"] = line.strip().split(": ")[1]
                else:
                    continue 



    def getqcinfo(self):
        """
        get qc info from the json file
        """
        qcfile = self.path +"/../qc/"+self.name +"_qc.json"
        with open(qcfile,'r') as f1:
            qcinfo = json.load(f1)
            self.dict.update({"before_"+x:y for x,y in qcinfo["summary"]["before_filtering"].items()})
            self.dict.update({"after_"+x:y for x,y in qcinfo["summary"]["after_filtering"].items()})
            self.dict["Pass_filter_rate"] = self.dict['after_total_reads']/self.dict['before_total_reads']
            self.dict.update(qcinfo['filtering_result'])

        return self.dict
     
parser = argparse.ArgumentParser(description='sample name and result path')
parser.add_argument("--samples" ,'-s',help="the sample information to collect",required=True)
parser.add_argument("--path",'-p' ,help = "the path of the project ",required=True)
argv = vars(parser.parse_args())


path    = argv['path']
samples = argv['samples']

sinfo = SampleInfo(samples,path)
#print(testdd.name,testdd.path,testdd.dict)
#sinfo.getdupinfo()
#print(testdd.dict,len(testdd.dict))
#print(dupfile)
#print(testdd.dict)
sinfo.getqcinfo()
sinfo.getmapinfo()
sinfo.getdupinfo()
#print(testdd.dict,len(testdd.dict))
sinfo.getmethyinfo()
#print(testdd.dict,len(testdd.dict))
sinfo.getcytotate()

with open(path + "/"+samples+"_summary.txt",'w') as  f1:
    f1.write("Item" +  "\t" +  samples + "\n")
    for item,value in sinfo.dict.items():
        if re.search("Unknown context",item) or re.search("unknown context",item):
            continue 
        f1.write(item + "\t" + str(value) + "\n")
