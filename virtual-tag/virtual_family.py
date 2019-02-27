# coding: utf-8
import time
import pysam
import collections 
import sys
import re
import os
import fnmatch
import pandas as pd
import pickle
import itertools
import numpy as np
import argparse
import ast    

def get_family(samfile,chromosome,position,alt,switch1=True,switch2=True,minquality=30,mq=30,withtag=''):
    templen=''
    start=''
    tag=''
    base=collections.defaultdict(list)
    strand=collections.defaultdict(list)
    qname =collections.defaultdict(list)
    for pileupcolumn in samfile.pileup(chromosome, position-1, position,max_depth=100000,truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                if pileupread.alignment.query_qualities[pileupread.query_position]>=minquality:
                    qname[pileupread.alignment.qname].append(pileupread.alignment.query_sequence[pileupread.query_position])

    for pileupcolumn in samfile.pileup(chromosome, position-1, position,max_depth=100000,truncate=True):
        for pileupread in pileupcolumn.pileups:
            if withtag:
                templen=str(abs(pileupread.alignment.template_length))
                start=min(pileupread.alignment.reference_start,pileupread.alignment.next_reference_start)
                tag=pileupread.alignment.query_name.split('|')[0]
            else:
                start=min(pileupread.alignment.reference_start,pileupread.alignment.next_reference_start)
                templen=str(abs(pileupread.alignment.template_length))
                tag = ''
            if pileupread.alignment.is_reverse is pileupread.alignment.is_read1:
                base_strand='-'
            else:
                base_strand='+'
            if  not pileupread.is_del and not pileupread.is_refskip:
                if pileupread.alignment.query_qualities[pileupread.query_position]< minquality:
                    continue
                if pileupread.alignment.mapping_quality <mq:
                    continue
                if not pileupread.alignment.is_proper_pair:
                        continue
                if (switch1 and switch2):#switch1 is the length filter, switch2 is the overlapping filter
                    if pileupcolumn.pos <= (start+int(templen)-1):
                        if withtag:
                            if len(qname[pileupread.alignment.qname]) ==1:
                                base[str(start)+','+templen+','+tag].append(pileupread.alignment.query_sequence[pileupread.query_position])
                                strand[str(start)+','+templen+','+tag].append(base_strand)
                            elif len(qname[pileupread.alignment.qname]) >1:
                                if len(np.unique(qname[pileupread.alignment.qname]))>1:
                                    continue   
                                elif len(np.unique(qname[pileupread.alignment.qname])) ==1:
                                    base[str(start)+','+templen+','+tag].append(pileupread.alignment.query_sequence[pileupread.query_position])
                                    strand[str(start)+','+templen+','+tag].append(base_strand)
                                    del qname[pileupread.alignment.qname]
                                else:
                                    print "other information"
                                    continue
                        else:
                            if len(qname[pileupread.alignment.qname]) ==1:
                                base[str(start)+','+templen+','+tag].append(pileupread.alignment.query_sequence[pileupread.query_position])
                                strand[str(start)+','+templen+','+tag].append(base_strand)
                            elif len(qname[pileupread.alignment.qname]) >1:
                                if len(np.unique(qname[pileupread.alignment.qname]))>1:
                                    continue
                                elif len(np.unique(qname[pileupread.alignment.qname])) == 1:
                                    base[str(start)+','+templen+','+tag].append(pileupread.alignment.query_sequence[pileupread.query_position])
                                    strand[str(start)+','+templen+','+tag].append(base_strand)
                                    del qname[pileupread.alignment.qname]
                                else:
                                    continue
                elif switch1:#switch 1 open
                    if (pileupcolumn.pos <= start+int(templen)-1):
                        base[str(start)+','+templen+','+tag].append(pileupread.alignment.query_sequence[pileupread.query_position])
                        strand[str(start)+','+templen+','+tag].append(base_strand)
                elif switch2: #switch2 open
                    if len(qname[pileupread.alignment.qname]) ==1:
                        base[str(start)+','+templen+','+tag].append(pileupread.alignment.query_sequence[pileupread.query_position])
                        strand[str(start)+','+templen+','+tag].append(base_strand)
                    elif len(qname[pileupread.alignment.qname]) >1:
                        if len(np.unique(qname[pileupread.alignment.qname]))>1:
                            continue
                        elif len(np.unique(qname[pileupread.alignment.qname])) == 1:
                            base[str(start)+','+templen+','+tag].append(pileupread.alignment.query_sequence[pileupread.query_position])
                            strand[str(start)+','+templen+','+tag].append(base_strand)
                            del qname[pileupread.alignment.qname]
                        else:
                            continue
    return [base,strand] 
    
def get_var_info(family,alt,min_family_size=2,cutoff=1.0,SUMMARY=True):
    var_single_plus=0
    var_single_minus=0
    non_var_single_plus=0
    non_var_single_minus=0
    var_single=0
    var_raw=0
    depth_raw=0
    family_alt =0 
    family_has_alt=0 
    all_family_number = len(family[0])
    family_strand_distribution = collections.defaultdict(list)
    family_plus=0
    family_minus=0
    family_unint_plus = 0
    family_unint_mins =0
    family_unint_double = 0
    alt_family_plus=0
    alt_family_mins =0
    alt_family_double = 0
    family_alt_percentage = collections.defaultdict(list)
    family_strand = collections.defaultdict(list)
    family_size = collections.defaultdict(list)
    non_alt_strand_plus=0
    non_alt_strand_mins=0
    ratio_1_plus=0
    ratio_1_minus=0
    ratio_1_double=0
    ratio = 0
    ratio_1_number = 0
    family_number = 0
    ratio_1_size=collections.defaultdict(list)
    ratio_variant =0
    if not alt:
        #print "Total family info for  non-alt genomic site"
        for family_tag,baselist in family[0].items():
            depth_raw +=len(baselist)
            if len(baselist)>= min_family_size:
                family_number = family_number+1
                family_size[family_tag].append(len(baselist))
        return [all_family_number,family_number,family_size,depth_raw]

    for family_tag,baselist in family[0].items():
        #print "family info for a mutation candidate site"
        if len(baselist)>= min_family_size:
            family_number = family_number+1
            family_size[family_tag].append(len(baselist))
            family_alt = baselist.count(alt)
            var_raw += family_alt  ###this var_raw is the number of alt raw reads
            depth_raw +=len(baselist) ##this depth_raw is the number of all reads in this position
            if family_alt:## calculate every family inner feature:family_alt_percentage & family_strand
                ratio = float(family_alt)/len(baselist)
                family_alt_percentage[family_tag].append(ratio)
                family_alt_percentage[family_tag].append(family_alt)
                family_alt_percentage[family_tag].append(len(baselist))
                family_has_alt +=1 
                family_strand_distribution = [family[1][family_tag][x] for x,y in enumerate(baselist) if y==alt]
                family_strand[family_tag].append(family_strand_distribution )
                family_plus+=family_strand_distribution.count("+")
                family_minus+=family_strand_distribution.count("-")
                ###family level plus and minus distribution
                if ratio >= cutoff:
                    ratio_1_number =ratio_1_number+1
                    ratio_1_size[family_tag].append(len(baselist))
                else:
                    ratio_variant = ratio_variant+1
                if set(family[1][family_tag]).issubset(["+"]):#1 represents plus 2 represent minus 3 represent double
                    family_alt_percentage[family_tag].append(1)
                    alt_family_plus = alt_family_plus+1
                    if ratio >= cutoff:
                        ratio_1_plus = ratio_1_plus+1
                elif set(family[1][family_tag]).issubset(["-"]):
                    family_alt_percentage[family_tag].append(2)
                    alt_family_mins = alt_family_mins+1
                    if ratio >= cutoff:
                        ratio_1_minus = ratio_1_minus+1
                else:
                    family_alt_percentage[family_tag].append(3)
                    alt_family_double = alt_family_double+1
                    if ratio >= cutoff:
                        ratio_1_double = ratio_1_double+1
            else:
                if set(family[1][family_tag]).issubset(["+"]):
                    family_unint_plus = family_unint_plus+1
                elif set(family[1][family_tag]).issubset(["-"]):
                    family_unint_mins = family_unint_mins+1
                else:
                    family_unint_double = family_unint_double+1
                non_alt_strand_plus += family[1][family_tag].count("+")
                non_alt_strand_mins +=family[1][family_tag].count("-")
                #family_context[family_tag].append(family[3][family_tag][pos])   
            #to see whether is a var family
        else:
            depth_raw+=len(baselist)
            var_raw+=baselist.count(alt)
            if baselist[0]==alt:
                var_single+=1
                if set(family[1][family_tag]).issubset(['+']):
                    var_single_plus+=1
                else:
                    var_single_minus+=1 
            else:
                non_var_single_plus += family[1][family_tag].count("+")
                non_var_single_minus+= family[1][family_tag].count("-")
                non_alt_strand_plus +=family[1][family_tag].count("+")
                non_alt_strand_mins +=family[1][family_tag].count("-")
    if SUMMARY:
        return [var_raw,depth_raw,ratio_1_number,ratio_1_plus,ratio_1_minus,ratio_1_double,ratio_variant,var_single,family_number]
    else:
        return[var_raw,depth_raw,family_has_alt,var_single,all_family_number,family_plus,family_minus,var_single_plus,var_single_minus,family_plus+var_single_plus,family_minus+var_single_minus,non_alt_strand_plus,non_alt_strand_mins,family_strand,family_alt_percentage,family_unint_plus,family_unint_mins,family_unint_double,alt_family_plus,alt_family_mins,alt_family_double,non_var_single_plus,non_var_single_minus,ratio_1_plus,ratio_1_minus,ratio_1_double,ratio_1_number,family_size,family_number,ratio_1_size,ratio_variant]

def xopen(filename, mode='r'):
    assert isinstance(filename, basestring)

    if filename == '-':
        if 'r' in mode:
            return sys.stdin
    return open(filename, mode)

def main():
    ap = argparse.ArgumentParser()
     
    ap.add_argument("-s1","--switch1",help="optinal length filter;default is open",type=ast.literal_eval,nargs="?",const=True,default=True)
    ap.add_argument("-s2","--switch2",help="optinal R1,R2 base filter;default is open",type=ast.literal_eval,nargs="?",const=True,default=True)
    ap.add_argument("-bq","--basequality",help="mininum base quality;default is Q30",nargs="?", const=30, type=int,default=30)
    ap.add_argument("-mq","--mappingquality",help="mininum base quality;default is Q30",nargs="?", const=30, type=int,default=30)
    ap.add_argument("-fsize","--familysize",help="mininum family size;default is 2",nargs="?",const=2,type=int,default=2)
    ap.add_argument("-fratio","--familyratio",help="mininum family ratio;default is 1.0",nargs="?",const=1.0,type=float,default=1.0)
    ap.add_argument("-sum","--summary",help="get summary of family statistics;default is True",type=ast.literal_eval,nargs="?",const=True,default=True)
    ap.add_argument("-tag","--UMI",help="UMI data;default is false",type=ast.literal_eval,nargs="?",const=False,default=False)
    ap.add_argument("bamfilepath",help="bamfile_path",type=str)
    if "-I" in sys.argv[1:] or "--Infile" in sys.argv[1:]:
        ap.add_argument("-I","--Infile",help="VCF format file",type=str)
        args=ap.parse_args()
    elif "-ih" in sys.argv[1:] or "--ihgvs" in sys.argv[1]:
        ap.add_argument("-ih","--ihgvs",help="file following HGVS standards:7:55259515T>C",type=str)
        args = ap.parse_args()
    else:
        ap.add_argument("-is","--isite", help = "single hgvs site",type=str)
        args = ap.parse_args()

    samfile = pysam.AlignmentFile(args.bamfilepath, "rb" )
    print "\t".join(["var_raw","depth_raw","f=1.0","+","-","+/-","f<1.0","var_single","family_number"])

    if 'Infile' in vars(args).keys():
       with xopen(args.Infile,'r') as f:
           for line in f:
               if line.startswith("#"):
                   continue
               data = line.strip().split("\t")
               family_info = get_family(samfile,data[0],int(data[1]),data[4],args.switch1,args.switch2,args.basequality,args.mappingquality,args.UMI)
               result = get_var_info(family_info,data[4],args.familysize,args.familyratio,args.summary)
               print "\t".join([str(x) for x in result]) 
             
    elif "ihgvs" in vars(args).keys():
       pattern = re.compile("[0-9]+|[ATCG]$")
       with xopen(args.ihgvs,'r') as f:
           for line in f:
               in_data =pattern.findall(line.strip())
               family_info = get_family(samfile,in_data[0],int(in_data[1]),in_data[2],args.switch1,args.switch2,args.basequality,args.UMI)
               result = get_var_info(family_info,in_data[2],args.familysize,args.familyratio,args.summary)
               print "\t".join([str(x) for x in result]) 
    else:
       pattern = re.compile("[0-9]+|[ATCG]$")
       in_data =pattern.findall(args.isite)
       family_info = get_family(samfile,in_data[0],int(in_data[1]),in_data[2],args.switch1,args.switch2,args.basequality,args.UMI)
       result = get_var_info(family_info,in_data[2],args.familysize,args.familyratio,args.summary)
       print "\t".join([str(x) for x in result]) 
    samfile.close()
    
if __name__=='__main__':
    main()
