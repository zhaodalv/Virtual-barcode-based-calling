
python 2.7

####This file is the instructions for use virtual family step and polishing step
1.1: python virtual_family.py -is '1:115256536C>T'(hgvs format) path_to_bam


2.1: For file input both in stdin or as a input file(format:VCF OR HGVS): python virtual_family.py -ih hgvs_file OR cat hgvs_file | python virtual_family.py -ih - path_to_bam |||| python virtual_family.py -I VCF_file path_to_bam OR cat VCF_file | python virtual_family -I - path_to_bam


1.1 Required arguments
-is ISITE, --isite ISITE single hgvs site


2.1 Required arguments for 
-ih, --ihgvs, file following HGVS standards:7:55259515T>C
OR -I, --Infile , VCF format file

all required arguments
bamfilepath           bamfile_path

Optional arguments:

-s1 [SWITCH1], --switch1 [SWITCH1]
                        optinal length filter;default is open

-s2 [SWITCH2], --switch2 [SWITCH2]
                        optinal R1,R2 base filter;default is open

-bq [BASEQUALITY], --basequality [BASEQUALITY]
                        mininum base quality;default is Q30

-fsize [FAMILYSIZE], --familysize [FAMILYSIZE]
                        mininum family size;default is 2

-fratio [FAMILYRATIO], --familyratio [FAMILYRATIO]
                        mininum family ratio;default is 1.0

-sum [SUMMARY], --summary [SUMMARY]
                        get summary of family statistics;default is True

-tag [UMI], --UMI [UMI]
                        UMI data;default is false

-mq [MAPPINGQUALITY], --mappingquality [MAPPINGQUALITY]
                        mininum base quality;default is Q30

Output: default summary output (-sum): 
columns	definition
var_raw	 variant read numbers
depth_raw	total sequence depth
f=1.0	f=1.0 family numbers
+	number of f=1.0 family consisting of reads from DNA plus strand
-	number of f=1.0 family consisting of reads from DNA minus strand
+/-	number of f=1.0 family consisting of reads both from DNA plus and minus strand
f<1.0	f<1.0 family numbers
var_single	variant singleton numbers
family_number total family numbers with size>2 (default and this value can change by argument -fsize)

detailed output (-sum False):
var_raw          variant read numbers
depth_raw        total sequence depth
family_has_alt   f=1.0 familynumber + f<1.0 family number
var_single       variant singleton numbers 
all_family_number    total family number for any size
family_plus    number of alt reads from DNA plus strand
family_minus   number of alt reads from DNA minus strand
var_single_plus  varint singleton from DNA plus strand
var_single_minus  variant singleton from DNA minus strand
non_alt_strand_plus  number of wild reads from plus strand 
non_alt_strand_mins  number of wild reads from minus strand
family_strand        dict for every family_has_alt family strand consititution: '55259356,185,': [['-', '+', '+', '+', '+', '-', '-', '-', '+', '+', '+', '+', '-', '+']
family_alt_percentage  dict for every family_has_alt family f value,alt number,family size,strand constitution (1 represents plus 2 represent minus 3 represent double): 55259356,185,': [1.0, 14, 14, 3]
family_unint_plus  number of f=0 family consisting of reads from DNA plus strand
family_unint_mins  number of f=0 family consisting of reads from DNA minus strand
family_unint_double number of f=0 family consisting of reads both from DNA plus and minus
alt_family_plus   number of f>0 family consisting of reads from DNA plus strand
alt_family_mins   number of f>0 family consisting of reads from DNA minus strand
alt_family_double number of f>0 family consisting of reads both from DNA plus and minus
non_var_single_plus  singleton from DNA plus strand
non_var_single_minus  variant singleton from DNA minus strand
ratio_1_plus    number of f=1.0(or f>familyratio)  consisting of reads from DNA plus strand 
ratio_1_minus   number of f=1.0(or f>familyratio)  consisting of reads from DNA minus strand
ratio_1_double  number of f=1.0(or f>familyratio)  consisting of reads both from DNA plus and minus strand
ratio_1_number  number of f=1.0(or f>familyratio) family
family_size     dict for size of every family:'55259347,173,': [5] 
family_number   number of family with size >=2((default and this value can change by argument -fsize))
ratio_1_size    dict for size of every f=1.0 famiy (or f>familyratio)
ratio_variant   f<1.0 family numbers


######polishing #####
usage: polishing_step.py [-h] AF_file

AF_file format: dataframe or matrix: row is polishing site AF value, columns is every sample.
###S1,S2,S3,S4###
###0.1,0.2,0.3,0.4###

output:polishing site; cutoff value;max value
('11', '533838', 'T') 	0.219766800011 	0.0178346861781

