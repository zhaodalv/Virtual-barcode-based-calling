#!/usr/bin/env sh



###Variables description###
##PATH_PICARD is the path of picard.jar; GATK_PATH is the path of GenomeAnalysisTK.jar; HUMAN_REF is the path of human reference file (human_g1k_v37_decoy.fasta); WKDIR is the directory for saving result; DPSNP is the file path of dbsnp_138.b37.vcf; INDEL_1000G is the path of 1000G_phase1.indels.b37.vcf; INDEL_1000G_GD is the path of Mills_and_1000G_gold_standard.indels.b37.vcf;BAM_FILE is the input bam file path; OUT_PUT_BAM is the out put bam file name;${PATH_TO_LIB} is the path of libIntelDeflater.so; PATH_readcount is the path of bam-readcount;ABRA2 is path of ABRA2; SAMPLE_NAME is the name of sample with mutations; NORMAL_NAME is the corresponding sample without mutations;

if [ -d ${WKDIR}/MarkDuplicates ]; then
mkdir -p ${WKDIR}/MarkDuplicates/${SAMPLE_NAME}/bam
mkdir -p ${WKDIR}/MarkDuplicates/${SAMPLE_NAME}/tmp
mkdir -p ${WKDIR}/MarkDuplicates/${NORMAL_NAME}/bam
mkdir -p ${WKDIR}/MarkDuplicates/${NORMAL_NAME}/tmp
fi


java -Xmx20g -XX:ParallelGCThreads=2 -jar ${PATH_PICARD} MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${BAM_FILE} OUTPUT=${WKDIR}/MarkDuplicates/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.MKD.bam METRICS_FILE=${WKDIR}/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.metrics TMP_DIR=${WKDIR}/${SAMPLE_NAME}/tmp

samtools.1.3 index ${WKDIR}/MarkDuplicates/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.MKD.bam
###second step of preprocessing is BaseRecalibrator
if [ -d ${WKDIR}/preprosessing_BRSQ ]; then
mkdir -p ${WKDIR}/preprosessing_BRSQ/${SAMPLE_NAME}/bam
mkdir -p ${WKDIR}/preprosessing_BRSQ/${SAMPLE_NAME}/tmp
mkdir -p ${WKDIR}/preprosessing_BRSQ/${NORMAL_NAME}/bam
mkdir -p ${WKDIR}/preprosessing_BRSQ/${NORMAL_NAME}/tmp
fi


java -Xmx20g -XX:ParallelGCThreads=2 -jar ${GATK_PATH} -T BaseRecalibrator -R ${HUMAN_REF} -I ${WKDIR}/MarkDuplicates/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.MKD.bam -knownSites ${DPSNP} -knownSites ${INDEL_1000G} -knownSites ${INDEL_1000G_GD} -o ${WKDIR}/preprosessing_BRSQ/${SAMPLE_NAME}/tmp/${SAMPLE_NAME}.table 


java -Xmx20g -XX:ParallelGCThreads=2 -jar ${GATK_PATH}  -T PrintReads -R ${HUMAN_REF} -I ${WKDIR}/MarkDuplicates/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.MKD.bam -BQSR ${WKDIR}/preprosessing_BRSQ/${SAMPLE_NAME}/tmp/${SAMPLE_NAME}.table -o ${WKDIR}/preprosessing_BRSQ/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.BRSQ.bam

#third step of preprocessing is realignment 
if [ -d ${WKDIR}/preprosessing_realignment ]; then
mkdir -p ${WKDIR}/preprosessing_realignment/${SAMPLE_NAME}/bam
mkdir -p ${WKDIR}/preprosessing_realignment/${SAMPLE_NAME}/tmp
mkdir -p ${WKDIR}/preprosessing_realignment/${NORMAL_NAME}/bam
mkdir -p ${WKDIR}/preprosessing_realignment/${NORMAL_NAME}/tmp
fi
java -Xmx4g -XX:ParallelGCThreads=2 -jar ${GATK_PATH} -T IndelRealigner -I ${WKDIR}/preprosessing_BRSQ/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.BRSQ.bam -o ${WKDIR}/preprosessing_realignment/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.BRSQ.REAlig.bam  -targetIntervals ${WKDIR}/preprosessing_realignment/tmp/${SAMPLE_NAME}.realignment.interval  -R ${HUMAN_REF} -known ${INDEL_1000G} -known ${INDEL_1000G_GD}


##MUTECT pipleline mutect-1.1.7.jar
if [ -d ${WKDIR}/Mutect ]; then
mkdir -p ${WKDIR}/Mutect
fi

java -jar ${MUTECT_PATH} --analysis_type MuTect --reference_sequence ${HUMAN_REF} --cosmic ${CosmicCodingMuts} --dbsnp ${DPSNP} --intervals ${BED}  --input_file:normal ${WKDIR}/preprosessing_realignment/${NORMAL_NAME}/bam/${NORMAL_NAME}.BRSQ.REAlig.bam --input_file:tumor ${WKDIR}/preprosessing_realignment/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.BRSQ.REAlig.bam --out ${WKDIR}/Mutect/${SAMPLE_NAME}.txt 

##MUTECT2 pipleline 
if [ -d ${WKDIR}/Mutect2 ]; then
mkdir -p ${WKDIR}/Mutect2
fi

java -jar ${GATK_PATH} -T MuTect2 -R ${HUMAN_REF} --cosmic ${CosmicCodingMuts} --dbsnp ${DPSNP} -L ${BED} --max_alt_alleles_in_normal_count 20 -mbq 30 -I:normal ${WKDIR}/preprosessing_realignment/${NORMAL_NAME}/bam/${NORMAL_NAME}.BRSQ.REAlig.bam -I:tumor ${WKDIR}/preprosessing_realignment/${SAMPLE_NAME}/bam/${SAMPLE_NAME}.BRSQ.REAlig.bam -o ${WKDIR}/Mutect2/${SAMPLE_NAME}.vcf


#VARSCAN2 pipleline
if [ -d ${WKDIR}/Varscan2 ]; then
mkdir -p ${WKDIR}/Varscan2/pileup/${SAMPLE_NAME}
mkdir -p ${WKDIR}/Varscan2/normal_pileup/${NORMAL_NAME}
mkdir -p ${WKDIR}/Varscan2/result
fi
samtools.1.3 mpileup -q 1 -Q 30 -d 100000 -f ${HUMAN_REF} -l ${BED} ${BAM_FILE}  -o ${WKDIR}/Varscan2/pileup/${SAMPLE_NAME}/${SAMPLE_NAME}.pileup
java -jar VarScan.v2.3.9.jar somatic ${WKDIR}/Varscan2/normal_pileup/${NORMAL_NAME}/${NORMAL_NAME}.pileup ${WKDIR}/Varscan2/pileup/${SAMPLE_NAME}/${SAMPLE_NAME}.pileup ${WKDIR}/Varscan2/result -min-var-freq 0.0001 -min-coverage 3 -min-coverage-normal 3 -min-coverage-tumor 3 -output-vcf 1

#SINVICT pipleline
if [ -d ${WKDIR}/SiNVICT ]; then
mkdir -p ${WKDIR}/SiNVICT/abra2/${SAMPLE_NAME}
mkdir -p ${WKDIR}/SiNVICT/readcount/${SAMPLE_NAME}
mkdir -p ${WKDIR}/Varscan2/result
fi

java -Dsamjdk.intel_deflater_so_path=${PATH_TO_LIB}/libIntelDeflater.so -Xmx20g -jar ${ABRA2}/abra2-2.06.jar --in ${BAM_FILE} --out ${WKDIR}/SiNVICT/abra2/${SAMPLE_NAME}/${SAMPLE_NAME}.abra2.bam  --ref ${HUMAN_REF} --targets ${BED}  --threads 8

samtools.1.3 index ${WKDIR}/SiNVICT/abra2/${SAMPLE_NAME}/${SAMPLE_NAME}.abra2.bam

${PATH_readcount}/bam-readcount -b 30 -d 100000 -w 0 -f ${HUMAN_REF} -l ${BED} ${WKDIR}/SiNVICT/abra2/${SAMPLE_NAME}/${SAMPLE_NAME}.abra2.bam >${WKDIR}/SiNVICT/readcount/${SAMPLE_NAME}/${SAMPLE_NAME}.readcount

sinvict -m 1000 --qscore-cutoff 20 -t ${WKDIR}/SiNVICT/readcount/${SAMPLE_NAME}/${SAMPLE_NAME}.readcount -o ${WKDIR}/SiNVICT/result/${SAMPLE_NAME}
