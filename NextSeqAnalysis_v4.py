__author__ = 'HOI QIANGZE'
import argparse
import os
import glob
import sys
import SoftwareConfiguration
import datetime
import re

###########VERSION ###################
version=4.0
######################################


# time=datetime.date.today()
# date=time.strftime("%y%m%d")
timenow=datetime.datetime.now().strftime("%y%m%d.%H%M")
version= str(SoftwareConfiguration.version)
SOFTWARE_FOLDER=str(SoftwareConfiguration.software_directory)
REFERENCE_FOLDER=str(SoftwareConfiguration.reference_directory)
BWA_EXE=str(SoftwareConfiguration.bwa_exe)
PICARD_JAR=str(SoftwareConfiguration.picard_jar)
GATK_JAR=str(SoftwareConfiguration.gatk_jar)
SNPEFF_JAR=str(SoftwareConfiguration.snpeff_jar)
SNPEFF_CONF=str(SoftwareConfiguration.snpEff_config)
SNPSIFT_JAR=str(SoftwareConfiguration.snpsift_jar)
emission_conf=str(SoftwareConfiguration.emission_conf)
calling_conf=str(SoftwareConfiguration.calling_conf)
bamUtil=str(SoftwareConfiguration.bamUtil)
BSQR_VCF=str(SoftwareConfiguration.ThirtySamples)
omni=REFERENCE_FOLDER+"/"+str(SoftwareConfiguration.omni)
vcf_1000g_snp=REFERENCE_FOLDER+"/"+str(SoftwareConfiguration.oneKg_snp)
vcf_1000g_indel=REFERENCE_FOLDER+"/"+str(SoftwareConfiguration.oneKg_indel)
dbsnp=REFERENCE_FOLDER+"/"+str(SoftwareConfiguration.dbsnp)
hapmap=REFERENCE_FOLDER+"/"+str(SoftwareConfiguration.hapmap)
mills=REFERENCE_FOLDER+"/"+str(SoftwareConfiguration.mills)
#fetch options

#declare some variables
parser = argparse.ArgumentParser(description="Use the following options to run the preparation step for Next Seq Analysis")
parser.add_argument("--name" , '-n', help="Name of experiment. This will be used to create the folder")
parser.add_argument("--id" , '-i', help="Run Number #.", default=1)
parser.add_argument("--input" , '--in', help="Input folder containing fastq.")
parser.add_argument("--ref" , '-r', help="reference. eg hg19")
parser.add_argument("--genelist" , '-g', help="genes of interest to filter. have to be a file path")
parser.add_argument("--bed" , '-b', help="bed file with chroms start end coordinates")

args = vars(parser.parse_args())

experimentName=args['name']
runNo=args['id']
reference=args['ref']
if not reference:
    reference='hg19'
if not args['genelist']:
    print("Please provide genelist with -g option")
    sys.exit()
else:
    gene_list=args['genelist']

homefolder='/home/ubuntu'
referencefolder='/mnt1/references'
workingdir='/mnt2/output' # for now use this
#workingdir=os.getcwd() # for now use this
experimentPath=workingdir+"/"+experimentName
rootpath=""
#glob all the fastq files in the given input folder
fastq_folder=args['input']
fastq_folder_name=os.path.basename(fastq_folder)
#fastqfile=glob.glob(fastq_folder+"/*fastq.gz")

bedfile=args['bed']

#fix to get read 1 in the array first
fastq_read1=" ".join(glob.glob(fastq_folder+"/*R1.fastq.gz"))
fastq_read2=" ".join(glob.glob(fastq_folder+"/*R2.fastq.gz"))
fastq_prefix = []
m = re.match("(.*)\.fastq.gz",os.path.basename(fastq_read1))
fastq_prefix.append(m.group(1))
m = re.match("(.*)\.fastq.gz",os.path.basename(fastq_read2))
fastq_prefix.append(m.group(1))
print (fastq_prefix)
#print (fastqfile)
#for files in fastqfile:
#    print("found this file %s from input folder %s" %(files,fastq_folder))
#    m = re.match("(.*)\.fastq.gz",os.path.basename(files))
#    if m:
#        fastq_prefix.append(m.group(1))
#        print (m.group(1))

if not os.path.exists(experimentPath):
    os.makedirs(experimentPath)

rootpath=experimentPath+"/"+fastq_folder_name+"_pipeline_"+version
if not os.path.exists(rootpath):
    os.makedirs(rootpath)
    os.makedirs(rootpath+"/INPUT_FILES")
    os.makedirs(rootpath+"/CONFIG_FILES")
    os.makedirs(rootpath+"/MAPPING")
    os.makedirs(rootpath+"/GATK_PROCESSING")
    os.makedirs(rootpath+"/FILTERING")
    os.makedirs(rootpath+"/METRICS")


runlog = workingdir+"/"+"run.log"
runlogwriter=open(runlog,"a")
runlogwriter.write("nohup snakemake -s %s/NextSeq.snakefile > %s/logger.txt &\n" %(rootpath,rootpath))
print("Root output folder at %s" %rootpath)
runlogwriter.write("nohup snakemake -s %s/NextSeq.snakefile > %s/logger.txt &\n" %(rootpath,rootpath))
print("Root output folder at %s" %rootpath)
    # sys.exit("This folder exist. Please make use of another name or check %s" %experimentPath)
#here i copy the fastq to INPUT FILES FOLDER
INPUT=rootpath+"/INPUT_FILES"
#DEBUG: Commend out for debugging
if os.path.exists(fastq_folder):
    os.system('ln -s '+fastq_folder+'/*fastq.gz '+INPUT)
else:
    print("FASTQ FILES NOT DETECTED. PLEASE ENSURE YOU DOWNLOAD THEM BEFORE YOU START")
print('ln -s '+fastq_folder+'/*fastq.gz '+INPUT)
#print ('cp '+fastq_folder+'/*fastq.gz '+INPUT) 
print('make file at %s/NextSeq.snakefile' %rootpath)
# here we write the make file
snakefile=rootpath+"/NextSeq.snakefile"
f=open(snakefile,"w")

#declare some variables. will be easier to write
make_VERSION="VERSION=\'%s'\n" %version
make_SAMPLES='SAMPLES = \''+" ".join(fastq_prefix) + "\'.split() \n"
make_FIRSTSAMPLENAME= 'SAMPLENAME = \''+fastq_prefix[0]+"\'\n"
make_ROOT_FOLDER= 'ROOT_FOLDER = \''+rootpath + "\'\n"
make_PICARD_TMP= 'PICARD_TMP = \''+rootpath + "/tmp\'\n"
make_REF_NAME= "REF_NAME = \'%s\'\n" %reference
make_REFERENCE_DIR="REFERENCE_DIR = \'%s/%s/ucsc.hg19.fasta\'\n" %(referencefolder,reference)
make_SOFTWARE_FOLDER_PATH= 'SOFTWARE_FOLDER =  \''+SOFTWARE_FOLDER + "\'\n"
make_GETONTARGETDEPTH= "ONTARGET_DEPTH=  \'~/scripts/calculate_on_target.py\'\n"
make_FASTQ_DIR= 'FASTQ_DIR = \''+INPUT+ "\'\n"
make_BWA_OUTPUT='BWA_OUTPUT = \''+rootpath+"/MAPPING\'\n"
make_GATK_OUTPUT='GATK_OUTPUT = \''+rootpath+"/GATK_PROCESSING\'\n"
make_FILTERING='FILTER_OUTPUT = \''+rootpath+"/FILTERING\'\n"
make_KNOWN_INDELS="KNOWN_INDELS=\'%s\'\n" %mills
make_DBSNP="DBSNP=\'%s'\n" %dbsnp
make_1000Gsnp="vcf_1000g_snp=\'%s'\n" %vcf_1000g_snp
make_1000Gindel="vcf_1000g_indel=\'%s'\n" %vcf_1000g_indel
make_OMNI="OMNI=\'%s'\n" %omni
make_HAPMAP="HAPMAP=\'%s'\n" %hapmap
make_bamUtil="BAMUTIL=\'%s'\n" %bamUtil
make_BSQR_VCF="BSQR_VCF=\'%s'\n" %BSQR_VCF
make_bed="BED=\'%s'\n" %bedfile

f.write(make_VERSION)
f.write(make_SAMPLES)
f.write(make_FIRSTSAMPLENAME)
f.write(make_ROOT_FOLDER)
f.write(make_PICARD_TMP)
f.write(make_REF_NAME)
f.write(make_REFERENCE_DIR)
f.write(make_SOFTWARE_FOLDER_PATH)
f.write(make_GETONTARGETDEPTH)
f.write(make_FASTQ_DIR)
f.write(make_BWA_OUTPUT)
f.write(make_GATK_OUTPUT)
f.write(make_FILTERING)
f.write(make_KNOWN_INDELS)
f.write(make_DBSNP)
f.write(make_1000Gsnp)
f.write(make_1000Gindel)
f.write(make_OMNI)
f.write(make_HAPMAP)
f.write(make_GETONTARGETDEPTH)
f.write(make_bamUtil)
f.write(make_BSQR_VCF)
f.write(make_bed)


rule_all="""rule all:
  input: FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.IMPACT.DP20.GENESET.DBSNP.DBNSFP.csv" \n """


rule_make_picard_tmp="""rule make_picard_tmp :
  output: PICARD_TMP
  shell: "mkdir -p {output}"\n"""

rule_map_reads="""rule map_reads :
  input: ref = REFERENCE_DIR, reads1 = FASTQ_DIR + "/"+SAMPLES[0]+".fastq.gz", reads2= FASTQ_DIR+"/"+SAMPLES[1]+".fastq.gz", dummy=rules.make_picard_tmp.output
  output: temp(BWA_OUTPUT + "/"+SAMPLES[0] + ".bam")
  shell: "%s mem -M -R'@RG\\\\tID:group1\\\\tSM:{SAMPLENAME}\\\\tPL:illumina' {input.ref} {input.reads1} {input.reads2} > {output}"\n""" %BWA_EXE

rule_sort_sam_convert_bam= """rule sort_sam_convert_bam :
  input: mapped=rules.map_reads.output ,tmp=PICARD_TMP
  output: bam=temp(BWA_OUTPUT + "/"+SAMPLES[0] + ".sorted.bam")
  run:
    shell("java -Djava.io.tmpdir={input.tmp} -Xmx4g -jar %s SortSam INPUT={input.mapped} OUTPUT={output.bam} SORT_ORDER=coordinate")\n"""%PICARD_JAR

rule_dedup = """rule dedup:
  input: bamfile=rules.sort_sam_convert_bam.output.bam,tmp=PICARD_TMP
  output: dedup=temp(BWA_OUTPUT + "/"+SAMPLES[0] + ".sorted.dedup.bam"),metrics=BWA_OUTPUT+"/metrics"
  run:
    shell("java -Djava.io.tmpdir={input.tmp} -Xmx4g -jar %s MarkDuplicates INPUT={input.bamfile} OUTPUT={output.dedup} METRICS_FILE={output.metrics}")\n""" %PICARD_JAR

rule_index_bam= """rule index_bam:
  input: dedup=rules.dedup.output.dedup,tmp=PICARD_TMP
  output: index=BWA_OUTPUT + "/"+SAMPLES[0] + ".sorted.dedup.bam.bai"
  run:
    shell("java -Djava.io.tmpdir={input.tmp} -Xmx4g -jar %s BuildBamIndex INPUT={input.dedup} OUTPUT={output.index} ")\n""" %PICARD_JAR

rule_realignertargetcreator= """rule realignertargetcreator:
  input: ref=REFERENCE_DIR , known=KNOWN_INDELS, dummy=rules.index_bam.output.index
  output: GATK_OUTPUT +"/output.intervals"
  shell: "java -Xmx4g -jar %s -T RealignerTargetCreator -R {input.ref} -o {output} -known {input.known} "\n""" %GATK_JAR

rule_indel_realigner= """rule indel_realigner:
  input: intervals=rules.realignertargetcreator.output , bam=rules.dedup.output.dedup , ref=REFERENCE_DIR , known=KNOWN_INDELS
  output: temp(GATK_OUTPUT+ "/"+SAMPLES[0] + ".sorted.dedup.realigned.bam")
  shell: "java -Xmx4g -jar %s -T IndelRealigner -I {input.bam} -R {input.ref} -targetIntervals {input.intervals} -o {output} -known {input.known} --consensusDeterminationModel KNOWNS_ONLY"\n""" %GATK_JAR

rule_base_recalibration_first= """rule base_recalibration_first:
  input: bam=rules.indel_realigner.output , ref=REFERENCE_DIR , known_indel=KNOWN_INDELS, known_snp=DBSNP
  output: GATK_OUTPUT+ "/recal_data.table"
  shell: "java -Xmx4g -jar %s -T BaseRecalibrator -I {input.bam} -R {input.ref} -knownSites {BSQR_VCF} -knownSites {input.known_snp} -knownSites {input.known_indel} -o {output}"\n""" %GATK_JAR

rule_base_recalibration_second= """rule base_recalibration_second:
  input: bam=rules.indel_realigner.output , ref=REFERENCE_DIR , known=KNOWN_INDELS, known_sites=DBSNP, bqsr= rules.base_recalibration_first.output
  output: GATK_OUTPUT+ "/post_recal_data.table"
  shell: "java -Xmx4g -jar %s -T BaseRecalibrator -I {input.bam} -R {input.ref} -knownSites {BSQR_VCF} -knownSites {input.known} -knownSites {input.known_sites} -BQSR {input.bqsr} -o {output}"\n""" %GATK_JAR

rule_analyze_cov= """rule analyze_cov :
  input: before=rules.base_recalibration_first.output, after = rules.base_recalibration_second.output , ref=REFERENCE_DIR 
  output: GATK_OUTPUT+ "/recalibration_plots.pdf"
  shell: "java -Xmx4g -jar %s -T AnalyzeCovariates -R {input.ref} -before {input.before} -after {input.after} -plots {output}"\n""" %GATK_JAR

rule_apply_recal= """rule apply_recal:
  input: bqsr=rules.base_recalibration_first.output, bam=rules.indel_realigner.output , ref=REFERENCE_DIR, dependency=rules.base_recalibration_second.output 
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + ".sorted.dedup.realigned.recal.bam" 
  shell: "java -Xmx4g -jar %s -T PrintReads -R {input.ref} -I {input.bam} -BQSR {input.bqsr} -o {output}"\n""" %GATK_JAR

#python3 ~/scripts/calculate_on_target.py /mnt2/output/LB_study/bamfiles/NexteraRapidCapture_Exome_TargetedRegions_v1.2.bed /mnt2/output/LB_study/bamfiles/LB175_S10/LB175_S10_combined_R1.sorted.dedup.realigned.recal.bam LB175_S10_combined_R1.sorted.dedup.realigned.recal.bai outfile ~/installations/bamUtil_1.0.13/bamUtil/bin/bam


rule_get_ontargetdepth= """rule get_ontargetdepth:
  input: bam=rules.apply_recal.output, bamindex=rules.apply_recal.output+".bai" 
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + "ontarget.depth"
  shell: "python {ONTARGET_DEPTH} {BED} {input.bam} {input.bamindex} {output} {BAMUTIL}"\n"""

rule_haplotype_caller= """rule haplotype_caller:
  input: bam=rules.apply_recal.output , ref=REFERENCE_DIR
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + "raw_variants.vcf" 
  shell: "java -Xmx4g -jar %s -T HaplotypeCaller -R {input.ref} -I {input.bam} -variant_index_type LINEAR -variant_index_parameter 128000 --genotyping_mode DISCOVERY -stand_emit_conf %s -stand_call_conf %s --emitRefConfidence GVCF -o {output}"\n""" %(GATK_JAR, emission_conf, calling_conf)

rule_VariantAnnotator= """rule VariantAnnotator:
  input: vcf=rules.haplotype_caller.output
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + "raw_variants.annotated.vcf" 
  shell: "java -Xmx4g -jar %s -T VariantAnnotator -R {REFERENCE_DIR} --variant {input.vcf} -A StrandOddsRatio -A FisherStrand -o {output}"\n""" %GATK_JAR

rule_VariantRecalibrator= """rule VariantRecalibrator:
  input: HC_vcf=rules.VariantAnnotator.output, ref=REFERENCE_DIR, hapmap=HAPMAP,omni=OMNI, oneK_snp=vcf_1000g_snp, dbsnp=DBSNP
  output: recal_file=GATK_OUTPUT+ "/"+"recalibrate_SNP.recal", tranches_file=GATK_OUTPUT+ "/"+"recalibrate_SNP.tranches", rscriptfile=GATK_OUTPUT+ "/"+"recalibrate_SNP_plots.R"
  shell: "java -Xmx4g -jar %s -T VariantRecalibrator -R {input.ref} -input {input.HC_vcf} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} -resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} -resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.oneK_snp} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} -an DP -an QD -an FS -an MQ -an ReadPosRankSum -an InbreedingCoeff -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile {output.recal_file} -tranchesFile {output.tranches_file} -rscriptFile {output.rscriptfile} "\n""" %GATK_JAR

rule_apply_recalibration= """rule apply_recalibration:
  input: recalfile=rules.VariantRecalibrator.output.recal_file, vcf=rules.VariantAnnotator.output, tranches_file=rules.VariantRecalibrator.output.tranches_file
  output: applied_recal=GATK_OUTPUT+ "/"+SAMPLES[0] + "rcalibrated.vcf"
  shell: "java -Xmx4g -jar %s -T ApplyRecalibration -R {REFERENCE_DIR} -input {input.vcf} -mode SNP --ts_filter_level 99.0 -recalFile {input.recalfile} -tranchesFile {input.tranches_file} -o {output.applied_recal}"\n""" %GATK_JAR

rule_VariantRecalibrator_indel= """rule VariantRecalibrator_indel:
  input: ref=REFERENCE_DIR, mills=KNOWN_INDELS, HC_vcf=rules.apply_recalibration.output
  output: recal_file=GATK_OUTPUT+ "/"+"recalibrate_INDEL.recal", tranches_file=GATK_OUTPUT+ "/"+"recalibrate_INDEL.tranches", rscriptfile=GATK_OUTPUT+ "/"+"recalibrate_INDEL_plots.R"
  shell: "java -Xmx4g -jar %s -T VariantRecalibrator -R {input.ref} -input {input.HC_vcf} -resource:mills,known=true,training=true,truth=true,prior=12.0 {input.mills} -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile {output.recal_file} -tranchesFile {output.tranches_file} -rscriptFile {output.rscriptfile} "\n""" %GATK_JAR

rule_apply_recalibration_indel= """rule apply_recalibration_indel:
  input: vcf=rules.apply_recalibration.output , ref=REFERENCE_DIR, recalfile=rules.VariantRecalibrator_indel.output.recal_file, tranches_file=rules.VariantRecalibrator_indel.output.tranches_file
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + "recalibrated_variants.vcf" 
  shell: "java -Xmx4g -jar %s -T ApplyRecalibration -R {input.ref} -input {input.vcf} -mode INDEL --ts_filter_level 99.0 -recalFile {input.recalfile} -tranchesFile {input.tranches_file} -o {output}"\n""" %GATK_JAR

#Variant annotation section


rule_snpEff_annotate= """rule snpEff_annotate:
  input: vcf=rules.apply_recalibration.output.applied_recal
  #input: vcf=rules.VariantAnnotator.output
  output: vcf=FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.vcf",stats=FILTER_OUTPUT+"/"+SAMPLES[0]+".stats"
  shell: "java -Xmx4g -jar %s -c %s  -stats {output.stats} %s {input.vcf} > {output.vcf}"\n""" %(SNPEFF_JAR,SNPEFF_CONF,reference)

rule_filter_impact_variants= """rule filter_impact_variants:
  input: vcf=rules.snpEff_annotate.output.vcf
  output: vcf=FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.IMPACTFUL.DP20.GENESET.vcf"
  shell: "java -Xmx4g -jar %s filter --set %s \\"((ANN[*].IMPACT =  'MODERATE')|(ANN[*].IMPACT = 'HIGH'))  &(DP>=20) & (ANN[*].GENE in SET[0])\\" -f {input.vcf} >{output.vcf}"\n""" %(SNPSIFT_JAR,gene_list)

rule_filter_moderate_impact_variants= """rule filter_moderate_impact_variants:
  input: vcf=rules.snpEff_annotate.output.vcf,  check=rules.filter_high_impact_variants.output.vcf
  output: vcf=FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.MODERATE.DP20.GENESET.vcf"
  shell: "java -Xmx4g -jar %s filter --set %s \\"(ANN[*].IMPACT = 'MODERATE')  &(DP>=20) & (ANN[*].GENE in SET[0])\\" -f {input.vcf} >{output.vcf}"\n""" %(SNPSIFT_JAR,gene_list)

rule_annotate_impact_with_dbsnp= """rule annotate_impact_with_dbsnp:
  input: vcf=rules.filter_impact_variants.output, dbsnp=DBSNP
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.IMPACTFUL.DP20.GENESET.DBSNP.vcf"
  shell: "java -Xmx4g -jar %s annotate -id {input.dbsnp} {input.vcf} > {output}"\n""" %SNPSIFT_JAR

rule_annotate_moderate_impact_with_dbsnp= """rule annotate_moderate_impact_with_dbsnp:
  input: vcf=rules.filter_moderate_impact_variants.output, dbsnp=DBSNP
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.MODERATE.DP20.GENESET.DBSNP.vcf"
  shell: "java -Xmx4g -jar %s annotate -id {input.dbsnp} {input.vcf} > {output}"\n""" %SNPSIFT_JAR

rule_annotate_moderate_with_dbnsfp= """rule annotate_moderate_with_dbnsfp:
  input: vcf=rules.annotate_moderate_impact_with_dbsnp.output
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] +  ".SNPEFF.MODERATE.DP20.GENESET.DBSNP.DBNSFP.vcf"
  shell: "java -Xmx4g -jar %s dbnsfp -v {input.vcf} > {output}"\n""" %SNPSIFT_JAR

rule_annotate_with_dbnsfp= """rule annotate_with_dbnsfp:
  input: vcf=rules.annotate_impact_with_dbsnp.output
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] +  ".SNPEFF.IMPACT.DP20.GENESET.DBSNP.DBNSFP.vcf"
  shell: "java -Xmx4g -jar %s dbnsfp -v {input.vcf} > {output}"\n""" %SNPSIFT_JAR

rule_extract_moderate_to_csv= """rule extract_moderate_to_csv:
  input: vcf=rules.annotate_moderate_with_dbnsfp.output
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] +  ".SNPEFF.MODERATE.DP20.GENESET.DBSNP.DBNSFP.csv"
  shell: "java -Xmx4g -jar %s extractFields -s ',' -e '.' {input.vcf} 'CHROM' 'POS' 'REF' 'ALT' 'ID' 'DP' 'EFF[*].IMPACT' 'EFF[*].GENE' 'dbNSFP_SIFT_pred' 'dbNSFP_Polyphen2_HVAR_pred' 'dbNSFP_MutationTaster_pred' 'dbNSFP_1000Gp1_ASN_AF''dbNSFP_Uniprot_acc' > {output}"\n""" %SNPSIFT_JAR

rule_extract_to_csv= """rule extract_to_csv:
  input: vcf=rules.annotate_with_dbnsfp.output
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] +  ".SNPEFF.IMPACT.DP20.GENESET.DBSNP.DBNSFP.csv"
  shell: "java -Xmx4g -jar %s extractFields -s ',' -e '.' {input.vcf} 'CHROM' 'POS' 'REF' 'ALT' 'ID' 'DP' 'EFF[*].IMPACT' 'EFF[*].EFFECT' 'EFF[*].GENE' 'EFF[*].HGVS_P' 'EFF[*].HGVS_C' 'dbNSFP_SIFT_pred' 'dbNSFP_Polyphen2_HVAR_pred' 'dbNSFP_MutationTaster_pred' 'dbNSFP_1000Gp1_ASN_AF' > {output}"\n""" %SNPSIFT_JAR

rule_post_sneff_gatk_annotate= """rule post_sneff_gatk_annotate:
  input: snpeff_vcf=rules.snpEff_annotate.output.vcf , ref=REFERENCE_DIR , raw_vcf= rules.VariantAnnotator.output
  output: vcf=FILTER_OUTPUT+ "/"+SAMPLES[0] + ".final.annotated.vcf"
  shell: "java -Xmx4g -jar %s -T VariantAnnotator -R {input.ref} -A SnpEff --variant {input.raw_vcf} --snpEffFile {input.snpeff_vcf} -L {input.raw_vcf} -o {output.vcf}"\n""" %GATK_JAR

f.write(rule_all+"\n")
f.write(rule_make_picard_tmp+"\n")
f.write(rule_map_reads+"\n")
f.write(rule_sort_sam_convert_bam+"\n")
f.write(rule_dedup+"\n")
f.write(rule_index_bam+"\n")
f.write(rule_realignertargetcreator+"\n")
f.write(rule_indel_realigner+"\n")
f.write(rule_base_recalibration_first+"\n")
f.write(rule_base_recalibration_second+"\n")
f.write(rule_analyze_cov+"\n")
f.write(rule_apply_recal+"\n")
f.write(rule_haplotype_caller+"\n")
f.write(rule_VariantAnnotator+"\n")
f.write(rule_VariantRecalibrator+"\n")
f.write(rule_apply_recalibration+"\n")
f.write(rule_VariantRecalibrator_indel+"\n")
f.write(rule_apply_recalibration_indel+"\n")
#f.write(rule_annotate_with_dbsnp+"\n")
f.write(rule_snpEff_annotate+"\n")
f.write(rule_filter_impact_variants+"\n")
#f.write(rule_filter_moderate_impact_variants+"\n")
f.write(rule_annotate_impact_with_dbsnp+"\n")
#f.write(rule_annotate_moderate_impact_with_dbsnp+"\n")
f.write(rule_annotate_with_dbnsfp+"\n")
#f.write(rule_annotate_moderate_with_dbnsfp+"\n")
f.write(rule_extract_to_csv+"\n")
#f.write(rule_extract_moderate_to_csv+"\n")
#f.write(rule_post_sneff_gatk_annotate+"\n")
f.close()
# f.write()




# print (experimentName)



#print (SoftwareConfiguration.version)

