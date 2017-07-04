__author__ = 'HOI QIANGZE'
import argparse
import os
import glob
import sys
import SoftwareConfiguration
import datetime
import re

# time=datetime.date.today()
# date=time.strftime("%y%m%d")
timenow=datetime.datetime.now().strftime("%y%m%d.%H%M")
version= str(SoftwareConfiguration.version)
SOFTWARE_FOLDER=str(SoftwareConfiguration.software_directory)
BWA_EXE=str(SoftwareConfiguration.bwa_exe)

#fetch options

#declare some variables
parser = argparse.ArgumentParser(description="Use the following options to run the preparation step for Next Seq Analysis")
parser.add_argument("--name" , '-n', help="Name of experiment. This will be used to create the folder")
parser.add_argument("--id" , '-i', help="Run Number #.", default=1)
parser.add_argument("--input" , '--in', help="Input folder containing fastq.")
parser.add_argument("--ref" , '--r', help="reference. eg hg19")

args = vars(parser.parse_args())

experimentName=args['name']
runNo=args['id']
reference=args['ref']
if not reference:
    reference='hg19'

homefolder='/home/ubuntu'
workingdir='/home/ubuntu/output' # for now use this
#workingdir=os.getcwd() # for now use this
experimentPath=workingdir+"/"+experimentName
rootpath=""
#glob all the fastq files in the given input folder
fastq_folder=args['input']
fastqfile=glob.glob(fastq_folder+"/*fastq.gz")

fastq_prefix = []
print (fastqfile)
for files in fastqfile:
    print("found this file %s from input folder %s" %(files,fastq_folder))
    m = re.match("(.*)\.fastq.gz",os.path.basename(files))
    if m:
        fastq_prefix.append(m.group(1))
        print (m.group(1))

if not os.path.exists(experimentPath):
    os.makedirs(experimentPath)

rootpath=experimentPath+"/Run."+ timenow+ "pipeline"+version
os.makedirs(experimentPath+"/Run."+ timenow+ "pipeline "+version)
os.makedirs(rootpath+"/INPUT_FILES")
os.makedirs(rootpath+"/CONFIG_FILES")
os.makedirs(rootpath+"/MAPPING")
os.makedirs(rootpath+"/GATK_PROCESSING")
os.makedirs(rootpath+"/FILTERING")
os.makedirs(rootpath+"/METRICS")
print("Root output folder at %s" %rootpath)
    # sys.exit("This folder exist. Please make use of another name or check %s" %experimentPath)
#here i copy the fastq to INPUT FILES FOLDER
INPUT=rootpath+"/INPUT_FILES"
#DEBUG: Commend out for debugging
#os.system('cp '+fastq_folder+'/*fastq.gz '+INPUT)
print ('cp '+fastq_folder+'/*fastq.gz '+INPUT) 
# here we write the make file
snakefile=rootpath+"/NextSeq.snakefile"
f=open(snakefile,"w")

#declare some variables. will be easier to write
SAMPLES='SAMPLES = \''+" ".join(fastq_prefix) + "\'\n"
ROOT_FOLDER= 'ROOT_FOLDER = \''+rootpath + "\'\n"
REF_NAME= "REF_NAME = \'%s\'\n" %reference
REFERENCE_DIR="REFERENCE_DIR = \'%s/%s/ucsc.hg19.fa\'\n" %(homefolder,reference)
SOFTWARE_FOLDER_PATH= 'SOFTWARE_FOLDER =  \''+SOFTWARE_FOLDER + "\'\n"
FASTQ_DIR= 'FASTQ_DIR = \''+INPUT+ "\'\n"
BWA_OUTPUT='BWA_OUTPUT = \''+rootpath+"/MAPPING\'\n"
GATK_OUTPUT='GATK_OUTPUT = \''+rootpath+"/GATK_PROCESSING\'\n"

f.write(SAMPLES)
f.write(ROOT_FOLDER)
f.write(REF_NAME)
f.write(REFERENCE_DIR)
f.write(SOFTWARE_FOLDER_PATH)
f.write(FASTQ_DIR)
f.write(BWA_OUTPUT)
f.write(GATK_OUTPUT)

rule_all="""rule all:
  input: RESULTS_DIR + PREFIX + REF_NAME + ".bwa.40.fas" \n"""

rule_map_reads="""rule map_reads :
  input: ref = REFERENCE_DIR, reads = FASTQ_DIR + "{sample}.fastq.gz"
  output: temp(BWA_OUTPUT + "{sample[0]}" + REF_NAME +  ".bam")
  shell: "%s mem -M -R'@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' -p {input.ref} {input.reads} > {output}\n""" %BWA_EXE

rule_sort_bam= """input: RESULTS_DIR + "{sample}_vs_" + REF_NAME + ".bam"
  output: protected(RESULTS_DIR + "{sample}_vs_" + REF_NAME + ".sorted.bam")
  run:
    shell("/usr/local/samtools/samtools sort {input} " + RESULTS_DIR + "{wildcards.sample}_vs_%s.sorted" % REF_NAME)\n"""

rule_index_bam = """rule index_bam:
  input: RESULTS_DIR + "{sample}_vs_" + REF_NAME + ".sorted.bam"
  output: protected(RESULTS_DIR + "{sample}_vs_" + REF_NAME + ".sorted.bam.bai")
  shell: "/usr/local/samtools/samtools index {input}" \n"""

rule_call_GATK= """input: ref = REF, bam = expand(RESULTS_DIR + "{sample}_vs_" + REF_NAME + ".sorted.bam", sample=SAMPLES), bam_indicies = expand(RESULTS_DIR + "{sample}_vs_" + REF_NAME + ".sorted.bam.bai", sample=SAMPLES)
  output: RESULTS_DIR + PREFIX + REF_NAME + ".bwa.vcf"
  run:
    inputs_with_prefix = ["-I " + bam_file for bam_file in expand(RESULTS_DIR + "{sample}_vs_" + REF_NAME + ".sorted.bam", sample=SAMPLES)]
    shell("java -Xmx30g -jar /usr/local/gatk2/GenomeAnalysisTK.jar -T UnifiedGenotyper -nt 12 -R  {input.ref} --sample_ploidy 2 --genotype_likelihoods_model BOTH -rf BadCigar {inputs_with_prefix} -o {output}") \n"""
f.write(rule_all)
f.write(rule_map_reads)

# f.write()




# print (experimentName)



#print (SoftwareConfiguration.version)

