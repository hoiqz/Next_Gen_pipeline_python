SAMPLES = 'LB50_R1 LB50_R2'.split() 
SAMPLENAME = 'LB50_R1'
ROOT_FOLDER = '/mnt2/output/LB_study4/LB50_pipeline_0.01'
PICARD_TMP = '/mnt2/output/LB_study4/LB50_pipeline_0.01/tmp'
REF_NAME = 'hg19'
REFERENCE_DIR = '/mnt1/references/hg19/ucsc.hg19.fasta'
SOFTWARE_FOLDER =  '/home/ubuntu/installations'
FASTQ_DIR = '/mnt2/output/LB_study4/LB50_pipeline_0.01/INPUT_FILES'
BWA_OUTPUT = '/mnt2/output/LB_study4/LB50_pipeline_0.01/MAPPING'
GATK_OUTPUT = '/mnt2/output/LB_study4/LB50_pipeline_0.01/GATK_PROCESSING'
FILTER_OUTPUT = '/mnt2/output/LB_study4/LB50_pipeline_0.01/FILTERING'
KNOWN_INDELS='/mnt1/references/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'
DBSNP='/mnt1/references/dbsnp_138.hg19.vcf'
vcf_1000g_snp='/mnt1/references/1000G_phase1.snps.high_confidence.hg19.sites.vcf'
vcf_1000g_indel='/mnt1/references/1000G_phase1.indels.hg19.sites.vcf'
OMNI='/mnt1/references/1000G_omni2.5.hg19.sites.vcf'
HAPMAP='/mnt1/references/hapmap_3.3.hg19.sites.vcf'
rule all:
  input: FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.IMPACT.DP20.GENESET.DBSNP.DBNSFP.csv" 
 
rule make_picard_tmp :
  output: PICARD_TMP
  shell: "mkdir -p {output}"

rule map_reads :
  input: ref = REFERENCE_DIR, reads1 = FASTQ_DIR + "/"+SAMPLES[0]+".fastq.gz", reads2= FASTQ_DIR+"/"+SAMPLES[1]+".fastq.gz", dummy=rules.make_picard_tmp.output
  output: temp(BWA_OUTPUT + "/"+SAMPLES[0] + ".bam")
  shell: "/home/ubuntu/installations/bwa-0.7.12/bwa mem -M -R'@RG\\tID:group1\\tSM:{SAMPLENAME}\\tPL:illumina' {input.ref} {input.reads1} {input.reads2} > {output}"

rule sort_sam_convert_bam :
  input: mapped=rules.map_reads.output ,tmp=PICARD_TMP
  output: bam=temp(BWA_OUTPUT + "/"+SAMPLES[0] + ".sorted.bam")
  run:
    shell("java -Djava.io.tmpdir={input.tmp} -Xmx4g -jar ~/installations/picard-tools-1.130/picard.jar SortSam INPUT={input.mapped} OUTPUT={output.bam} SORT_ORDER=coordinate")

rule dedup:
  input: bamfile=rules.sort_sam_convert_bam.output.bam,tmp=PICARD_TMP
  output: dedup=temp(BWA_OUTPUT + "/"+SAMPLES[0] + ".sorted.dedup.bam"),metrics=BWA_OUTPUT+"/metrics"
  run:
    shell("java -Djava.io.tmpdir={input.tmp} -Xmx4g -jar ~/installations/picard-tools-1.130/picard.jar MarkDuplicates INPUT={input.bamfile} OUTPUT={output.dedup} METRICS_FILE={output.metrics}")

rule index_bam:
  input: dedup=rules.dedup.output.dedup,tmp=PICARD_TMP
  output: index=BWA_OUTPUT + "/"+SAMPLES[0] + ".sorted.dedup.bam.bai"
  run:
    shell("java -Djava.io.tmpdir={input.tmp} -Xmx4g -jar ~/installations/picard-tools-1.130/picard.jar BuildBamIndex INPUT={input.dedup} OUTPUT={output.index} ")

rule realignertargetcreator:
  input: ref=REFERENCE_DIR , known=KNOWN_INDELS, dummy=rules.index_bam.output.index
  output: GATK_OUTPUT +"/output.intervals"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {input.ref} -o {output} -known {input.known} "

rule indel_realigner:
  input: intervals=rules.realignertargetcreator.output , bam=rules.dedup.output.dedup , ref=REFERENCE_DIR , known=KNOWN_INDELS
  output: temp(GATK_OUTPUT+ "/"+SAMPLES[0] + ".sorted.dedup.realigned.bam")
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T IndelRealigner -I {input.bam} -R {input.ref} -targetIntervals {input.intervals} -o {output} -known {input.known} --consensusDeterminationModel KNOWNS_ONLY"

rule base_recalibration_first:
  input: bam=rules.indel_realigner.output , ref=REFERENCE_DIR , known=KNOWN_INDELS, known_sites=DBSNP
  output: GATK_OUTPUT+ "/recal_data.table"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T BaseRecalibrator -I {input.bam} -R {input.ref} -knownSites {input.known} -knownSites {input.known_sites} -o {output}"

rule base_recalibration_second:
  input: bam=rules.indel_realigner.output , ref=REFERENCE_DIR , known=KNOWN_INDELS, known_sites=DBSNP, bqsr= rules.base_recalibration_first.output
  output: GATK_OUTPUT+ "/post_recal_data.table"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T BaseRecalibrator -I {input.bam} -R {input.ref} -knownSites {input.known} -knownSites {input.known_sites} -BQSR {input.bqsr} -o {output}"

rule analyze_cov :
  input: before=rules.base_recalibration_first.output, after = rules.base_recalibration_second.output , ref=REFERENCE_DIR 
  output: GATK_OUTPUT+ "/recalibration_plots.pdf"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T AnalyzeCovariates -R {input.ref} -before {input.before} -after {input.after} -plots {output}"

rule apply_recal:
  input: bqsr=rules.base_recalibration_first.output, bam=rules.indel_realigner.output , ref=REFERENCE_DIR, dependency=rules.base_recalibration_second.output 
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + ".sorted.dedup.realigned.recal.bam" 
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T PrintReads -R {input.ref} -I {input.bam} -BQSR {input.bqsr} -o {output}"

rule haplotype_caller:
  input: bam=rules.apply_recal.output , ref=REFERENCE_DIR
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + "raw_variants.vcf" 
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T HaplotypeCaller -R {input.ref} -I {input.bam} -variant_index_type LINEAR -variant_index_parameter 128000 --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF -o {output}"

rule VariantAnnotator:
  input: vcf=rules.haplotype_caller.output, ref=REFERENCE_DIR
  output: GATK_OUTPUT+ "/"+SAMPLES[0] + "raw_variants.annotated.vcf" 
  shell: "java -Xmx4g -jar /home/ubuntu/installations/gatk-3.2.2/GenomeAnalysisTK.jar -T VariantAnnotator -R {input.ref} --variant {input.vcf} -A StrandOddsRatio -A FisherStrand -o {output}"

rule snpEff_annotate:
  input: vcf=rules.VariantAnnotator.output
  output: vcf=FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.vcf",stats=FILTER_OUTPUT+"/"+SAMPLES[0]+".stats"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/snpEff/snpEff.jar -c /home/ubuntu/installations/snpEff/snpEff.config  -stats {output.stats} hg19 {input.vcf} > {output.vcf}"

rule filter_impact_variants:
  input: vcf=rules.snpEff_annotate.output.vcf
  output: vcf=FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.IMPACTFUL.DP20.GENESET.vcf"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/snpEff/SnpSift.jar filter --set /mnt1/references/gene_list/obesity_genes.txt \"((ANN[*].IMPACT =  'MODERATE')|(ANN[*].IMPACT = 'HIGH'))  &(DP>=20) & (ANN[*].GENE in SET[0])\" -f {input.vcf} >{output.vcf}"

rule annotate_impact_with_dbsnp:
  input: vcf=rules.filter_impact_variants.output, dbsnp=DBSNP
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] + ".SNPEFF.IMPACTFUL.DP20.GENESET.DBSNP.vcf"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/snpEff/SnpSift.jar annotate -id {input.dbsnp} {input.vcf} > {output}"

rule annotate_with_dbnsfp:
  input: vcf=rules.annotate_impact_with_dbsnp.output
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] +  ".SNPEFF.IMPACT.DP20.GENESET.DBSNP.DBNSFP.vcf"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/snpEff/SnpSift.jar dbnsfp -v {input.vcf} > {output}"

rule extract_to_csv:
  input: vcf=rules.annotate_with_dbnsfp.output
  output: FILTER_OUTPUT+ "/"+SAMPLES[0] +  ".SNPEFF.IMPACT.DP20.GENESET.DBSNP.DBNSFP.csv"
  shell: "java -Xmx4g -jar /home/ubuntu/installations/snpEff/SnpSift.jar extractFields -s ',' -e '.' {input.vcf} 'CHROM' 'POS' 'REF' 'ALT' 'ID' 'DP' 'EFF[*].IMPACT' 'EFF[*].GENE' 'dbNSFP_SIFT_pred' 'dbNSFP_Polyphen2_HVAR_pred' 'dbNSFP_MutationTaster_pred' 'dbNSFP_1000Gp1_ASN_AF' > {output}"

