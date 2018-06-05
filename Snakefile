# Snakefile to run FASTQ -> Processed BAM -> basic somatic variant calling (SNV, SV, and CNV)
# 2018.06.03 Jongsoo Yoon

configfile: 'pathConfig.yaml'
configfile: 'sampleConfig.yaml'

rule all:
    input:
        'done'


rule bwa_normal:
    input:
        bwa = config['bwa'],
        samtools = config['samtools'],
        ref = config['reference'],
        fq1 = lambda wildcards: config['samples'][wildcards.sample]['normal_fastq'] + '1.fastq.gz',
        fq2 = lambda wildcards: config['samples'][wildcards.sample]['normal_fastq'] + '2.fastq.gz'

    output:
        bam = temp("temp_dna/{sample}_N.temp.bam"),
        sortedbam = temp('temp_dna/{sample}_N.temp.sorted.bam')
    params:
        rg = "@RG\tID:{sample}_N\tSM:{sample}_N\tPL:Illumina"
    log:
        "logs/{sample}_N.bwa.log"
    threads: 4
    shell:
        "{input.bwa} mem -t {threads} -R '{params.rg}' {input.ref} {input.fq1} {input.fq2} |"
        "{input.samtools} view -Sb - > {output.bam}; {input.samtools} sort -@ {threads} -o {output.sortedbam} {output.bam}"



rule markdup:
    input:
        java = config['java8'], 
        picard = config['picard'],
        samtools = config['samtools'],
        sortedbam = "temp_dna/{sample}_{tn}.temp.sorted.bam", 
    output:
        mdbam = temp("temp_dna/{sample}_{tn}.temp.sorted.md.bam"),
        mdbai = temp("temp_dna/{sample}_{tn}.temp.sorted.md.bam.bai")

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}_{tn}.md.log"
    shell:
        "{input.java} -XX:ParallelGCThreads={threads} -Xmx8g -jar {input.picard} MarkDuplicates "
        "REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true I={input.sortedbam} O={output.mdbam} "
        "M={output.mdbam}.metric VALIDATION_STRINGENCY=LENIENT TMP_DIR=md_temp QUIET=true; "
        "{input.samtools} index {output.mdbam}"

rule realign:
    input:
        java = config['java8'],
        gatk = config['gatk'],
        samtools = config['samtools'],
        ref = config['reference'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}_{tn}.temp.sorted.md.bam", 

    output:
        realignedbam = temp("temp_dna/{sample}_{tn}.temp.sorted.md.ir.bam"),
        realignedbai = temp("temp_dna/{sample}_{tn}.temp.sorted.md.ir.bam.bai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    shell:
        "{input.java} -Xmx4g -jar {input.gatk} -T RealignerTargetCreator -R {input.ref} "
        " -I {input.bam} --known {input.knownindel} -o {output.realignedbam}.intervals; "
        "{input.java} -Xmx4g -jar {input.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} "
        "-targetIntervals {output.realignedbam}.intervals -o {output.realignedbam}; "
        "{input.samtools} index {output.realignedbam}"


rule baserecal:
    input:
        java = config['java8'], 
        gatk = config['gatk'], 
        samtools = config['samtools'],
        ref = config['reference'], 
        dbsnp = config['dbsnp'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}_{tn}.temp.sorted.md.ir.bam"
    output:
        recalTable = temp('temp_dna/{sample}_{tn}.recaltable'),
        recalbam = "dna_bam/{sample}_{tn}.bam"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    shell:
        "{input.java} -Xmx4g -jar {input.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -knownSites {input.dbsnp} --knownSites {input.knownindel} -o {output.recalTable}; "
        "{input.java} -Xmx4g -jar {input.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {output.recalTable} -o {output.recalbam} -nct {threads}; "
        "{input.samtools} index {output.recalbam}"

rule freebayes:
    input:
        inputbams=expand("dna_bam/{sample}_{tn}.bam", sample=config["samples"], tn=['N']), 
        freebayes = config['freebayes'], 
        ref = config['reference'], 
    params:
        project = config['project']
    output:
        combined_vcf = "variants/{params.project}.everyone.freebayes.vcf"
    threads: 1

    run:
        input_string = ['-b' + bam for bam in input.inputbams]
        input_string = " ".join(input_string)
        shell("{input.freebayes} -f {input.ref} -v {output.combined_vcf} \
        --standard-filters -j --min-coverage 10 -F 0.2 -C 2 \
        --read-snp-limit 3 --read-mismatch-limit 3 --ploidy 2 {inputstring}")

rule decompose:
    # decomposes multiallelic primitives into SNPs and indels
    input:
        combined_vcf = "variants/{params.project}.everyone.freebayes.vcf", 
        vcfallelicprimitives = config['vcfallelicprimitives'], 

    params:
        project = config['project']
    output:
        decomposed_vcf = "variants/{params.project}.everyone.freebayes.dc.vcf"
    shell:
        "{input.vcfallelicprimitives} -kg {input.combined_vcf} > {output.decomposed_vcf}"

rule leftnorm:
    input:
        decomposed_vcf = "variants/{params.project}.everyone.freebayes.dc.vcf", 
        ref = config['reference'], 
        vt = config['vt'], 
    params:
        project = config['project']
    output:
        normalized_vcf = "variants/{params.project}.everyone.freebayes.dc.norm.vcf"
    shell:
        "{input.vt} normalize -r {input.ref} {input.decomposed_vcf} > {output.normalized_vcf}"


rule finish:
    input:
        final_vcf = "variants/{params.project}.everyone.freebayes.dc.norm.vcf"
    params:
        project = config['project']
    output:
        "done"
    shell:
        "touch done"

