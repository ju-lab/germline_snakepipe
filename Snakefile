# Snakefile to run FASTQ -> Processed BAM -> basic germline calling usiung `Freebayes`
# Reused Snakefile from somatic pipeline
# 2018.06.14 Jongsoo Yoon

configfile: 'pathConfig.yaml'
configfile: 'sampleConfig.yaml'

# import pysam

# # define major contigs in the reference to parallelize Freebayes per chromosome
# ref = pysam.FastaFile(config['reference'])
# major_contigs = [str(i) for i in range(1, 23)] + ['chr' + str(i) for i in range(1, 23)] + \
#                  ['X', 'chrX', 'Y', 'chrY', 'MT', 'chrM']

# chromosomes = [chrom for chrom in ref.references if chrom in major_contigs]

rule all:
    input:
        'done'

rule bwa_align:
    input:
        bwa = config['bwa'],
        samtools = config['samtools'],
        ref = config['reference'],
        fq1 = lambda wildcards: config['samples'][wildcards.sample]['fq1'],
        fq2 = lambda wildcards: config['samples'][wildcards.sample]['fq2']

    output:
        bam = temp("temp_dna/{sample}.temp.bam"),
        sortedbam = temp('temp_dna/{sample}.temp.sorted.bam')
    params:
        rg = "@RG\tID:{sample}\tSM:{sample}\tPL:Illumina"
    log:
        "logs/{sample}.bwa.log"
    threads: 4
    shell:
        "({input.bwa} mem -t {threads} -R '{params.rg}' {input.ref} {input.fq1} {input.fq2} |"
        "{input.samtools} view -Sb - > {output.bam}; {input.samtools} sort -@ {threads} -o {output.sortedbam} {output.bam}) "
        " &> {log}"



rule markdup:
    input:
        java = config['java8'], 
        picard = config['picard'],
        samtools = config['samtools'],
        sortedbam = "temp_dna/{sample}.temp.sorted.bam", 
    output:
        mdbam = temp("temp_dna/{sample}.temp.sorted.md.bam"),
        mdbai = temp("temp_dna/{sample}.temp.sorted.md.bam.bai")

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}.md.log"
    shell:
        "({input.java} -XX:ParallelGCThreads={threads} -Xmx8g -jar {input.picard} MarkDuplicates "
        "REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true I={input.sortedbam} O={output.mdbam} "
        "M={output.mdbam}.metric VALIDATION_STRINGENCY=LENIENT TMP_DIR=md_temp QUIET=true; "
        "{input.samtools} index {output.mdbam})"
        " &> {log}"

rule realign:
    input:
        java = config['java8'],
        gatk = config['gatk'],
        samtools = config['samtools'],
        ref = config['reference'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}.temp.sorted.md.bam", 

    output:
        realignedbam = temp("temp_dna/{sample}.temp.sorted.md.ir.bam"),
        realignedbai = temp("temp_dna/{sample}.temp.sorted.md.ir.bam.bai")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}.realign.log"
    shell:
        "({input.java} -Xmx4g -jar {input.gatk} -T RealignerTargetCreator -R {input.ref} "
        " -I {input.bam} --known {input.knownindel} -o {output.realignedbam}.intervals; "
        "{input.java} -Xmx4g -jar {input.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} "
        "-targetIntervals {output.realignedbam}.intervals -o {output.realignedbam}; "
        "{input.samtools} index {output.realignedbam}) "
        " &> {log}"

rule baserecal:
    input:
        java = config['java8'], 
        gatk = config['gatk'], 
        samtools = config['samtools'],
        ref = config['reference'], 
        dbsnp = config['dbsnp'], 
        knownindel = config['knownindel'], 
        bam = "temp_dna/{sample}.temp.sorted.md.ir.bam"
    output:
        recalTable = temp('temp_dna/{sample}.recaltable'),
        recalbam = "dna_bam/{sample}.bam"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "logs/{sample}.baserecal.log"    
    shell:
        "({input.java} -Xmx4g -jar {input.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -knownSites {input.dbsnp} --knownSites {input.knownindel} -o {output.recalTable}; "
        "{input.java} -Xmx4g -jar {input.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {output.recalTable} -o {output.recalbam} -nct {threads}; "
        "{input.samtools} index {output.recalbam})"
        " &> {log}"

rule freebayes_parallel:
    input:
        inputbams=expand("dna_bam/{sample}.bam", sample=config["samples"]), 
        freebayes_parallel = config['freebayes_parallel'], 
        ref = config['reference'], 
        fasta_generate_regions = config['fasta_generate_regions']

    output:
        combined_vcf = "variants/everyone.freebayes.vcf",
    threads: 12
    log:
        "logs/freebayes.log"
    run:
        input_string = ['-b' + bam for bam in input.inputbams]
        input_string = " ".join(input_string)
        shell("({input.freebayes_parallel} <({input.fasta_generate_regions} {input.ref} 100000) {threads} -f {input.ref} \
        --standard-filters -j --min-coverage 10 -F 0.2 -C 2 \
        --read-snp-limit 3 --read-mismatch-limit 3 --ploidy 2 {inputstring} > {output.comvined_vcf}) \
        &> {log}")


rule vt_decompose_normalize:
    # decomposes multiallelic primitives into SNPs and indels
    input:
        combined_vcf = "variants/everyone.freebayes.vcf", 
        ref = config['reference'], 
        vt = config['vt']
    params:
        project = config['project']
    log:
        "logs/{params.project}.vt.log"
    output:
        vt_vcf = "variants/{params.project}.everyone.freebayes.vt.vcf"
    shell:
        "({input.vt} decompose -s {input.combined_vcf} | {input.vt} normalize -r {input.ref} - > {output.vt_vcf}) "
        "&> {log}"


rule finish:
    input:
        final_vcf = "variants/{params.project}.everyone.freebayes.vt.vcf"
    params:
        project = config['project']
    output:
        "done"
    shell:
        "touch done"

