
rule input_prep:
    input:
        units=config["run"]["units"]
    params:
        outdir = temp("fastq/"),
        sort_check = True,
        old_cram_ref = config["ref"]["old_cram_ref"],
        new_cram_ref = config["ref"]["new_cram_ref"]
    output:
        fastq1 = temp("fastq/{family}_{sample}_R1.fastq.gz"),
        fastq2 = temp("fastq/{family}_{sample}_R2.fastq.gz")
    log:
        "logs/input_prep/{family}_{sample}.log"
    threads:
        4
    wrapper:
        get_wrapper_path("samtools", "fastq")
    
        
rule map_reads:
    input:
        reads = ["fastq/{family}_{sample}_R1.fastq.gz", "fastq/{family}_{sample}_R2.fastq.gz"]
    output:
        temp("mapped/{family}_{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{family}_{sample}.log"
    params:
        index = config["ref"]["genome"],
        extra = get_read_group,
        verbosity = config["params"]["bwa"]["verbosity"],
        maxMem = config["params"]["bwa"]["maxMem"],
        markSplitReads = config["params"]["bwa"]["markSplitReads"],
        sort = "samtools",
        sort_order = "coordinate"
    threads: 16
    resources: 
        mem=lambda wildcards, threads: threads * 2
    wrapper:
        get_wrapper_path("bwa", "mem")


rule mark_duplicates:
    input:
        "mapped/{family}_{sample}.sorted.bam"
    output:
        bam = temp("dedup/{family}_{sample}.bam"),
        metrics = "qc/dedup/{family}_{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{family}_{sample}.log"
    params:
        markDuplicates = config["params"]["picard"]["MarkDuplicates"],
        java_opts = config["params"]["picard"]["java_opts"],
    wrapper:
        get_wrapper_path("picard", "markduplicates")


rule recalibrate_base_qualities:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"],
        bed = "mapped/{family}_{sample}-callable.bed" 
    output:
        bam = protected("recal/{family}_{sample}.bam")
    params:
        extra = get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk/bqsr/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("gatk", "baserecalibrator")

rule realignertargetcreator:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"]
    output:
        intervals="recal/gatk3/realignertargetcreator/{family}_{sample}.intervals",
    log:
        "logs/gatk3/realignertargetcreator/{family}_{sample}.log"
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["RealignerTargetCreator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    threads: 8
    #resources: #cannot access threads here; fix later
     #   mem=lambda wildcards, threads: {threads} * 3
     #using already installed haplotypecaller conda-env, otherwise conda takes forever; 
     #todo: change to a common gatk3.yaml instead of seperate conda for sub-commands
    conda: 
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "realignertargetcreator")

rule indelrealigner:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"],
        target_intervals="recal/gatk3/realignertargetcreator/{family}_{sample}.intervals",
    output:
        bam = protected("recal/gatk3/indelrealigner/{family}_{sample}-realign.bam"),
    log:
        "logs/gatk3/indelrealigner/{family}_{sample}.log"
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["IndelRealigner"],
        java_opts = config["params"]["gatk"]["java_opts"],
    conda:
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "indelrealigner")

rule mosdepth:
    input:
        bam = get_recal_input(),
        bai = get_recal_input(bai=True),
    output:
        qbed = "mapped/{family}_{sample}.quantized.bed.gz",
        bed = "mapped/{family}_{sample}-callable.bed"
    params:
        prefix = f"mapped/{{family}}_{{sample}}",
        by = config["ref"]["canon_bed"],
        quantiles = "0:1:4:"
    log:
        "logs/mosdepth/{family}_{sample}.log"
    wrapper:
        get_wrapper_path("mosdepth")
       
rule baserecalibrator:
    input:
        bam = "recal/gatk3/indelrealigner/{family}_{sample}-realign.bam",
        bai = "recal/gatk3/indelrealigner/{family}_{sample}-realign.bam.bai",
        bed = "mapped/{family}_{sample}-callable.bed",
        ref = config["ref"]["genome"],
        known = config["ref"]["known-variants"]
    output:
        "recal/gatk3/baserecalibrator/{family}_{sample}-recal.grp"
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["BaseRecalibrator"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk3/baserecalibrator/{family}_{sample}.log"
    threads: 8
    conda:
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "baserecalibrator")

rule printreads:
    input:
        bam = "recal/gatk3/indelrealigner/{family}_{sample}-realign.bam",
        bai = "recal/gatk3/indelrealigner/{family}_{sample}-realign.bam.bai",
        ref = config["ref"]["genome"],
        recal_data = "recal/gatk3/baserecalibrator/{family}_{sample}-recal.grp"
    output:
        protected("recal/gatk3/{family}_{sample}.bam")
    params:
        extra = get_regions_param() + config["params"]["gatk3"]["PrintReads"],
        java_opts = config["params"]["gatk"]["java_opts"],
    log:
        "logs/gatk3/printreads/{family}_{sample}.log"
    conda:
        "../wrappers/gatk3/haplotypecaller/environment.yaml"
    wrapper:
        get_wrapper_path("gatk3", "printreads")
rule md5:
    input: 
        bam = "recal/{family}_{sample}.bam"
    output:
        md5 = "recal/{family}_{sample}.bam.md5"
    shell:
        """
        md5sum {input.bam} > {output.md5}
        """

rule remove_decoy:
    #redirect first samtools command to a temp output file "decoy.bam" 
    #otherwise it floods the log file with binary stream of decoy reads
    input:
        bam = "recal/{family}_{sample}.bam",
        canon = config["ref"]["canon_bed"],
    output:
        out_f = temp("decoy_rm/{family}_{sample}.no_decoy_reads.bam")
    log:
        "logs/remove_decoys/{family}_{sample}.log"
    threads: 8
    resources: mem_mb = 10000
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view --threads {threads} -h -L {input.canon} {input.bam} | egrep -v "hs37d5|NC_007605" | samtools view --threads {threads} - -bh > {output.out_f}
        """

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        get_wrapper_path("samtools", "index")


rule concatenate_callable_bed:
    input:
        expand("mapped/{family}_{sample}-callable.bed",family=project,sample=samples.index)
    output: 
        "mapped/{family}-concat-sort.bed"
    log:
        "logs/bash/{family}-callable-concat.log"
    shell:
        '''
            cat {input} | sort -k1,1 -k2,2n > {output} 2>{log}
        '''
    
rule merge_bed:
    input: 
        "mapped/{family}-concat-sort.bed"
    output:
        "mapped/{family}-sort-callable.bed"
    log:
        "logs/bedtools/{family}-callable-merge.log"
    wrapper:
        get_wrapper_path("bedtools", "merge")

rule contigwise_bed:
    input:
        "mapped/{family}-sort-callable.bed"
    output:
        "mapped/bed/{family}-sort-callable-{contig}.bed"
    log:
        "logs/bash/{family}.{contig}.log"    
    shell:
        """
            awk '{{ if($1=="{wildcards.contig}") print $0; }}' {input} > {output} 2>{log}
        """
    
