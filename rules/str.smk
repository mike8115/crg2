rule EHdn:
    input:
        reads = "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str/EHdn/{sample}-{unit}.str_profile.json",
        motif = "str/EHdn/{sample}-{unit}.locus.tsv",
        locus = "str/EHdn/{sample}-{unit}.motif.tsv"
    params:
        ref = config["ref"]["genome"]
    log:
        "logs/str/{sample}-{unit}-EHdn.log"
    wrapper:
        get_wrapper_path("EHdn")   

rule EH:
    input:
        reads = "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str/EH/{sample}-{unit}.json",
        vcf = "str/EH/{sample}-{unit}.vcf",
        bam = "str/EH/{sample}-{unit}_realigned.bam"
    params:
        ref = config["ref"]["genome"],
        sex = lambda w: "`sh {}/scripts/str_helper.sh mapped/{}-{}.sorted.bam`".format(workflow.basedir, w.sample, w.unit)
    log:
        "logs/str/{sample}-{unit}-EH.log"
    wrapper:
        get_wrapper_path("EH")