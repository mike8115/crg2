rule EHdn:
    input: "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str/EHdn/{sample}-{unit}.str_profile.json",
        motif = "str/EHdn/{sample}-{unit}.locus.tsv",
        locus = "str/EHdn/{sample}-{unit}.motif.tsv",
        manifest = "str/EHdn/{project}_manifest.txt"
    params:
        ref = config["ref"]["genome"]
    log:
        "logs/str/{sample}-{unit}-EHdn.log"
    wrapper:
        get_wrapper_path("EHdn", "profile")   

rule EHdn_combine:
   input:
        json = expand("str/EHdn/{sample}-{unit}.json",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)
        manifest = "str/EHdn/{project}_manifest.txt"
    output:
        "str/EHdn/{project}_multisample_profile.json"
    params:
        ref = config["ref"]["genome"]
    log:
        "logs/str/{sample}-{unit}-EHdn.log"
    wrapper:
        get_wrapper_path("EHdn", "profile")   

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