rule EHdn:
    input:
        reads = "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str_profiles/{sample}-{unit}.str_profile.json",
        motif = "str_profiles/{sample}-{unit}.locus.tsv",
        locus = "str_profiles/{sample}-{unit}.motif.tsv",
    params:
        ref = config["ref"]["genome"]
    log:
        "logs/str/{sample}-{unit}.log"
    wrapper:
        get_wrapper_path("EHdn")   