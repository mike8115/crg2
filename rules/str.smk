rule EHdn_profile:
    input: 
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str/EHdn/{sample}-{unit}.str_profile.json",
        motif = "str/EHdn/{sample}-{unit}.locus.tsv",
        locus = "str/EHdn/{sample}-{unit}.motif.tsv",
    params:
        ref = config["ref"]["genome"]
    log:
        "logs/str/{sample}-{unit}-EHdn.log"
    wrapper:
        get_wrapper_path("EHdn", "profile")   

def get_manifest(wildcards):
    json = expand("str/EHdn/{sample}-{unit}.str_profile.json", sample=wildcards.sample, unit=units.loc[wildcards.sample].unit)
    with open("str/EHdn/manifest.txt", "w") as f:
        for i,j in zip(wildcards.sample, json):
            f.writelines("{}\tcase\t{}\n".format(i,j))



rule EHdn_merge:
    input:
        #json = expand("str/EHdn/{sample}-{unit}.json",sample=wildcards.sample, unit=units.loc[wildcards.sample].unit),
        manifest = get_manifest
    output:
        "str/EHdn/all/{project}_multisample_profile.json"
    params:
        ref = config["ref"]["genome"]
    log:
        "logs/str/{project}-EHdn.log"
    wrapper:
        get_wrapper_path("EHdn", "merge")   

rule EH:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        json = "str/EH/{sample}-{unit}.json",
        vcf = "str/EH/{sample}-{unit}.vcf",
        bam = "str/EH/{sample}-{unit}_realigned.bam"
    params:
        ref = config["ref"]["genome"],
        sex = lambda w: "`sh {}/scripts/str_helper.sh mapped/{}-{}.sorted.bam`".format(workflow.basedir, w.sample, w.unit),
        catalog = config["annotation"]["eh"]["catalog"]
    log:
        "logs/str/{sample}-{unit}-EH.log"
    wrapper:
        get_wrapper_path("EH")

rule EH_report:
    input:
        json = get_eh_json()
    output:
        tsv = "str/EH/{project}_EH_str.tsv",
        annot = "str/EH/{project}_EH_str.annot.tsv",
        xlsx = "str/EH/{project}_EH_v1.1.xlsx"
    log:
        "logs/str/{project}-eh-report.log"
    params:
        trf = config["annotation"]["eh"]["trf"]
    conda:
        "../envs/eh-report.yaml"
    shell:
        '''
        python ../scripts/generate_EH_genotype_table.generic.py str/EH > {output}
        python ../scripts/add_gene+threshold_to_EH_column_headings2.py {output.tsv} {params.trf} > {output.annot}
        python ../scripts/eh_sample_report.py {output.trf} {output.xlsx}

        '''

