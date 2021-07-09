#gatk VCF is missing because it is symlinked to ensemble VCF
callers = [ "samtools", "freebayes", "platypus" ] 

rule allsnvreport:
    input:
        db="annotated/coding/{family}-gemini.db",
        vcf="annotated/coding/vcfanno/{family}-ensemble-decomposed.vep.vcfanno.vcf",
        caller_vcfs = expand("filtered/{family}-{caller}-vt-sf-pass.vcf.gz", family=project, caller=callers)
    output:
        directory("report/coding/{family}")
    conda:
        "../../envs/cre.yaml"
    log:
        "logs/report/{family}/cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         ref=config["ref"]["genome"],
         family = config["run"]["project"]
    shell:
         '''
         mkdir -p {output}
         cd {output}
         ln -s ../../../{input.db} {params.family}-ensemble.db
         for i in {input.caller_vcfs}; do
            caller=`basename ${{i}} .vcf.gz | cut -d "-" -f2`;
            ln -s ../../../${{i}} {params.family}-${{caller}}-annotated-decomposed.vcf.gz;
         done
         bgzip ../../../{input.vcf} -c > {params.family}-gatk-haplotype-annotated-decomposed.vcf.gz
         tabix {params.family}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s {params.family}-gatk-haplotype-annotated-decomposed.vcf.gz {params.family}-ensemble-annotated-decomposed.vcf.gz
         ln -s {params.family}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi {params.family}-ensemble-annotated-decomposed.vcf.gz.tbi
         cd ../
         {params.cre}/cre.sh {params.family} 
         type=wes.synonymous {params.cre}/cre.sh {params.family}
         unset type
         '''


