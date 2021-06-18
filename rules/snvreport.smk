rule allsnvreport:
    input:
        db="annotated/{p}/{family}-gemini.db",
        vcf="annotated/{p}/vcfanno/{family}.{p}.vep.vcfanno.vcf"
    output:
        directory("report/{p}/{family}")
    conda:
        "../envs/cre.yaml"
    log:
        "logs/report/{p}/{family}.cre.log"
    resources:
         mem_mb=40000
    params:
         cre=config["tools"]["cre"],
         ref=config["ref"]["genome"]
    shell:
         '''
         mkdir -p {output}
         cd {output}
         ln -s ../../../{input.db} {project}-ensemble.db
         bgzip ../../../{input.vcf} -c > {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         tabix {project}-gatk-haplotype-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz {project}-ensemble-annotated-decomposed.vcf.gz
         ln -s {project}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi {project}-ensemble-annotated-decomposed.vcf.gz.tbi
         cd ../
         if [ {wildcards.p} == "coding" ]; then  
         {params.cre}/cre.sh {project} 
         elif [ {wildcards.p} == "denovo" ]; then  
         type=denovo {params.cre}/cre.sh {project} 
         else
         type=wgs {params.cre}/cre.sh {project}
         unset type
         fi;
         '''
if config["run"]["hpo"]:

    def get_panel(wildcards):
        if not config["run"]["panel"]:
            return "genes/genes.bed"
        return config["run"]["panel"]

    def get_bed(wildcards):
        if wildcards.p == "panel-flank":
            return "{}-flank-{}k.bed".format(config["run"]["project"], int(config["run"]["flank"]/1000))
        return get_panel()


    rule hpo_to_panel:
        input: 
            hpo=config["run"]["hpo"],
            ensembl=config["genes"]["ensembl"],
            refseq=config["genes"]["refseq"],
            hgnc=config["genes"]["hgnc"]
        params: 
            crg2=config["tools"]["crg2"],
            cre=config["tools"]["cre"]
        output: 
            genes="genes/genes.bed"
        conda: "../envs/hpo_to_panel.yaml"
        log:
            "logs/hpo_to_panel/genes.log"
        script:
            "../scripts/hpo_to_panel.py"

    rule add_flank:
        input: config["run"]["panel"]
        output: "{}-flank-{}k.bed".format(config["run"]["project"], int(config["run"]["flank"]/1000))
        params: config["run"]["flank"]
        shell:
            '''
            cat {input} | awk -F "\t" '{{print $1"\t"$2-{params}"\t"$3+{params}}}' | sed 's/-[0-9]*/0/g' | bedtools sort | bedtools merge > {output}
            '''

    rule intersect:
        input: 
            left="filtered/{family}.vcf.gz",
            right=get_bed
        output:
            vcf="filtered/{p}/{family}.{p}.vcf.gz"
        params:
            extra="-header"
        log: "logs/report/bedtools-{family}-{p}.log"
        wrapper:
            get_wrapper_path("bedtools", "intersect")
            
        