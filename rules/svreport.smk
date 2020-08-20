rule svreport:
    input:
        get_annotated_sv_vcf()
    output:
        dir=directory("report/sv")
    log:
        "logs/report/sv.log"
    params:
        PIPELINE_VERSION=PIPELINE_VERSION,
        project=config["run"]["project"],
        hgmd=config["annotation"]["svreport"]["hgmd"],
        hpo=config["run"]["hpo"],
        protein_coding_genes=config["annotation"]["svreport"]["protein_coding_genes"],
        exon_bed=config["annotation"]["svreport"]["exon_bed"],
        exac=config["annotation"]["svreport"]["exac"],
        omim=config["annotation"]["svreport"]["omim"],
        gnomad=config["annotation"]["svreport"]["gnomad"],
        biomart=config["annotation"]["svreport"]["biomart"],
        mssng_manta_counts=config["annotation"]["svreport"]["mssng_manta_counts"],
        mssng_lumpy_counts=config["annotation"]["svreport"]["mssng_lumpy_counts"],
    conda:
        "../envs/crg.yaml"
    script:
        os.path.join(config["tools"]["crg"], "crg.intersect_sv_vcfs.py")
