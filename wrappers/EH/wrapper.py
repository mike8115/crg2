from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
prefix = os.path.splitext(snakemake.output.json)[0].split(".")

# must be moved to config
catalog = "/hpf/largeprojects/ccmbio/arun/C4Rare/Tandem_repeat_disease_loci.hg19.json"

shell(
    "sex={snakemake.params.sex} && "
    "( ExpansionHunter"
    " --reference {snakemake.params.ref}"
    " --reads {snakemake.input.reads}"
    " --output-prefix {prefix}"
    " --variant-catalog {catalog}"
    " --sex $sex )"
    " {log}"
)

