from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
prefix = os.path.splitext(snakemake.output.json)[0].split(".")

shell(
    "(ExpansionHunterDenovo profile --reference {snakemake.params.ref} --reads {snakemake.input.reads} --output-prefix {prefix} --min-anchor-mapq 50 --max-irr-mapq 40 ) {log}"
)

