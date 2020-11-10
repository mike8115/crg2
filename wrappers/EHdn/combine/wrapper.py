from snakemake.shell import shell
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
prefix = os.path.splitext(snakemake.output.json)[0].split(".")


shell(
    "( ExpansionHunterDenovo merge"
    " --reference {snakemake.params.ref}"
    " --manifest {manifest}"
    " --output-prefix {prefix} )"
    " {log}"
)
