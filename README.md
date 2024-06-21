# simple_vcf_comparison

A Simple way to compare differences between vcfs with human readable output. (Done for a work project, not intended to be useful outside this scope)

To execute the process, rename the ILLUMINA files to match the NANOPORE files:
- If the used technology is DEVYSER, leave them just as "BRCA_xx_year"
- If the used technology is HC SOPHIA, rename them as "BRCA_SOPHIA_xx_year"

Separate the ILLUMINA and NANOPORE vcfs in 2 separate folders (ILLUMINA and NANOPORE, or if given different names, change the parameter in the first script)
Make sure that every ILLUMINA file name as an analogous in the NANOPORE folder.
