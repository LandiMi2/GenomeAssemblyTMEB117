## Haplotype-resolved genome of heterozygous African cassava TMEB117 _(Manihot esculenta)_



This repository contains scripts of the methodology for assembling the genome of the TMEB117 cassava cultivar. Utilizing PacBio HiFi reads, we constructed a haplotype-resolved genome with an accuracy of QV > 64, N50 > 35 Mbp, and 98.9% BUSCO completeness. This represents the most precise cassava genome sequence to date. The genome is characterized by heterozygosity, encompassing over 60% repetitive elements. Our gene prediction identified more than 45,000 gene models for both haplotypes.

#### Scripts for methodology section
- BashScripts
  - QC_to_Assembly.sh: Bash commands covering quality control to genome assembly.
  - annotation.sh: Commands for transposable elements (TE) and the funannotate singularity script used in the gene annotation section.
- R scripts: scripts for TE and gene annotation visualizations

For more details and methods, please take a look at the published paper [here](https://www.nature.com/articles/s41597-023-02800-0).







