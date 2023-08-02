#!bin/bash

#Transposable elements annotation by EDTA 
EDTA.pl --genome asm.hap1.fasta --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 32 2> edta.log &

#filtering of the EDTA output using Rscript (find the script in the Rscript folder)

#Gene annotation using the funannotate tool
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -w compute06
#SBATCH -J fun-annotate
#SBATCH -n 8
#SBATCH -o /home/mkofia/TMEB117Anno/slurm_out
#SBATCH -e /home/mkofia/TMEB117Anno/slurm_errors

workdir=/home/mkofia/TMEB117Anno/fun
genome=/home/mkofia/TMEB117Anno/genome
protein=/home/mkofia/TMEB117Anno/protein_evidence
trinity=/home/mkofia/TMEB117Anno/rna_asm
rna=/home/mkofia/TMEB117Anno/rna_seq
anno=/home/mkofia/TMEB117Anno/anno

#module load funannotate/1.8.9 
#run train 
singularity run \
	--bind "$workdir:$workdir","$genome:$genome","$trinity:$trinity","$rna:$rna" \
	/export/apps/funannotate/1.8.11/funannotate.sif \
	funannotate train -i "${genome}/asm.hap1.fasta" -o "$workdir" \
  --left "${rna}/R1.fq.gz" \
  --right "${rna}/R2.fq.gz" \
  --species "Manihot esculenta" \
  --cpus 8 --max_introlen 10000 --jaccard_clip \
  --trinity "${trinity}/Trinity.fasta"

#run predict
singularity run \
	--bind "$workdir:$workdir","$genome:$genome","$protein:$protein" \
	/export/apps/funannotate/1.8.11/funannotate.sif \
	funannotate predict -i "${genome}/asm.hap1.fasta" -o "$workdir" \
  --species "Manihot esculenta" --cpus 8 --max_introlen 10000 \
  --busco_db embryophyta --repeats2evm --organism other \
  --optimize_augustus --protein_evidence "${protein}/UniprotRevised.fasta"

#run update
singularity run \
	--bind "$workdir:$workdir" \ 
  /export/apps/funannotate/1.8.11/funannotate.sif \
  funannotate update -i "$workdir" --cpus 8 --max_intronlen 10000

#run eggnog separately to use in the last step of annotate 
emapper.py -i Manihot_esculenta.proteins.fa --output eggnog.emapper -m diamond --cpu 10

#run annotate 
singularity run \
	--bind "$workdir:$workdir","$anno:$anno" \
	/export/apps/funannotate/1.8.11/funannotate.sif \
	funannotate annotate -i "${workdir}" --cpus 8 \
	--eggnog "${anno}/eggnog.emapper.annotations" \
	--busco_db embryophyta \
	--species "Manihot esculenta" \
	--out "${workdir}/annotate"

