#!bin/bash

#Transposable elements annotation by EDTA 
EDTA.pl --genome asm.hap1.fasta --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 32 2> edta.log &

#filtering of the EDTA output using Rscript (find script on folder)

#Gene annotation using the funannotate tool

