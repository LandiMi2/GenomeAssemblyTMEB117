#!bin/bash

#Note these are commands used without specific files paths

#Fastqc version - 0.11.5
fastqc -o QC <raw_reads>.fastq.gz

#Filtering reads of length 10k to 30K using fastp

fastp.0.23.1 -i \
    <raw_reads>.fastq.gz \
    -o filtered<raw_reads>.fastq.gz --length_required 10000 --length_limit 30000

#mapping HiFi raw reads to a genome using pbmm2 aligner (https://github.com/PacificBiosciences/pbmm2)
#Pbmm2 version - 1.10.0; minimap2 2.15
Pbmm2 index <ref.fa> <ref.mmi>
pbmm2 align ref.index/<ref>.mmi <raw_reads>.fastq.gz <output>.align.bam --unmapped --log-level INFO --log-file align.log -j 10

#extracting unmapped reads 
samtools fastq -f 4 <output>.align.bam > unmapped.fastq #-f (unmapped) -F (mapped) reads

###### Assembly using hifiasm version - 0.16.1-r375
hifiasm -o asm -t 32 filtered<raw_reads>.fastq.gz 2> hifiasm.cell1.log

#HiCanu Version - 2.3
canu -d hicanu.asm -p TME117 -pacbio-hifi filtered<raw_reads>.fastq.gz  genomeSize=750m -useGrid=false -merylThreads=4 -merylMemory=8 corOverlapper=ovl 2> hicanu.log

##Flye Version - 2.9-b1768
flye --pacbio-hifi filtered<raw_reads>.fastq.gz -o flyasm --genome-size 750m --threads 32  2> flye.log

########### Downstream with hifiasm assembly 
#get primary contigs in FASTA from Hifiasm assembly gfa file 
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa 

#Run stats - n50 
n50 -x asm.hap1.fa asm.hap2.fa > assembly.stats

##Quast version - 5.1.0rc1
quast.py asm.hap1.fa asm.hap2.fa -o quast.asm.hap1.and.2.report

#Busco version - 5.3.2 ;  lineage = eudicots_odb10
busco  -m genome -i asm.hap1.fa -o busco.asm.hap1 -l eudicots_odb10 -c 10

#Evaluating the genome using merqury 
git clone https://github.com/marbl/merqury.git
export MERQURY=/path_to_/merqury
ln -s $MERQURY/merqury.sh

meryl count k=21 filtered<raw_reads>.fastq.gz output reads.meryl
./merqury.sh reads.meryl asm.hap1.fasta asm.hap2.fasta asm.merqury

#plot the kmers 
Rscript /data01/mlandi/software/merqury/plot/plot_spectra_cn.R -f asm.merqury.spectra-asm.hist -o <out> -m 250 -t fill






