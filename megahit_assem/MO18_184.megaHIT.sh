#!/bin/bash
#$ -N MO18_184.MEGA
#$ -ckpt blcr
#$ -j y
#$ -q bio,pub64
#$ -pe openmp 16

module load megahit/1.1.1
module load BBMap/37.50

REF=MO18_184
REFBASE=/dfs3/bio/abchase/moorea
OUTDIR=$REFBASE/megahit

# megahit allows for cross-assembly of the replicates A,B,C

READ1=$REFBASE/${REF}_A.filter.clean.R1.fq.gz,$REFBASE/${REF}_B.filter.clean.R1.fq.gz,$REFBASE/${REF}_C.filter.clean.R1.fq.gz
READ2=$REFBASE/${REF}_A.filter.clean.R2.fq.gz,$REFBASE/${REF}_B.filter.clean.R2.fq.gz,$REFBASE/${REF}_C.filter.clean.R2.fq.gz

cd $OUTDIR

# Megahit - much less resource-intensive than metaSPades or IDBA with comparable results

megahit \
-1 $READ1 \
-2 $READ2 \
-t 16 \
--min-count 3 \
--k-list 31,41,51,61,71,81,91,95,101,105,111 \
--kmin-1pass \
--min-contig-len 1000 \
--memory 0.95 \
--out-dir ${REF} \
--continue


## extra precaution to gets reads into suitable format for binning steps
bbduk.sh in=${REF}/final.contigs.fa out=$OUTDIR/${REF}.contigs.L1kbp.temp.fna minlen=1000 ow=t

rename.sh in=$OUTDIR/${REF}.contigs.L1kbp.temp.fna \
out=$OUTDIR/${REF}.contigs.L1kbp.fna prefix=${REF} addprefix=t ow=t

rm -rf $OUTDIR/${REF}/
rm -f $OUTDIR/${REF}.contigs.L1kbp.temp.fna

	
