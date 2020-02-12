#!/bin/bash
#$ -N MO18_074.IDBA
#$ -m a
#$ -ckpt restart
#$ -j y
#$ -q bio,pub64
#$ -pe openmp 16

module load joannlp/miniconda2
module load BBMap/37.50

REF=MO18_074
REFBASE=/dfs3/bio/abchase/moorea
OUTDIR=$REFBASE/idba

cd $REFBASE

source activate test

## IDBA needs the reads to be in fastA format 
## we also have technical replicates of the sampleIDs in 3 different copies
## combine these 3 technical replicates and assemble from there

fq2fa --merge <(zcat ${REF}_A.filter.clean.R1.fq.gz) <(zcat ${REF}_A.filter.clean.R2.fq.gz) ${REF}.A.temp.fas
fq2fa --merge <(zcat ${REF}_B.filter.clean.R1.fq.gz) <(zcat ${REF}_B.filter.clean.R2.fq.gz) ${REF}.B.temp.fas
fq2fa --merge <(zcat ${REF}_C.filter.clean.R1.fq.gz) <(zcat ${REF}_C.filter.clean.R2.fq.gz) ${REF}.C.temp.fas

cat ${REF}.A.temp.fas ${REF}.B.temp.fas ${REF}.C.temp.fas > $OUTDIR/${REF}.fas
rm -f ${REF}.*.temp.fas

cd $OUTDIR
rm -rf $OUTDIR/${REF}/

idba_ud -r ${REF}.fas --pre_correction \
--mink 30 --maxk 200 --step 10 --num_threads 16 \
--min_contig 1000 --out ${REF}

source deactivate test

## extra precaution to gets reads into suitable format for binning steps
bbduk.sh in=${REF}/contig.fa out=$OUTDIR/${REF}.contigs.L1kbp.temp.fna minlen=1000 ow=t

rename.sh in=$OUTDIR/${REF}.contigs.L1kbp.temp.fna \
out=$OUTDIR/${REF}.contigs.L1kbp.fna prefix=${REF} addprefix=t ow=t

rm -rf $OUTDIR/${REF}/

	
