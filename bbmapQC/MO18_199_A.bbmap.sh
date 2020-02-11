#!/bin/bash
#$ -N MO18_199_A_QC
#$ -m a
#$ -ckpt restart
#$ -j y
#$ -q bio,pub8i,pub64
#$ -pe openmp 2

module load BBMap/37.50
module load picard-tools/1.96

REF=MO18_199_A
REFBASE=/dfs3/bio/abchase/moorea
FFILE=$REFBASE/rawdata/11922_MO18_199_A_S17_L001_R1_001.fastq.gz
RFILE=$REFBASE/rawdata/11922_MO18_199_A_S17_L001_R2_001.fastq.gz

cd $REFBASE

rm -f $REF.filter.total.fa
rm -f $REF.filter.total.faa
rm -f $REF.gff

bbduk.sh qtrim=rl trimq=10 threads=2 \
in1=$FFILE in2=$RFILE \
out1=$REF.clean1.fq out2=$REF.clean2.fq

repair.sh in=$REF.clean1.fq in2=$REF.clean2.fq \
out=$REF.filter.clean.R1.fq.gz out2=$REF.filter.clean.R2.fq.gz

rm -f $REF.clean1.fq
rm -f $REF.clean2.fq

bbmerge.sh \
in1=$REF.filter.clean.R1.fq.gz in2=$REF.filter.clean.R2.fq.gz \
out=$REF.filter.clean.merged.fq.gz outu=$REF.filter.clean.unmerged.fq.gz

reformat.sh in=$REF.filter.clean.merged.fq.gz out=$REF.filter.clean.merged.fa
reformat.sh in=$REF.filter.clean.unmerged.fq.gz \
out=$REF.filter.clean.unmerged1.fa out2=$REF.filter.clean.unmerged2.fa

## only take the forward read and see if that will annotate, do not want duplicates
cat $REF.filter.clean.merged.fa $REF.filter.clean.unmerged1.fa > $REF.filter.total.fa
rm -f $REF.filter.clean.merged.fq.gz
rm -f $REF.filter.clean.unmerged.fq.gz
rm -f $REF.filter.clean.merged.fa
rm -f $REF.filter.clean.unmerged1.fa
rm -f $REF.filter.clean.unmerged2.fa

prodigal -i $REF.filter.total.fa \
-a $REF.filter.total.faa -q \
-f gff -p meta > $REF.gff

rm -f $REF.gff
rm -f $REF.filter.total.fa

	
