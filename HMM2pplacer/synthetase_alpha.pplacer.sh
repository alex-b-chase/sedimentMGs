#!/bin/bash
#$ -N synthetase_alpha_pplacer
#$ -q bio,pub8i
#$ -j y
#$ -pe openmp 8

module load enthought_python/7.3.2
module load ea-utils/r811
module load phylosift/1.0.1
module load clcbio/8.5.1
module load BBMap/37.50
module load hmmer/3.1b2

REFDIR=/dfs3/bio/abchase/moorea/comm_markers
HMMPROF=/dfs3/bio/abchase/refDB/hmmprofiles
REFDB=/dfs3/bio/abchase/refDB/refpkg

OUTDIR=$REFDIR/pplacer 

cd $REFDIR

protein=synthetase_alpha

for minID in {15,20,25,30}
do

	cd $REFDIR

	hmmsearch --tblout ${protein}.hmm.txt -E 1e-${minID} \
	--cpu 8 $HMMPROF/${protein}p.hmm ${protein}.blat.faa > ${protein}.log.txt

	cat ${protein}.hmm.txt | cut -f1 -d' ' | sort | uniq > ${protein}.temp.txt

	# subset new HMMer filtered reads
	filterbyname.sh \
	in=${protein}.blat.faa \
	out=total_${protein}.${minID}.hmm.faa \
	names=${protein}.temp.txt ow=t include=t 2>/dev/null

	rm -f ${protein}.temp.txt
	rm -f ${protein}.hmm.txt
	rm -f ${protein}.log.txt

	cat total_${protein}.${minID}.hmm.faa | tr -d '*' > $OUTDIR/${protein}.${minID}.temp.hmm.faa

	cd $OUTDIR

	rm -f $OUTDIR/${protein}.${minID}.fa

	clustalo-1.2.0 --profile1 $REFDB/${protein}p.refpkg/${protein}p.good.final.aln \
	-i ${protein}.${minID}.temp.hmm.faa -o $OUTDIR/${protein}.${minID}.fa

	rm ${protein}.${minID}.temp.hmm.faa

	### now can test input for pplacer 
	pplacer --pretend \
	-c $REFDB/${protein}p.refpkg \
	${protein}.${minID}.fa > ${protein}.${minID}.pplacer.log

	### check if pplacer worked
	if grep -Fxq "everything looks OK." ${protein}.${minID}.pplacer.log
	then
		echo -e "${protein}\t${minID}" >> $OUTDIR/hmmerfiltered.txt
		break

	else 
		rm -f $REFDIR/total_${protein}.${minID}.hmm.faa
		rm -f $OUTDIR/${protein}.${minID}.fa
		continue

	fi


done

cd $OUTDIR

rm -f ${protein}.*.pplacer.log

### now can run pplacer 
pplacer -c $REFDB/${protein}p.refpkg \
${protein}.${minID}.fa \
-p --keep-at-most 20 

guppy to_csv --point-mass --pp ${protein}.${minID}.jplace > ${protein}.${minID}.csv
guppy fat --node-numbers --point-mass --pp ${protein}.*.jplace


	
