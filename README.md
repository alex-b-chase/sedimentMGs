# sedimentMGs
analysis into the metagenomes derived from Moorea coral-algal gradient

# Metagenome analysis and pipeline

Samples were run by UCSD CMI seed grant - N=20 samples

Samples were run in 3 technical replicates (A, B, C) on an Illumina 

**1. QC filter**

Filtering was performed with [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/). Steps 1 and 2 can be found in [bbmapQC](bbmapQC/)

```
bbduk.sh \
qtrim=rl trimq=10 threads=2 \
minlen=25 ktrim=r k=25 ref=$BBMAPDIR/nextera.fa.gz hdist=1 \
stats=$REF.stats.txt \
in1=$FFILE in2=$RFILE \
out1=$REF.clean1.fq out2=$REF.clean2.fq 

repair.sh \
in=$REF.clean1.fq in2=$REF.clean2.fq \
out=$REF.filter.clean.R1.fq.gz out2=$REF.filter.clean.R2.fq.gz

### the *.filter.clean.R*.fq.gz files are the ones to be used for assembly input
```

**2. Gene-centric approach**

Again we will use BBMap (it literally has everything!). We can use the non-assembled, merged reads to perform a gene-centric approach for taxonomy and function

Some reads cannot be merged, so we will only take the forward reads of those, most of these will be short and probably will not annotate - but worth a shot!

```
bbmerge.sh \
in1=$REF.filter.clean.R1.fq.gz in2=$REF.filter.clean.R2.fq.gz \
out=$REF.filter.clean.merged.fq.gz outu=$REF.filter.clean.unmerged.fq.gz

reformat.sh \
in=$REF.filter.clean.merged.fq.gz out=$REF.filter.clean.merged.fa
reformat.sh \
in=$REF.filter.clean.unmerged.fq.gz out=$REF.filter.clean.unmerged1.fa out2=$REF.filter.clean.unmerged2.fa

## only take the forward read and see if that will annotate, do not want duplicates
cat $REF.filter.clean.merged.fa $REF.filter.clean.unmerged1.fa > $REF.filter.total.fa

## the $REF.filter.total.fa can then be used to call ORFs
## for MGs we need to activate the -p meta option - this will only call ORFs, not annotate reads

prodigal \
-i $REF.filter.total.fa \
-a $REF.filter.total.faa -q \
-f gff -p meta > $REF.gff

## the $REF.filter.total.faa is the file we can parse for specific taxonomic marker genes or functional genes

```

We will pass the ORF predicted file to:
1. Taxonomic Characterization: we will use the taxonomic pipeline outlined [here](https://github.com/alex-b-chase/elevation-community) for community profiles. We can compare these to the 16S community reads.
2. Functional Characterization: NEED TO DO

__2A. Taxonomic Classification__

In brief, we will utilize the phylogenetic classification of 21 taxonomic marker genes against the database described in this paper:
>AB Chase, Z Gomez-Lunar, AE Lopez, J Li, SD Allison, AC Martiny, and JBH Martiny. 2018. Emergence of soil bacterial ecotypes along a climate gradient. Environmental Microbiology.

The database consists of 7392 publicly available genomes that are designated as ‘representative’ genomes by the PATRIC database, plus other soil derived genomic isolates. Processing includes:
1. Run all ORFs with a [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html) filter against the curated protein database
2. Secondary, more stringent filter with Hidden Markov Models ([HMMer](http://hmmer.org/)) against the aligned proteins
3. Phylogenetic placement of each marker gene read to the reference phylogeny with [pplacer](https://matsen.fhcrc.org/pplacer/)

```
blat -prot -fastMap -minIdentity=20 -out=blast8 \
$BLASTDB ${f} temp.20.${sampleID}.txt

# subset the giant output file for only relevant information (i.e. query sequence and protein match)
cut -f1-2 temp.20.${sampleID}.txt | \
awk 'BEGIN{FS="\t"; OFS="\t"} {gsub(/^[^_]*_/, "", $2); print}' | \
sort -u > $OUTDIR/${sampleID}.blat.total.txt

# now subet the marker gene reads from the MG library
cut -f1 $OUTDIR/${sampleID}.blat.total.txt | sort -u | sed 's/[ ]*$//' > $OUTDIR/${sampleID}.blat.temp.txt

# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
filterbyname.sh \
in=${sampleID}.filter.total.faa \
out=$OUTDIR/${sampleID}.blat.markers.faa \
names=$OUTDIR/${sampleID}.blat.temp.txt ow=t include=t 2>/dev/null

rm -f $OUTDIR/${sampleID}.blat.temp.txt
```

Next, we will aggregate the marker genes that were subset from each MG
```
for f in *.blat.total.txt
do

	metagenome=${f%.blat.total.txt}

	# rename each read to reflect the MGID on the fasta header and index each sequence
	cat ${metagenome}.blat.markers.faa | sed -e "s/^>/>${metagenome}_/" | sed '/^>/ s/ .*//' > ${metagenome}.temp.faa
	cat $f | sed -e "s/^/${metagenome}_/" | cut -f1-2 > $metagenome.blat.temp.txt

done

cat *.temp.faa > total.markers.faa 
cat *.blat.temp.txt > total.blatmarkers.txt 
rm -f *.temp.faa
rm -f *.blat.temp.txt

# now subset by each marker gene

while read protein
do
	echo "Processing ${protein}..."

	if [ ! -f "${protein}.blat.faa" ]
	then

		grep -w "${protein}p" total.blatmarkers.txt | cut -f1 | sort -u | sed 's/[ ]*$//' > ${protein}.temp.txt

		# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
		filterbyname.sh \
		in=total.markers.faa  \
		out=${protein}.temp.faa \
		names=${protein}.temp.txt ow=t include=t 2>/dev/null
		rm ${protein}.temp.txt

		# rename each read to index each sequence
		cat ${protein}.temp.faa | awk '/^>/ {$0=NR"_"$0} 1' | \
    sed 's/ //g; s/>//g' | awk '/^[1-9]/ {$0=">"$0} 1' > ${protein}.blat.faa
		rm ${protein}.temp.faa

	else 
		continue 
	fi

done < $BASEDIR/marker_genes.txt
```

Can now parse these through the HMM profiles and test whether the filter was stringent enough. Loop through various e-values for the HMM profiles and test pplacer. Scripts can be found [here](HMM2pplacer/).

Basically, use [hmmsearch](https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch) at various e-values, then align filtered reads to the refDB with [ClustalO](https://www.ebi.ac.uk/Tools/msa/clustalo/), then test [pplacer](https://matsen.fhcrc.org/pplacer/). If e-value threshold was good enough, run pplacer for real.

Reads will be assigned to the best node in the phylogeny, giving a more accurate, conservative representation of the taxonomic assignment. In understudied systems (i.e., soils and sediments), I would rather be more conservative with these assignments as known, reference genomes in most databases do NOT encompass the diversity of these microbial communities.

**3. Genome-centric approach**

For the assembly step, we will take some insights into this [recent paper](https://www.nature.com/articles/s41564-019-0449-y.pdf?origin=ppub) describing MAGs in another complex microbial community, soil. Soils and sediments both have issues with high species complexity, which can create problems in assembly and binning steps.

So, we will use the [IDBA assembler](https://www.ncbi.nlm.nih.gov/pubmed/22495754) which is optimized for SAGs and MAGs. This software requires specific file format as input so we need to do a little work prior to assembly.

IDBA used a lot of memory and wasn't feasible with the number of samples. Scripts can still be found in [ibda_assem](ibda_assem/), but ended up going with [megahit](https://github.com/voutcn/megahit) for this initial assembly.

```
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

```

All scripts to process this step are located in [mega_assem](mega_assem/)
