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

Reads will be assigned to the best node in the phylogeny, giving a more accurate, conservative representation of the taxonomic assignment. In understudied systems (i.e., soils and sediments), I would rather be more conservative with these assignments as known, reference genomes in most databases do NOT encompass the diversity of these microbial communities.

**3. Genome-centric approach**

For the assembly step, we will take some insights into this [recent paper](https://www.nature.com/articles/s41564-019-0449-y.pdf?origin=ppub) describing MAGs in another complex microbial community, soil. Soils and sediments both have issues with high species complexity, which can create problems in assembly and binning steps.

So, we will use the [IDBA assembler](https://www.ncbi.nlm.nih.gov/pubmed/22495754) which is optimized for SAGs and MAGs. This software requires specific file format as input so we need to do a little work prior to assembly.

```
## First, we have 3 technical replicates that need to be combined. convert each paired-end, post-QC sample to fastA format
fq2fa --merge <(zcat ${REF}_A.filter.clean.R1.fq.gz) <(zcat ${REF}_A.filter.clean.R2.fq.gz) ${REF}.A.temp.fas
fq2fa --merge <(zcat ${REF}_B.filter.clean.R1.fq.gz) <(zcat ${REF}_B.filter.clean.R2.fq.gz) ${REF}.B.temp.fas
fq2fa --merge <(zcat ${REF}_C.filter.clean.R1.fq.gz) <(zcat ${REF}_C.filter.clean.R2.fq.gz) ${REF}.C.temp.fas

## can then combine
cat ${REF}.A.temp.fas ${REF}.B.temp.fas ${REF}.C.temp.fas > $OUTDIR/${REF}.fas

## ready for IDBA to do its thing
idba_ud \
-r ${REF}.fas --pre_correction \
--mink 30 --maxk 200 --step 10 --num_threads 16 \
--min_contig 1000 --out ${REF}

## the assembly will be output in a folder called ${REF}/contigs.fna

## I like to rename the assembled contigs with the sampleID for easier processing later
## extra precaution to gets reads into suitable format for binning steps

## just make sure IDBA got rid of contigs <1000bp and rename fastA header to append sampleID
bbduk.sh \
in=${REF}/contig.fa \
out=${REF}.contigs.L1kbp.temp.fna minlen=1000 ow=t

rename.sh \
in=${REF}.contigs.L1kbp.temp.fna \
out=${REF}.contigs.L1kbp.fna prefix=${REF} addprefix=t ow=t

```

All scripts to process this step are located in [ibda_assem](ibda_assem/)
