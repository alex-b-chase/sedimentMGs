# sedimentMGs
analysis into the metagenomes derived from Moorea coral-algal gradient

# Metagenome analysis and pipeline

Samples were run by UCSD CMI seed grant - N=20 samples
Samples were run in 3 technical replicates (A, B, C) on an Illumina 

** 1. QC filter **
Filtering was performed with [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)

```
bbduk.sh \
qtrim=rl trimq=10 threads=2 \
in1=$FFILE in2=$RFILE \
out1=$REF.clean1.fq out2=$REF.clean2.fq

repair.sh \
in=$REF.clean1.fq in2=$REF.clean2.fq \
out=$REF.filter.clean.R1.fq.gz out2=$REF.filter.clean.R2.fq.gz

### the *.filter.clean.R*.fq.gz files are the ones to be used for assembly input
```

** 2. Gene-centric approach **
We can use the non-assembled, merged reads to perform a gene-centric approach for taxonomy and function
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

```

