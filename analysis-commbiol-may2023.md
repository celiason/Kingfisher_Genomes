---
title: kingfisher genomic analyses
output: html_document
editor_options: 
  chunk_output_type: console
---

## Setup

```{r}
setwd("~/Dropbox/Projects/King_genome")

library(ape)  # For phylo analysis, etc.
library(Biostrings)
library(cluster)
library(DOSE)  # for enrichment analysis
library(dplyr)
library(ecodist)
library(enrichplot)  # for enrichment analysis
library(evolqg)
library(geiger)
library(geometry)
library(geomorph)
library(ggmsa)
library(ggplot2)
library(ggridges)
library(ggtree)
library(gprofiler2)
library(grid)
library(igraph)
library(magrittr)
library(MCMCglmm)  # for phylo analysis, etc.
library(mvMORPH)  # for morphometric analysis
library(parallel)
library(pbapply)
library(pbmcapply)
library(phylolm)
library(phytools)
library(readxl)
library(RERconverge)  # for convergence analysis
library(rgl)
library(RPANDA)  # for penalized ou models
library(Rvcg)
library(scales)  # for transparent colors
library(seqinr)
library(stringr)
library(tidyr)
library(treeio)
library(treeplyr)
library(UpSetR)
library(ViSEAGO)
library(visreg)  # for regression visualizations
library(visreg)
library(viridis)
library(tidyverse)

devtools::load_all('~/genomeR')

# Load functions
source("R/gost_to_enrich.R")
source("R/genome-sensory-ms-functions.R")
source("~/Dropbox/Projects/Coraciiform morphometrics/R/tres.R")

# devtools::load_all('~/genomeR')  # <-- or just load this package/include w/paper on my github?

# Phylogenies
# [ ] TODO fix names so consistent throughout
tree <- read.tree("data/paml/kingTree_abbrev.phy")
phy_genomes <- read.tree("data/paml/kingTree_abbrev.phy")

# Trait data sets
Species <- readxl::read_xlsx("~/Dropbox/Projects/king_tongues/data/kingfisher_CT_scans_050817.xlsx")
DietData <- readxl::read_xlsx("~/Dropbox/Projects/king_tongues/data/kingdata_v2_edited.xlsx")
dnadat <- read_excel("~/Dropbox/Projects/king_tongues/Kingfisher DNA list 20190605 updated.xlsx")
brain <- read.csv("~/Dropbox/Projects/Coraciiform morphometrics/data/processed/species_rates.csv", row=1)
brain$abbrev <- sppAbbrev(brain$spp)

treA <- read.tree("~/Dropbox/Projects/King_CT/data/trees/trees_MAndersen/kingTree.BEAST.mean.age.newick.tre")
island_data <- readxl::read_excel("~/Dropbox/Projects/Coraciiform morphometrics/data/coraciiformes-biogeography.xlsx", sheet=5)
island_data$only_island <- as.numeric(island_data$Island==1 & island_data$Continent==0)
# table(island_data$only_island)
islands <- setNames(island_data$only_island, island_data$Species)

# plungeFg <- c('corLeu','corVin','corCri','ceyCya','ceyArg','alcSem','alcQua','alcAtt','megMax','megTor','megAlc','chlInd','chlAma','cerRud')

# Setup "foreground" traits (plunge-diving behavior, island-dwelling)
plungeFg <- c('corLeu','corVin','corCri','ceyCya','ceyArg','alcSem','alcQua',
			  'alcAtt','megMax','megTor','megAlc','chlInd','chlAma','cerRud')

islandFg <- c('corMad','corVin','ceyFal','ceyWeb','ceyMel','ceyCya','ceyArg',
			  'ceyMar','ceyLep', 'ceyCol','ceyNig','actHom','actLin','todAlb',
			  'todDio','todFar','todCin','todRuf')
```

## Genome resequencing

```sh
# todChl 10X reference - only ~50% smaller, so not many dups
hts_SuperDeduper -1 reads/todChl_reads1.fastq.gz -2 reads/todChl_reads2.fastq.gz --stats-file todChl.log -i reads/todChl-pe150-reads
```

```r
# Deduplicate new Illumina resequencing data
setwd("~/uce-alcedinidae/reads")
files <- list.files(".", pattern="pe150-reads\\.fastq")
lapply(files, function(f) {
	system(paste0("hts_SuperDeduper -I ", f, " --stats-file ", substr(f, 1, 6), ".log", " -i ", gsub(".fastq.gz", "", f)))
})
```

## Genome assembly

### Align

```r
spp <- list.files("genomes")

# ref <- "genomes/todChl/todChl.fasta"
# DONE - check this is masked reference-
ref <- "genomes/todChl.scaffolds.full_mask.fa"
# samtools tview  -p 1:700-750 alignments/actHom-to-todChl.rg.indelrealigner.bam repeats/Full_mask/todChl.scaffolds.full_mask.fa

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
reads <- list.files("reads", pattern="INTERLEAVED.fastq.gz", full=TRUE)
names(reads) <- substr(basename(reads), 1, 6)
m2 <- unlist(sapply(spp, grep, substr(basename(reads), 1, 6)))
picks <- spp[!spp %in% c("todChl")]
# picks <- spp[!spp %in% c("actHom", "todChl")]
newreads <- reads[m2[picks]]

# Try this (should work)
# alignReads(ref=ref, reads=newreads[1], cores=48, ram=150, suffix=2, test=T, force=T)

# Do the alignments - estimating ~8h per genome, so ~10 days -- UNCOMMENT TO RUN!
# lapply(newreads, alignReads, ref=ref, cores=48, ram=150)

#-------------------------------------------------------------------------------
# Get stats
#-------------------------------------------------------------------------------
spp <- list.files("genomes")
spp <- spp[spp!="todChl"]
files <- paste0("alignments/", spp, "-to-todChl.bam")
mclapply(files, mc.cores=length(files), function(f) {
	system(paste0("samtools flagstat ", f, " > ", f, ".stats"))
})
```

<!-- We performed INDEL realignment with GATK v. 3.8 (56) -->

```r
# alias gatk="java -jar /home/FM/celiason/gatk-3.8/GenomeAnalysisTK.jar"
# alias picard="java -jar /home/FM/celiason/picard.jar"

REF <- "/home/FM/celiason/uce-alcedinidae/repeats/Full_mask/todChl.scaffolds.full_mask.fa"
# BAM="alignments/alcAtt-to-todChl.bam"

bams <- list.files("alignments", pattern="todChl.bam$", full=T)
# bams <- bams[-c(1, 3)]

mclapply(bams, mc.cores=10, function(bam) {
	sample <- gsub(".bam", "", basename(bam))
	if (file.exists(paste0("alignments/", sample,".rg.bam"))) {
		stop("File exists")
	}
	# Add read groups (~20 mins)
	system(paste0("java -jar /home/FM/celiason/picard.jar AddOrReplaceReadGroups I=", bam, " O=alignments/", sample, ".rg.bam", " RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=Sample1"))
	# Index
	system(paste0("samtools index alignments/", sample, ".rg.bam"))
	
	# Get INDEL intervals (~30 mins)
	system(paste0("java -jar /home/FM/celiason/gatk-3.8/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ", REF, " -I alignments/", sample, ".rg.bam -o alignments/", sample, ".rg.bam.intervals"))
	# realign INDELs (~5 hours)
	system(paste0("java -Xmx8G -Djava.io.tmpdir=/tmp -jar /home/FM/celiason/gatk-3.8/GenomeAnalysisTK.jar -T IndelRealigner -R ", REF, " -targetIntervals alignments/", sample, ".rg.bam.intervals -I alignments/", sample, ".rg.bam -o alignments/", sample, ".rg.indelrealigner.bam"))
})
```

<!-- call genotypes and generate individual consensus sequences -->

```r
# Get list of files
files <- list.files("alignments", pattern="to-todChl.rg.indelrealigner.bam$", full=TRUE)

# Repeat-masked collared kingfisher reference we'll be using-
ref <- "repeats/Full_mask/todChl.scaffolds.full_mask.fa"

# Do the SNV calling and get concensus genome sequences
# Took ~24h for 30 VCF files (~1.3 GB each), consensus calling only 5m each
mclapply(files, mc.cores=length(files), getConsensus, ref=ref, onlyvar=TRUE, force=TRUE)

# getConsensus(files[1], ref="genomes/todChl/todChl.fasta", outpath="redo")
```

## Locating CNEEs

```sh
cd "/home/FM/celiason/uce-alcedinidae/cnee"

ref="/home/FM/celiason/CNEEs/galGal3_all.fa"
# query="/home/FM/celiason/uce-alcedinidae/genomes/todChl/todChl.fasta"
query="/home/FM/celiason/uce-alcedinidae/genomes/todChl/todChl.scaffolds.full_mask.fa"

nm="todChl"

alias lastal="/home/FM/celiason/last-1133/src/lastal"

# Build the lastdb index
lastdb -c galGal $ref

# Run last aligner
parallel-fasta -j 24 "lastal galGal" < $query | maf-convert psl > todChl.psl

# Change coordinates of .psl files to parent coordinate system
mkdir liftup

# Get contig/chromosome lengths
faSize $ref -detailed > galGal3.chr_length.txt
faSize $query -detailed > $nm.chr_length.txt

# Create a lift-format file for later
# column 1: offset, 2: old name, 3: old length, 4: new name, 5: new length
mkdir lift
x=`cat todChl.chr_length.txt | wc -l`
yes 0 | head -n$x | paste - todChl.chr_length.txt todChl.chr_length.txt | column -s $'\t' -t > lift/todChl.lft

# Join close alignments (takes -3 hours)
mkdir chain_raw
axtChain -linearGap=medium -faQ -faT -psl $nm.psl $ref $query ./chain_raw/$nm.chain

# Merge and sort chain files
chainMergeSort chain_raw/*.chain | chainSplit chain_split stdin
  
# Make alignment nets from chain files
mkdir net
ls chain_split/*.chain | parallel --progress --eta -j 32 "chainNet {} galGal3.chr_length.txt $nm.chr_length.txt net/{/}.net /dev/null"

# Create liftOver chain file
mkdir over
for i in ./chain_split/*.chain
do
  tag=${i/\.\/chain_split\//} # getting base name- prob better way?
  netChainSubset net/$tag.net $i over/$tag
done

cat over/*.chain > galGal3-to-todChl.over.chain

# Make bed file to report converted coordinates. We can give the coordinates of our query regions (based on G1 assembly) in the input.bed file and liftOver will report the converted coordinates in conversion.bed file.
liftOver galGal3_top500.bed galGal3-to-todChl.over.chain todChl_top500.bed unMapped
liftOver cnee.bed galGal3-to-todChl.over.chain todChl_cnee.bed unMapped2

wc -l todChl_cnee.bed  # 235559/265983 (88.6%) mapped over
wc -l todChl_top500.bed  # 482/500 mapped over

# bedtools getfasta -fi $query -bed todChl_top500.bed

# Get CNEEs closest to annotation genes (GeMoMa)
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1,$2,$3,$4,$5,$6,$7,$8,$9}' annotations/gemoma/round2/filtered_predictions.gff | convert2bed --input=gff - > annotations/gemoma/round2/filtered_predictions.mRNA.bed

wc -l annotations/gemoma/round2/filtered_predictions.mRNA.bed  # 16540

# bedtools sort -i cnee_redo/todChl_top500.bed | bedtools closest -a - -b annotations/gemoma/round2/filtered_predictions.mRNA.bed > closest_gemoma.txt

# Filter out CNEEs overlapping transcripts in Todiramphus chloris
bedtools intersect -v -a cnee_redo/todChl_cnee.bed -b annotations/gemoma/round2/filtered_predictions.mRNA.bed > cnee_redo/todChl_ASHCE_nonoverlap.bed
wc -l cnee_redo/todChl_ASHCE_nonoverlap.bed # 194298 non-overlapping CNEEs
awk '{if(($3-$2)>=50) print}' cnee_redo/todChl_ASHCE_nonoverlap.bed > cnee_redo/todChl_ASHCE_nonoverlap_50bp.bed
wc -l cnee_redo/todChl_ASHCE_nonoverlap_50bp.bed  # 45613

inbed=cnee_redo/todChl_ASHCE_nonoverlap_50bp.bed

# Find CNEEs closest to each gene, ignore overlap (-io)
bedtools sort -i $inbed | bedtools closest -io -d -a annotations/gemoma/round2/filtered_predictions.mRNA.bed -b - | bedtools sort > cnee_redo/todChl_ASHCE_nonoverlap_closest.bed

# Filter out ones that didn't have a closest CNEE (using `grep -v`)
awk -v OFS="\t" '{print $11"\t"$12"\t"$13"\t"$14}' cnee_redo/todChl_ASHCE_nonoverlap_closest.bed | grep -v "\." | bedtools sort -i - > cnee_redo/todChl_ASHCE_nonoverlap_closest_filt.bed

wc -l cnee_redo/todChl_ASHCE_nonoverlap_closest_filt.bed  # 14344 filtered CNEEs

# awk '{if(($3-$2)>=50 && $4<=200) print $1"\t"$2"\t"$3"\t"$4}' cnee_redo/closest_gemoma_all_sort.bed | wc -l
# awk '{if($4<=200) print $1"\t"$2"\t"$3"\t"$4}' cnee_redo/todChl_ASHCE_nonoverlap_closest.bed | wc -l
# 499 within 200 bp of a coding region

# Closest...

```

<!-- This allowed us to use the reference genome annotations (GFF format) obtained from GeMoMa and extract coding sequences (CDS) for each individual. All CDS, introns, genome alignments, and scripts for data analyses are available at Dryad (DOI pending). -->

```r
# Getting coding sequences (CDS) from consensus pseudo-genomes

setwd("~/uce-alcedinidae")

#-------------------------------------------------------------------------------
# GEMOMA redo with repeat-masked genomes
#-------------------------------------------------------------------------------

# file.copy("repeats/Full_mask/todChl.scaffolds.full_mask.fa", "genomes/todChl/todChl.scaffolds.full_mask.cds.fa")

# Setup parameters
spp <- list.files("genomes")
spp <- spp[spp!="todChl"]
genomes <- paste0("genomes/", spp, "/", spp, "-to-todChl.rg.indelrealigner.consensus.masked.fa")
gff <- "annotations/gemoma/round2/filtered_predictions.gff"  # after running GeMoMa GAF to filter

#-------------------------------------------------------------------------------
# Extract CDS from genomes
#-------------------------------------------------------------------------------
# no early stop codons (-V)
# complete genes (-J)

pblapply(genomes, function(g) {
	prefix <- gsub("\\..*?$", "", g)	
	system(paste0("gffread -J -V -x ", prefix, ".rg.indelrealigner.consensus.masked.cds.complete.fa -g ", g, " ", gff))
	system(paste0("gffread -V -x ", prefix, ".rg.indelrealigner.consensus.masked.cds.fa -g ", g, " ", gff))
})

g <- "genomes/todChl/todChl.scaffolds.full_mask.fa"
prefix <- paste0("genomes/todChl/", gsub(".fa", "", basename(g)))

system(paste0("gffread -J -V -x ", prefix, ".cds.complete.fa -g ", g, " ", gff))
system(paste0("gffread -V -x ", prefix, ".cds.fa -g ", g, " ", gff))

cds <- paste0("genomes/", spp, "/", spp, "-to-todChl.rg.indelrealigner.consensus.masked.cds.fa")
cds <- c(cds, "genomes/todChl/todChl.scaffolds.full_mask.cds.fa")

#-------------------------------------------------------------------------------
# Append species name
#-------------------------------------------------------------------------------
pblapply(cds, function(f) {
	prefix <- substr(basename(f), 1, 6)
	system(paste0("sed -i -e 's/^>/>", prefix, "|/' ", f))
})

#-------------------------------------------------------------------------------
# Convert from species to gene-wise
#-------------------------------------------------------------------------------

# Quick split and merge of species-wise multi fasta files 
# Only takes 1-2 minutes
# Requirement!!! Each multi-fasta file must have the same number and order of sequence names
spp <- list.files("genomes")
spp <- spp[spp!="todChl"]
files <- paste0("genomes/", spp, "/", spp, "-to-todChl.rg.indelrealigner.consensus.masked.cds.fa")
files <- c(files, "genomes/todChl/todChl.scaffolds.full_mask.cds.fa")

# Load sequences
seqs <- pbmcapply::pbmclapply(cds, invisible(seqinr::read.fasta), mc.cores=length(cds))

# Get gene names
ids <- lapply(seqs, function(x) {
	stringr::str_extract(names(x), "(?<=\\|).*?$")
})

#-------------------------------------------------------------------------------
# Find overlapping genes
#-------------------------------------------------------------------------------
picks <- Reduce(intersect, ids)
length(picks) # 11938

#-------------------------------------------------------------------------------
# Subset them
#-------------------------------------------------------------------------------
newseqs <- lapply(seq_along(seqs), function(x) {
	seqs[[x]][ids[[x]] %in% picks]
})

# Check
length(newseqs[[1]]) == length(newseqs[[2]])

#-------------------------------------------------------------------------------
# Output per-gene FASTA files
#-------------------------------------------------------------------------------

outpath <- "output/cds_FINAL"
dir.create(outpath)
pbmcapply::pbmclapply(seq_along(picks), mc.cores=31, function(i) {
	subseq <- sapply(newseqs, "[", i)
	seqinr::write.fasta(subseq, names=names(subseq), file.out = gsub(":", "_", paste0(outpath, "/", picks[i], ".fa")))
})

files <- list.files(outpath, pattern=".fa", full=TRUE)
length(files)

#-------------------------------------------------------------------------------
# Remove stop codons (keeping files in place)
#-------------------------------------------------------------------------------
pbmcapply::pbmclapply(files, mc.cores=48, function(f) {
	seqs <- multi2single(fasta=f)
	seqs <- lapply(seqs, toupper)
	nms <- gsub(">|\\|.*?$", "", names(seqs))
	seqinr::write.fasta(cutstops(seqs), file.out=f, names=nms)
})

# Only complete genes
cds <- paste0("genomes/", spp, "/", spp, "-to-todChl.rg.indelrealigner.consensus.masked.cds.complete.fa")

# Append species name
pblapply(cds, function(f) {
	prefix <- substr(basename(f), 1, 6)
	system(paste0("sed -i -e 's/^>/>", prefix, "|/' ", f))
})

cds <- c(cds, "genomes/todChl/todChl.scaffolds.full_mask.cds.fa")  # this is the repeat-masked genome

# Convert from species to gene-wise

# Load sequences
seqs <- pbmcapply::pbmclapply(cds, invisible(seqinr::read.fasta), mc.cores=8)

# Get gene names
ids <- lapply(seqs, function(x) {
	stringr::str_extract(names(x), "(?<=\\|).*?$")
})

# Find overlapping genes
picks <- Reduce(intersect, ids)
length(picks) # 4143

write.csv(picks, file="output/gemoma_complete_overlapping.csv")
```

<!-- Sensory genes from exonerate -->

```r
gff <- "annotations/kifi_sensory.gff"  # after running GeMoMa GAF to filter

pblapply(genomes, function(g) {
	prefix <- gsub("\\..*?$", "", g)
	# system(paste0("gffread -V -x ", prefix, ".cds2.fa -g ", g, " ", gff)) # -V option removes CDS with stop codons
	system(paste0("gffread -V -x ", prefix, ".cds.sensory.fa -g ", g, " ", gff)) # -V option removes CDS with stop codons
	# system(paste0("bedtools getfasta -fi ", g, " -bed cnee/todChl_cnee.bed > ", prefix, ".noncoding.fa"))
})

# Append species name
cds <- paste0("genomes/", spp, "/", spp, "-to-todChl.cds.sensory.fa")
cds <- c(cds, "genomes/todChl/todChl.cds.sensory.fa")
pblapply(cds, function(f) {
	prefix <- substr(basename(f), 1, 6)
	system(paste0("sed -i -e 's/^>/>", prefix, "|/' ", f))
})

#-------------------------------------------------------------------------------
# 2. Convert from species to gene-wise
#-------------------------------------------------------------------------------

# Load sequences
seqs <- pbmcapply::pbmclapply(cds, invisible(seqinr::read.fasta), mc.cores=length(cds))

# Get gene names
ids <- lapply(seqs, function(x) {
	stringr::str_extract(names(x), "(?<=\\|).*?$")
})

# Find overlapping genes
picks <- Reduce(intersect, ids)
length(picks) # 122 without early stop codons


# Subset them
newseqs <- lapply(seq_along(seqs), function(x) {
	seqs[[x]][ids[[x]] %in% picks]
})

# Output per-gene FASTA files
outpath <- "output/cds_sensory"
dir.create(outpath)
pbmcapply::pbmclapply(seq_along(picks), mc.cores=31, function(i) {
	subseq <- sapply(newseqs, "[", i)
	seqinr::write.fasta(subseq, names=names(subseq), file.out=paste0(outpath, "/", picks[i], ".fa"))
})


files <- list.files(outpath, pattern=".fa", full=TRUE)
length(files)

# Remove stop codons (files in place)
pbmcapply::pbmclapply(files, mc.cores=48, function(f) {
	seqs <- multi2single(fasta=f)
	seqs <- lapply(seqs, toupper)
	nms <- gsub(">|\\|.*?$", "", names(seqs))
	seqinr::write.fasta(cutstops(seqs), file.out=f, names=nms)
})

#-------------------------------------------------------------------------------
# CNEEs for closest GeMoMa-annotated genes
#-------------------------------------------------------------------------------

# Setup parameters
setwd("~/uce-alcedinidae")

spp <- list.files("genomes")
spp <- spp[spp!="todChl"]
genomes <- paste0("genomes/", spp, "/", spp, "-to-todChl.rg.indelrealigner.consensus.masked.fa")
genomes <- c(genomes, "genomes/todChl/todChl.scaffolds.full_mask.fa")
spp <- c(spp, "todChl")

cnee <- "cnee_redo/todChl_ASHCE_nonoverlap_50bp.bed"
# nrow(read.table(cnee))
# cnee <- "cnee_redo/todChl_ASHCE_nonoverlap_closest_filt.bed"
# cnee <- "cnee_redo/closest_gemoma_all.bed"

dir.create("cnee_redo/fasta")

# Get CNEEs
pbmcapply::pbmclapply(seq_along(genomes), mc.cores=length(genomes), mc.preschedule=FALSE, function(x) {
	prefix <- spp[x]
	system(paste0("bedtools getfasta -fi ", genomes[x], " -bed ", cnee, " > cnee_redo/fasta/", prefix, "_ASHCE_gemoma.fa"))
})

# Convert from per-species to per-CNEE alignments for baseml
files <- list.files("cnee_redo/fasta", pattern=".fa$", full=T)
dir.create("cnee_redo/fasta_reshape")
genewise2sppwise(files = files, outpath = "cnee_redo/fasta_reshape", suffix = "_ASHCE_gemoma.fa", force = TRUE)
```

<!-- Prep for baseml models -->

```sh
cd ~/uce-alcedinidae/cnee_redo/fasta_reshape

# Replace colons in files names
# prename -n 's/:/_/g' * # test
prename 's/:/_/g' * # real

mkdir trimmed

cd ~/uce-alcedinidae/cnee_redo

mkdir fasta_reshape_trimmed

# Trim alignments
for f in `find fasta_reshape/`
do
	pre=`basename $f .fa`
	# trimal -resoverlap 0.75 -seqoverlap 75 -in $f -out cnee_redo/fasta_reshape/trimmed/"$pre"_trimmed.fa
	seqkit fx2tab $f | awk -F'N' 'NF<=2' | seqkit tab2fx > fasta_reshape_trimmed/"$pre"_trimmed.fa
done

# Delete files that are 0 bytes
find cnee_redo/fasta_reshape_trimmed -size 0 -print -delete

# Counts to make sure trimming of sequences with N's worked
# cat cnee_redo/fasta_reshape/100:10077415_10077505.fa | grep -c ">"  # 31 species
# cat cnee_redo/fasta_reshape_trimmed/100:10077415_10077505_trimmed.fa | grep -c ">"  # 21 species
```

<!-- Fit baseml models -->

```r

tmux
R

# Run baseml

# dir.create("~/uce-alceinidae/cnee_redo/baseml")

oldwd <- "~/uce-alcedinidae/cnee_redo/baseml"

setwd(oldwd)

# Load packages
library(stringr)
library(pbmcapply)
library(pbapply)
library(txtplot)
library(ape)
library(phylolm)
library(pbapply)
devtools::load_all('~/genomeR')

# List of files with noncoding sequences
files <- list.files("/home/FM/celiason/uce-alcedinidae/cnee_redo/fasta_reshape_trimmed", pattern=".fa", full=TRUE)
length(files)

# Create directories
dir.create("res_global")
dir.create("res_local")
dir.create("res_free")
dir.create("res_free_anc")

#-------------------------------------------------------------------------------
# Run global baseml models
#-------------------------------------------------------------------------------
phy <- ape::read.tree("~/uce-alcedinidae/paml/kingtree_nolabels.phy")
mclapply(files, mc.cores=48, function(f) {
	ctl <- readLines("/home/FM/celiason/uce-alcedinidae/cnee_redo/baseml/baseml_global.ctl")
	ctl <- str_replace(ctl, "seqfile = .*", paste0("seqfile = ", f))
	tips <- stringr::str_match(readLines(f), pattern=">([A-Za-z]+)")[, 2]
	tips <- tips[!is.na(tips)]
	outpath <- paste0("~/uce-alcedinidae/cnee_redo/baseml/res_global/", gsub(".fa", "", basename(f)))
	dir.create(outpath)
	setwd(outpath)
	write.tree(ape::keep.tip(phy, tips), file="king_trim.phy")
	ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = ", "king_trim.phy"))
	cat(ctl, file="temp.ctl", sep="\n")
	system("baseml temp.ctl", ignore.stdout=F)
	setwd(oldwd)
})

#-------------------------------------------------------------------------------
# Run local baseml models
#-------------------------------------------------------------------------------
phy <- ape::read.tree("/home/FM/celiason/uce-alcedinidae/baseml/kingtree_clade.phy") 
mclapply(files, mc.cores=48, function(f) {
	ctl <- readLines("/home/FM/celiason/uce-alcedinidae/cnee_redo/baseml/baseml_local.ctl")
	ctl <- str_replace(ctl, "seqfile = .*", paste0("seqfile = ", f))
	tips <- stringr::str_match(readLines(f), pattern=">([A-Za-z]+)")[, 2]
	tips <- tips[!is.na(tips)]
	outpath <- paste0("~/uce-alcedinidae/cnee_redo/baseml/res_local/", gsub(".fa", "", basename(f)))
	dir.create(outpath)
	setwd(outpath)
	phy$tip.label <- gsub("#1", "", phy$tip.label)
	subphy <- ape::keep.tip(phy, tips)
	subphy$tip.label <- ifelse(subphy$tip.label %in% c('cerRud','chlInd','chlAma','megAlc','megTor','megMax','alcAtt','alcQua','alcSem','ceyArg','ceyCya','corCri','corVin','corLeu'), paste0(subphy$tip.label, "#1"), subphy$tip.label)
	write.tree(subphy, file="king_trim.phy")
	ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = ", "king_trim.phy"))
	cat(ctl, file="temp.ctl", sep="\n")
	system("baseml temp.ctl", ignore.stdout=F)
	setwd(oldwd)
})

#-------------------------------------------------------------------------------
# Run local baseml models
#-------------------------------------------------------------------------------
phy <- ape::read.tree("~/uce-alcedinidae/paml/kingtree_nolabels.phy")
files[1:3]
mclapply(files, mc.cores=48, function(f) {
	ctl <- readLines("/home/FM/celiason/uce-alcedinidae/cnee_redo/baseml/baseml_free_anc.ctl")
	ctl <- str_replace(ctl, "seqfile = .*", paste0("seqfile = ", f))
	tips <- stringr::str_match(readLines(f), pattern=">([A-Za-z]+)")[, 2]
	tips <- tips[!is.na(tips)]
	outpath <- paste0("~/uce-alcedinidae/cnee_redo/baseml/res_free/", gsub(".fa", "", basename(f)))
	dir.create(outpath)
	setwd(outpath)
	write.tree(ape::keep.tip(phy, tips), file="king_trim.phy")
	ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = ", "king_trim.phy"))
	cat(ctl, file="temp.ctl", sep="\n")
	system("baseml temp.ctl", ignore.stdout=F)
	setwd(oldwd)
})

phy <- ape::read.tree("~/uce-alcedinidae/paml/kingtree_nolabels.phy")
mclapply(files, mc.cores=48, function(f) {
	# f=files[1]
	ctl <- readLines("/home/FM/celiason/uce-alcedinidae/cnee_redo/baseml/baseml_free_anc.ctl")
	ctl <- str_replace(ctl, "seqfile = .*", paste0("seqfile = ", f))
	tips <- stringr::str_match(readLines(f), pattern=">([A-Za-z]+)")[, 2]
	tips <- tips[!is.na(tips)]
	outpath <- paste0("~/uce-alcedinidae/cnee_redo/baseml/res_free_anc/", gsub(".fa", "", basename(f)))
	dir.create(outpath)
	setwd(outpath)
	write.tree(ape::keep.tip(phy, tips), file="king_trim.phy")
	ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = ", "king_trim.phy"))
	cat(ctl, file="temp.ctl", sep="\n")
	system("baseml temp.ctl", ignore.stdout=F)
	setwd(oldwd)
})


#-------------------------------------------------------------------------------
# Do likelihood ratio tests
#-------------------------------------------------------------------------------

# ML for global rates models
res1 <- list.files("res_global", pattern="mlb0", recursive=TRUE, full=TRUE)
lik1 <- pbapply::pblapply(res1, function(x) {
	try(raw <- readLines(x))
	res <- stringr::str_match(raw[grep("lnL", raw)], "lnL.*?np.*?(\\d+)\\)\\:\\s+([\\.\\-\\d]+)")[,2:3]
	try(setNames(res, c('np', 'lnL')))
})
lik1[grep("Error", lik1)] <- c(NA, NA)

# ML for local rates models
res2 <- list.files("res_local", pattern="mlb1", recursive=TRUE, full=TRUE)
lik2 <- pblapply(res2, function(x) {
	try(raw <- readLines(x))
	res <- str_match(raw[grep("lnL", raw)], "lnL.*?np.*?(\\d+)\\)\\:\\s+([\\.\\-\\d]+)")[,2:3]
	try(setNames(res, c('np', 'lnL')))
})
lik2[grep("Error", lik2)] <- c(NA, NA)

# Estimate rates on branches for local model
rates <- pbapply::pbsapply(seq_along(res2), function(x) {
	raw <- readLines(res2[[x]])
	ss <- raw[grep("rates for branches", raw)]
	as.numeric(stringr::str_match(ss, "rates for branches\\:.*?([\\d\\.]+).*?([\\d\\.]+)")[, 2:3])
})
rates[grep("numeric", rates)] <- c(NA,NA)
rates <- do.call(rbind, rates)

# ML for free branch rates models
res3 <- list.files("res_free", pattern="mlb2", recursive=TRUE, full=TRUE)
lik3 <- pblapply(res3, function(x) {
	try(raw <- readLines(x))
	res <- str_match(raw[grep("lnL", raw)], "lnL.*?np.*?(\\d+)\\)\\:\\s+([\\.\\-\\d]+)")[,2:3]
	try(setNames(res, c('np', 'lnL')))
})
lik3[grep("Error", lik3)] <- c(NA, NA)

#-------------------------------------------------------------------------------
# Combine into final data set
#-------------------------------------------------------------------------------
x1 <- cbind(file=basename(dirname(res1)), do.call(rbind, lik1))
x1 <- as.data.frame(x1)
x2 <- cbind(file=basename(dirname(res2)), do.call(rbind, lik2))
x2 <- as.data.frame(x2)
x3 <- cbind(file=basename(dirname(res3)), do.call(rbind, lik3))
x3 <- as.data.frame(x3)
x1$lnL <- as.numeric(as.character(x1$lnL))
x2$lnL <- as.numeric(as.character(x2$lnL))
x3$lnL <- as.numeric(as.character(x3$lnL))

#-------------------------------------------------------------------------------
# Calculate log likelihoods (global vs local models)
#-------------------------------------------------------------------------------
degf <- as.numeric(as.character(x2$np[1])) - as.numeric(as.character(x1$np[1]))

pvals <- 1-pchisq(2*(x2$lnL-x1$lnL), df=degf)
pvals <- p.adjust(pvals, method="fdr")
table(pvals<.05)  # 20 significant (used to be 11306!!)

# Prep dataset for plotting, etc.
df <- data.frame(file=x1$file, df=degf, l=-2*(x1$lnL-x2$lnL), p=pvals)
df$rate1 <- rates[,1]
df$rate2 <- rates[,2]
txtplot::txtboxplot(na.omit(df$rate2[df$p<.05]))
t.test(na.omit(df$rate2[df$p<.05]), mu=1) # t=15.23, P<0.001, mean = 1.27X faster in plunge-divers (hmm..)

range(df$rate2,na.rm=T)
# slow evolving-
head(df)
picks <- as.character(df[df$rate2 < 1 & df$p < .05, "file"])

# cat(sort(unique(na.omit(as.character(l$geneID[match(picks, l$cneeID)])))), file="out.txt", sep="\n")
# l <- read.csv("../gemoma_lookup.csv")

# "For each taxon branch lengths were summed from root to tip for each CNE genomic region."

# Output list of regions varying along branches??

# Read the trees
res3 <- list.files("res_free", pattern="mlb2", recursive=TRUE, full=TRUE)
trees <- pblapply(res3, function(x) {
	raw <- readLines(x)
	s <- raw[grep("Detailed output", raw)-2]  # where to find the tree in the output
    cat(s, file = "ex.tre", sep = "\n")
    read.tree("ex.tre")
})
names(trees) <- basename(dirname(res3))
saveRDS(trees, "noncoding_freeratio_trees.rds")

#-------------------------------------------------------------------------------
# Statistical summary
#-------------------------------------------------------------------------------
tips.fg <- c('cerRud','chlInd','chlAma','megAlc','megTor','megMax','alcAtt','alcQua','alcSem', 'ceyArg','ceyCya','corCri','corVin','corLeu')
tips.bg <- setdiff(trees[[1]]$tip, tips.fg)
id.fg <- which(trees[[1]]$tip %in% tips.fg)
id.bg <- which(trees[[1]]$tip %in% tips.bg)
res <- lapply(seq_along(trees), function(tr) {
	# tr=10
	# Get branch lengths for foreground/background species:
	fg <- trees[[tr]]$edge.length[trees[[tr]]$edge[, 2] %in% id.fg]
	bg <- trees[[tr]]$edge.length[trees[[tr]]$edge[, 2] %in% id.bg]
	# Calculate mean difference in rates
	mean(fg) - mean(bg)
})

#-------------------------------------------------------------------------------
# Statistical testing
#-------------------------------------------------------------------------------

# [ ] TODO - need to run this through the RERconverge pipeline (but which trees, free ratio ones?)

t.test(unlist(res))

txtboxplot(as.numeric(res))  # no diffs

trees <- readRDS("~/Downloads/noncoding_freeratio_trees.rds")
tips.fg <- c('cerRud','chlInd','chlAma','megAlc','megTor','megMax','alcAtt','alcQua','alcSem', 'ceyArg','ceyCya','corCri','corVin','corLeu')
id.fg <- which(trees[[1]]$tip %in% tips.fg)

sprates <- sapply(1:length(trees), function(x) {
	try(tres(trees[[x]])[trees[[x]]$tip])
})
names(sprates) <- names(trees)
sprates <- sprates[-which(sapply(sprates, class)=="try-error")]

phy <- read.tree("/Users/chadeliason/Dropbox (The Field Museum)/Projects/King_genome/data/paml/kingTree_abbrev.phy")

forage <- factor(ifelse(trees[[1]]$tip %in% tips.fg, "fish", "nofish"), levels=c('nofish', 'fish'))

fits <- pblapply(1:length(sprates), function(x) {
	dat <- data.frame(rate=sprates[[x]], forage=forage)
	phylolm(rate ~ forage, data = dat, phy = phy, method="lambda")
})

pvals <- sapply(fits, function(x) summary(x)$coef[1, "p.value"])
newpvals <- p.adjust(pvals, method="fdr")
table(newpvals < 0.05)

betas <- sapply(fits, function(x) summary(x)$coef[1, "Estimate"])
names(betas) <- names(sprates)

MASS::truehist(betas, xlim=c(-1, 1), h=.001)

keep <- betas > quantile(betas, .025) & betas < quantile(betas, .975)

newbetas <- betas[keep]
newpvals <- newpvals[keep]

# Plot some interesting ones
sort(newbetas,decreasing=TRUE)[1:25]
plot(trees[["73_6707049-6813084"]], edge.col=ifelse(trees[[1]]$edge[,2] %in% id.fg, "red", "black"), edge.width=2)

l <- read.csv("~/Downloads/gemoma_lookup.csv")

cat(unique(na.omit(l$geneID[match(names(sprates)[which(newpvals < 0.05)], l$cneeID)])), file="~/Desktop/cnee_gemoma_fish.txt", sep="\n")

sumtree <- trees[[1]]
x <- do.call(cbind, sapply(trees, "[[", "edge.length"))
sumtree$edge.length <- apply(x, 1, mean)
plot(sumtree, edge.col=ifelse(trees[[1]]$edge[,2] %in% id.fg, "red", "black"), edge.width=2)
```

<!-- # Output CNEE trees for RERconverge: -->

```r
files <- list.files("~/uce-alcedinidae/cnee_redo/baseml/res_free/", pattern="mlb2", recursive=TRUE, full=TRUE)
trees <- sapply(files, function(f) {
	raw <- readLines(f)
	raw[grep("tree length", raw)[1]+4]
})
trees <- trees[!is.na(trees)]
trees <- trees[!grepl("nan", trees)]
nms <- basename(dirname(names(trees)))
res <- cbind(nms, trees)
numsp <- sapply(res[, 2], str_count, "[a-z][a-zA-Z]+")
res2 <- res[numsp==31, ]
write.table(res2, file="~/uce-alcedinidae/cnee_redo/ASHCE_trees.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=F)
```

## Reference bias assessment

```r
devtools::load_all('~/genomeR')

#-------------------------------------------------------------------------------
# Use faidx to N-mask reference
#-------------------------------------------------------------------------------

# List of species we'll test mapping bias for, representing each of major clades
picks <- c('alcAtt', 'actLin', 'megAlc')

REF <- "genomes/todChl/todChl.fasta"

mclapply(picks, mc.cores=length(picks), function(i) {
	VCF <- paste0("genomes/", i, "/", i, "-to-todChl.calls.vcf.gz")
	filtfile <- paste0("genomes/", i, "/", i, ".variants.vcf.gz")
	maskfile <- paste0("genomes/", i, "/todChl_", i, "_SNP_Nmask.fasta")
	if (file.exists(filtfile) | file.exists(maskfile)) {
		warning("files exists")
	}
	# Nmask variable sites in reference
	system(paste0('bcftools view -H -i "QUAL>30 && DP>5" ', VCF, ' | gzip > ', filtfile))
	system(paste0("cp genomes/todChl/todChl.fasta ", maskfile))  # fast
	# Takes ~15 mins for 1GB VCF-
	system(paste0("zcat ", filtfile, " | vcf2bed --do-not-sort | faidx ", maskfile, " --bed - -m -s N"))
	# check it worked-
	# zcat variants.vcf.gz | head -n15400
	# samtools faidx todChl_Nmask.fasta 46:1690-1699
})

#-------------------------------------------------------------------------------
# Remap to N-masked reference
#-------------------------------------------------------------------------------

lapply(picks, function(i) {
	rf=paste0("genomes/", i, "/todChl_", i, "_SNP_Nmask.fasta")
	rd=paste0("reads/", i, "-pe150-reads_INTERLEAVED.fastq.gz")
	alignReads(reads=rd, ref=rf, cores=48, ram=150, suffix=".Nmasked")
})

#-------------------------------------------------------------------------------
# Use featureCounts to calculate read coverage at GENES and NONCODING regions
#	 before after N-masking)
#-------------------------------------------------------------------------------

lapply(picks, function(i) {
	GFF="annotations/todChl_rnd3.all.maker.noseq.renamed.putative_function.domain_added.gff"
	BAM1=paste0("alignments/", i, "-to-todChl.bam")
	BAM2=paste0("alignments/", i, "-to-todChl.Nmasked.bam")
	system(paste0("featureCounts -T 12 -t exon -g ID -F GFF -a ", GFF, " -o ", i, "1.txt ", BAM1))
	system(paste0("featureCounts -T 12 -t exon -g ID -F GFF -a ", GFF, " -o ", i, "2.txt ", BAM2))
})

#-------------------------------------------------------------------------------
# 5. Plot results in R
#-------------------------------------------------------------------------------

x1 <- read.table("~/uce-alcedinidae/counts1.txt", head=T)
x2 <- read.table("~/uce-alcedinidae/counts2.txt", head=T)

# plot(density(x1$actHom_run1.bam))
# diffs <- (x1$actHom_run1.bam-x2$actHom_run2.bam) / (x1$actHom_run1.bam+x2$actHom_run2.bam)
df <- cbind(x1$actHom_run1.bam, x2$actHom_run2.bam)
df <- df[complete.cases(df),]
df <- df[df[,1]!=0 & df[,2]!=0, ]
df <- log(df)
df <- as.data.frame(df)

ggplot(df, aes(x=V1, y=V2)) + geom_point() + scale_x_log10() + scale_y_log10() + annotation_logticks(side="bl") + geom_abline(slope=1, intercept=0, lty=2) + labs(x="Unmasked per-exon coverage", y="Masked per-exon coverage") + theme_minimal()

t.test(df[,1], df[,2], paired=T)

boxplot(diffs)
abline(h=0, lty=2)

hist(diffs, breaks=250)

# Here we're loading from "Downloads" folder but on phoebe there are in ~/uce-alcedinidae folder
x1 <- read.table("~/Downloads/alcAtt1.txt", head=T)
x2 <- read.table("~/Downloads/alcAtt2.txt", head=T)
# x1 <- read.table("~/uce-alcedinidae/alcAtt1.txt", head=T)
# x2 <- read.table("~/uce-alcedinidae/alcAtt2.txt", head=T)
x3 <- read.table("~/Downloads/actLin1.txt", head=T)
x4 <- read.table("~/Downloads/actLin2.txt", head=T)
x5 <- read.table("~/Downloads/megAlc1.txt", head=T)
x6 <- read.table("~/Downloads/megAlc2.txt", head=T)

# plot(density(x1$actHom_run1.bam))
# diffs <- (x1$actHom_run1.bam-x2$actHom_run2.bam) / (x1$actHom_run1.bam+x2$actHom_run2.bam)
df1 <- cbind(x1[,7], x2[,7])
df1 <- df1[complete.cases(df1),]
df1 <- df1[df1[,1]!=0 & df1[,2]!=0, ]
df1 <- log(df1)
df1 <- as.data.frame(df1)

df2 <- cbind(x3[,7], x4[,7])
df2 <- df2[complete.cases(df2),]
df2 <- df2[df2[,1]!=0 & df2[,2]!=0, ]
df2 <- log(df2)
df2 <- as.data.frame(df2)

df3 <- cbind(x5[,7], x6[,7])
df3 <- df3[complete.cases(df3),]
df3 <- df3[df3[,1]!=0 & df3[,2]!=0, ]
df3 <- log(df3)
df3 <- as.data.frame(df3)

ggplot(exp(df1), aes(x=V1, y=V2)) + geom_point() + scale_x_log10() + scale_y_log10() + annotation_logticks(side="bl") + geom_abline(slope=1, intercept=0, lty=2) + labs(x="Unmasked per-exon coverage", y="Masked per-exon coverage") + theme_minimal()
ggsave("covplot_alcAtt.jpg")

# Stats
t.test(df[,1], df[,2], paired=T) # P<0.01, mean diff=-0.019
wilcox.test(df[,1], df[,2], paired=T)
# lower coverage for non-masked?

#-------------------------------------------------------------------------------
# Generate plots
#-------------------------------------------------------------------------------

# setwd("~/Dropbox/Projects/King_genome")

jpg("figs/nmask_coverage.jpg", width=6.5, height=2)
par(mfrow=c(1, 3), mar=c(1,1,1,1), mgp=c(1.5,.5,0), mex=1, ps=10, oma=c(3,3,0,0))
MASS::truehist(exp(df1[,1]), xlim=xl, h=5, col=alpha('red', .5))
par(new=T)
MASS::truehist(exp(df1[,2]), xlim=xl, h=5, col=alpha('blue', .5), axes=F, xlab="")
legend("topright", bty="n", legend=c("unmasked", "masked"), col=c('red','blue'), pch=15)
MASS::truehist(exp(df2[,1]), xlim=xl, h=5, col=alpha('red', .5))
par(new=T)
MASS::truehist(exp(df2[,2]), xlim=xl, h=5, col=alpha('blue', .5), axes=F, xlab="")
legend("topright", bty="n", legend=c("unmasked", "masked"), col=c('red','blue'), pch=15)
MASS::truehist(exp(df3[,1]), xlim=xl, h=5, col=alpha('red', .5))
par(new=T)
MASS::truehist(exp(df3[,2]), xlim=xl, h=5, col=alpha('blue', .5), axes=F, xlab="")
legend("topright", bty="n", legend=c("unmasked", "masked"), col=c('red','blue'), pch=15)
mtext(side=1, outer=T, line=1.5, text="Reads per exon")
mtext(side=2, outer=T, line=1.5, text="Density")
dev.off()
```

## PSMC analyses

### Runs

```r
# NB: DO NOT RUN UNLESS YOU HAVE A LOT OF TIME!!!

# Setup PSMCFA parameters
REF <- "genomes/todChl/todChl.scaffolds.full_mask.fa"
spp <- list.files("genomes")
spp <- spp[spp!="todChl"]

#-------------------------------------------------------------------------------
# Generate PSMCFA files-
#-------------------------------------------------------------------------------
mclapply(spp, mc.cores=length(spp), function(x) {
	bam <- paste0("alignments/", x, "-to-todChl.rg.indelrealigner.bam")
	meancovg <- as.numeric(as.character(tail(read.table(paste0("mosdepth/", x, "-to-todChl.mosdepth.summary.txt")), 1)$V4))
	mincov <- round((1/3) * meancovg)
	maxcov <- round(meancovg * 2)
	system(paste0("samtools mpileup -q 20 -Q 20 -uf ", REF, " ", bam, " |\
		bcftools call -V indels -c |\
		vcfutils.pl vcf2fq -d ", mincov, " -D ", maxcov, " -Q 20 |\
		fq2psmcfa -g0 - > PSMC/", gsub(".bam", ".psmcfa", basename(bam))))
})

#-------------------------------------------------------------------------------
# Run PSMC (using bird genomes paper settings: N30 –t5 –r5 –p 4+30*2+4+6+10)
#-------------------------------------------------------------------------------
mclapply(spp, mc.cores=length(spp), function(x) {
	system(paste0("psmc -N30 -t5 -r5 -p'4+30*2+4+6+10' -o PSMC/", x, ".psmc PSMC/", x, "-to-todChl.rg.indelrealigner.psmcfa"))
})

#-------------------------------------------------------------------------------
# Run PSMC (using penguin genome paper settings: N25 -t15 -r5 -p 4+25*2+4+6)
#-------------------------------------------------------------------------------
mclapply(spp, mc.cores=length(spp), function(x) {
	system(paste0("psmc -N25 -t15 -r5 -p'4+25*2+4+6' -o PSMC/", x, "_penguin.psmc PSMC/", x, "-to-todChl.rg.indelrealigner.psmcfa"))
})

# Settings from REF XX-
# -N25 -t15 -r5 -p 4 + 25*2 + 4 + 6

# Reference genome PSMC
# psmc -d -N30 -t65 -r5 -p'3+2*17+15*1+1*12' -o PSMC/test.psmc PSMC/alcAtt-to-todChl.rg.indelrealigner.psmcfa

# Run PSMC bootstraps
mclapply(1:100, mc.cores=50, function(x) {
	run <- paste0("./psmc/psmc -N30 -t5 -r5 -b -p'4+30*2+4+6+10' -o todChl_round-", x, ".psmc todChl_split.psmcfa")
	system(run)
})

# AnAge databse - animal longetivity
# 2-5 yrs generation time probably ok
# x1 <- psmc.result(file="data/psmc/todChl_diploid.psmc", mu=4.6e-9, g=5, i.iteration=30)

# Mutation rate = 0.10719441 (Lanfear 2010 for Alcedo) / [ 24e6 (Prum 2015 tree) / 5 (gen. time from xx) ] = 2.23e-8
```

### Analyze results

```{r}
# Mutation rates inferred from dS trees across 16000+ genes-
trees <- readRDS("output/dS_trees.rds")
res <- pbsapply(trees, tres)  # this computes total branch lengths for each species
mus <- apply(res, 1, mean, na.rm=T) / (37.5e6/5) # 37.5 Ma from McCullough et al. 2019, all gen times = 5

# Loading results from PSMC analyses
files <- list.files("output/all_psmc",pattern=".psmc",full=T)
files <- files[grepl("penguin", files)]
names(files) <- gsub("_penguin.psmc|.psmc", "", basename(files))
files <- c(files, todChl="data/psmc/todChl_combined.psmc")
psmc_res <- lapply(seq_along(files), function(x) {
	psmc.result(file=files[x], mu=mus[names(files)[x]], g=5, i.iteration=25)
})
names(psmc_res) <- names(files)
psmc_res <- lapply(psmc_res, function(x) as.data.frame(x[[1]]))

# Create data set for graphing
dat_psmc <- plyr::ldply(psmc_res, .id="spp")
dat_psmc$gen <- substr(dat_psmc$spp, 1, 3)
dat_psmc$geo <- ifelse(dat_psmc$spp %in% c("corMad","corVin","ceyFal","ceyWeb","ceyMel","ceyCya","ceyArg","ceyMar","ceyLep", "ceyCol","ceyNig","actHom","actLin","todAlb","todDio","todFar","todCin","todRuf"), "island", "cont")
dat_psmc <- dat_psmc %>% group_by(spp) %>% dplyr::slice(9:n())
dat_psmc <- dat_psmc %>% mutate(clade = case_when(gen %in% c('alc','cey','cor','isp') ~ "pygmy", gen %in% c('act','hal','tod') ~ "forest", gen %in% c('cer','chl','meg') ~ "river"))

# Preliminary plots
ggplot(filter(dat_psmc, gen!="mom"), aes(x=YearsAgo, y=Ne, color=gen, group=spp, linetype=geo)) +
	geom_step() +
	scale_y_log10() +
	scale_x_log10() +
	facet_wrap(~clade, nr=1) +
	annotation_logticks(side="bl") +
	theme_minimal() +
	coord_cartesian(x=c(1e4, 5e6), y=c(1e3, 1e6), expand=F)
ggsave("figs/psmc_by_clade.jpg", width=11, height=4)

ggplot(data=dat_psmc, aes(x=YearsAgo, y=Ne)) +
	geom_step() +
	scale_y_log10() +
	scale_x_log10() +
	facet_wrap(~spp) +
	annotation_logticks(side="bl") +
	coord_cartesian(x=c(1e4, 5e6), y=c(1e3, 1e6), expand=F)
ggsave("figs/psmc_all_species2.jpg",width=11,height=7.5)

# Create data set summarizing max Ne for statistical analysis
dat_psmc_sum <- dat_psmc %>% group_by(spp,geo) %>% dplyr::summarize(Ne.max = max(Ne))
dat_psmc_sum$spp <- as.character(dat_psmc_sum$spp)
dat_psmc_sum$brainsize <- setNames(brain$ln.brainsize, brain$abbrev)[dat_psmc_sum$spp]
dat_psmc_sum$lnmass <- setNames(brain$lnmass, brain$abbrev)[dat_psmc_sum$spp]
dat_psmc_sum$brainPC1 <- setNames(brain$brainPC1, brain$abbrev)[dat_psmc_sum$spp]
dat_psmc_sum$brainPC2 <- setNames(brain$brainPC2, brain$abbrev)[dat_psmc_sum$spp]
dat_psmc_sum$brainPC3 <- setNames(brain$brainPC3, brain$abbrev)[dat_psmc_sum$spp]
dat_psmc_sum$specrate <- setNames(brain$ln.speciationCLADS, brain$abbrev)[dat_psmc_sum$spp]
dat_psmc_sum$diver <- setNames(brain$dives, brain$abbrev)[as.character(dat_psmc_sum$spp)]

# Plot some relationships with phenotypes
# pairs(dat_psmc_sum[c('brainsize','lnmass','Ne.max')])
# ggplot(dat_psmc_sum, aes(x=brainPC1, y=log(Ne.max))) + geom_point() + stat_smooth(method="lm")
# ggplot(dat_psmc_sum, aes(x=brainPC2, y=log(Ne.max))) + geom_point() + stat_smooth(method="lm")
# ggplot(dat_psmc_sum, aes(x=brainPC3, y=log(Ne.max))) + geom_point() + stat_smooth(method="lm")
# ggplot(dat_psmc_sum, aes(x=brainsize, y=log(Ne.max))) + geom_point() + stat_smooth(method="lm")
# ggplot(dat_psmc_sum, aes(x=specrate, y=log(Ne.max))) + geom_point() + stat_smooth(method="lm")
# ggplot(dat_psmc_sum, aes(x=factor(diver), y=log(Ne.max))) + geom_boxplot()
# summary(lm(log(Ne.max) ~ brainPC1, data = dat_psmc_sum))

# Stats
td <- dat_psmc_sum %>% na.omit %>% make.treedata(tree=phy)
dat <- as.data.frame(td$dat)
rownames(dat) <- td$phy$tip
lm1 <- phylolm(log(Ne.max) ~ lnmass, data = dat, phy = td$phy, model="lambda")
summary(lm1)
# phytools::phylomorphospace(dat[c('lnmass','Ne.max')], tree=td$phy)

#-------------------------------------------------------------------------------
# Clustering analysis of time series
#-------------------------------------------------------------------------------

spp <- unique(dat_psmc$spp)
fits <- pblapply(unique(dat_psmc$spp), function(x) {
	ss <- subset(dat_psmc, spp==x)
	approx(ss$YearsAgo, ss$Ne, xout=seq(1e4, 1e6, length=500))
})
ts <- sapply(fits, "[[", "y")
ts <- log(sapply(fits, "[[", "y"))
ts <- apply(ts, 2, minmax,na.rm=T)
ts <- cbind(time=fits[[1]]$x, ts)
colnames(ts)[-1] <- as.character(unique(dat_psmc$spp))

dmat <- daisy(t(ts[,-1]), metric = "euclidean")  # allows for missing values
hc <- hclust(dmat)
# using frey bc they came up with kmeans clustering
nb <- NbClust::NbClust(diss=dmat, distance=NULL, method="ward.D2", index="silhouette")
grp <- nb$Best.partition

# Plot dendrogram + PSMC curves
df2 <- as.data.frame(ts)
df2[,-1] <- apply(df2[,-1], 2, minmax, na.rm=TRUE)
df2 <- df2 %>% gather(id, Ne, -time)
# Do we need this so that order will match up between tree and curves??
# df2$id <- factor(df2$id, levels=names_order)
df2 <- df2 %>% select(id, time, Ne)
df2$island <- ifelse(df2$id %in% c("corMad","corVin","ceyFal","ceyWeb","ceyMel","ceyCya","ceyArg","ceyMar","ceyLep", "ceyCol","ceyNig","actHom","actLin","todAlb","todDio","todFar","todCin","todRuf"), "island", "cont")
df2 <- df2[!is.na(df2$Ne), ]
df2$gen <- substr(df2$id, 1, 3)
df2$diver <- setNames(brain$dives, brain$abbrev)[df2$id]
df2 <- df2 %>% mutate(clade = case_when(gen %in% c('alc','cey','cor','isp') ~ "pygmy", gen %in% c('act','hal','tod') ~ "forest", gen %in% c('cer','chl','meg') ~ "river"))

p1 <- ggtree(phy) + geom_tiplab(cex=2)
p2 <- facet_plot(p1, 'Trait', data=df2, geom=geom_ridgeline, scale=1, mapping=aes(x=time, height=Ne, group=label, alpha=.25, color=island, fill=NULL)) + xlim_expand(xlim=c(1e4, 1e6), panel="Trait") + theme_bw()
p2 + xlim_expand(xlim=c(0, 40), panel="Tree")

# p2 <- facet_plot(p1, data=df2, geom=geom_point, panel="hi", mapping=aes(x=time, group=id, y=Ne)) + theme_void()

# This isn't working now- delete it?
p2 + scale_y_continuous(limits = c(0, 31), expand = c(0,0), oob = function(x, ...) x) +
	geom_text(aes(label=lab), data=d) + 
	coord_cartesian(clip='off')  + 
	theme(plot.margin=margin(6, 6, 40, 6))
ggsave("figs/psmctree.pdf", width=6.5, height=8)


# Another way to plot
facet_plot(p1, 'Trait', data=df2, geom=geom_ridgeline, scale=1.05,
           mapping=aes(x=time, height=Ne, group=label,
                       alpha=0,
                       color=factor(diver),
                       fill=factor(diver))) +
	xlim_expand(xlim=c(1e4, 1e6), panel="Trait") +
  xlim_expand(xlim=c(0, 35), panel="Tree") +
	scale_x_log10() +
	annotation_logticks(side="b") +
	theme(legend.position="none", text=element_text(size=8), panel.grid=element_blank(), strip.background=element_blank(), strip.text=element_blank())

# a base plot way-

# Create plunge-diving behavior vector
fish <- brain$dives[match(phy_genomes$tip, brain$abbrev)]
names(fish) <- phy_genomes$tip
fish[c('momMom','todWin')] <- 0
fish[c('alcSem')] <- 1

island <- ifelse(phy_genomes$tip %in% c("corMad","corVin","ceyFal","ceyWeb","ceyMel","ceyCya","ceyArg","ceyMar","ceyLep", "ceyCol","ceyNig","actHom","actLin","todAlb","todDio","todFar","todCin","todRuf"), "island", "cont")
names(island) <- phy_genomes$tip.label


# Do stochastic character mapping
tree <- drop.tip(treA, which(!treA$tip.label %in% names(islands)))
islands <- islands[tree$tip]

nms <- sppAbbrev(treA$tip.label)

sm <- make.simmap(tree=phy_genomes, x=fish[phy_genomes$tip], model="ARD", nsim=100)

isl2 <- setNames(islands, sppAbbrev(names(islands)))
paste0(isl2[names(fish)], fish)
sort(names(isl2))

obj = densityMap(sm)
# plotSimmap(obj)
obj$cols[1:1001] <- colorRampPalette(c("gray","blue"), space="Lab")(1001)
# plot(obj)

tipnames <- gsub("_", " ", treA$tip.label[match(sm[[1]]$tip.label, nms)])

# Color theme we'll use
pal <- RColorBrewer::brewer.pal(8, "Set2")

#-------------------------------------------------------------------------------
# Create figure for manuscript
#-------------------------------------------------------------------------------

nodes <- sapply(list(c('corLeu','corCri'), c('ceyCya','ceyArg'), c('alcSem','alcAtt'), c('cerRud','megMax')), getMRCA, phy=phy_genomes)

pal2 <- c(brewer.pal(8, "Pastel2")[c(8,4)], brewer.pal(8, "Dark2")[c(8,4)])

plot(1:7, col=brewer.pal(7, "YlGnBu"), pch=16, cex=3)
pal2=c('lightgray','gray25', brewer.pal(7, "YlGnBu")[c(2,6)])
plot(1:4, col=pal2, pch=16, cex=3)

# Paint branches
phy_mapped <- phy_genomes
phy_mapped <- paintSubTree(phy_mapped, node=getMRCA(phy_mapped, c('corLeu','corCri')), state=3, stem=TRUE)
phy_mapped <- paintBranches(phy_mapped, edge=which(phy_mapped$tip=="corMad"), state=2)
phy_mapped <- paintBranches(phy_mapped, edge=which(phy_mapped$tip=="corVin"), state=4)
phy_mapped <- paintSubTree(phy_mapped, node=getMRCA(phy_mapped, c('ceyArg','ceyMel')), state=2, stem=TRUE)
phy_mapped <- paintSubTree(phy_mapped, node=getMRCA(phy_mapped, c('ceyCya','ceyArg')), state=4, stem=TRUE)
phy_mapped <- paintSubTree(phy_mapped, node=getMRCA(phy_mapped, c('cerRud','megMax')), state=3, stem=TRUE)
phy_mapped <- paintSubTree(phy_mapped, node=getMRCA(phy_mapped, c('alcSem','alcAtt')), state=3, stem=TRUE)
phy_mapped <- paintBranches(phy_mapped, edge=which(phy_mapped$tip=="todCin"), state=2)
phy_mapped <- paintSubTree(phy_mapped, node=getMRCA(phy_mapped, c('actLin','actHom')), state=2, stem=TRUE)

# fish island
# testpal <- c('#00afc0', '#78c3cf', '#b5d8dd', '#ececec', '#b5b5b5', '#717171', '#333333')
# testpal <- testpal[c(4,7,3,1)]

# testpal <- c('#00a9da', '#7cbfe0', '#b7d5e6', '#ececec', '#b5b5b5', '#717171', '#333333')
# testpal <- testpal[c(4,7,3,1)]

pal2 <- c('lightgray', 'gray25', brewer.pal(9, "YlGnBu")[c(4,7)])
plot(1:length(pal2), col=pal2, cex=3, pch=15)

pal2 <- c('lightgray', 'gold', 'blue', 'green3')
# pal2 <- c('lightgray', 'gray25', 'blue', 'green3')
plot(1:4, col=pal2, pch=16, cex=2)

pdf("figs/focal_tree.pdf", width=1.5, height=4)
plot(phy_mapped, col=setNames(pal2, 1:4), fsize=0.5, mar=c(1,0,0,0), direction="leftwards", lwd=1.5)
# plot(phy_mapped, col=setNames(pal2, 1:4), fsize=0.5, mar=c(1,0,0,0), direction="leftwards", lwd=3)
# axis(side=1, cex=.5, lwd=.5, lwd.tick=.5, cex.axis=.5)
par(mgp=c(1.5,.5,0), mex=.75, ps=8)
axisPhylo(lwd=0.5, lwd.ticks=.25, gap.axis=1, cex.axis=.5, tcl=-.15, padj=-3.5)
# mtext(side=1, "Time (Ma)", line=2, cex=)
dev.off()

tipcols <- pal2[factor(paste0(fish, island))]

tipnames <- gsub("_", " ", treA$tip.label[match(phy_mapped$tip.label, sppAbbrev(treA$tip))])

pdf("figs/fig1_psmc_islandsNEW.pdf", width=4, height=6)
layout(cbind(1,2,3), widths=c(2, 1, 1))
# plot(obj, mar=c(2,0,0,0), ftype="off", offset=0.5, lwd=1, legend=F)
plot(phy_mapped, col=setNames(pal2, 1:4), offset=0.5, lwd=2, mar=c(2,0,0,0), ftype="off")
nodelabels(text=1:4, node=nodes, frame='circ', bg='white', cex=.65)
# legend("bottomleft", bty="n", legend=c('cont. non-plunge-diving', 'insular non-plunge-diving', 'cont. plunge-diving', 'insular plunge-diving'), col=pal2, pch=16, cex=0.65)
legend("bottomleft", bty="n", legend=c('background', 'island-dwelling', 'mainland plunge-diver', 'island plunge-diver'), col=pal2, pch=16, cex=0.65)
tiplabels(pch=16, cex = 3e7*mus, col = pal2[factor(paste0(fish, island))])
par(mar=c(2,0,0,0))
space <- 0.75 # scaling factor to give space between lines
for (i in seq_along(phy_genomes$tip)) {
  pick <- phy_genomes$tip[i]
  dat <- subset(dat_psmc, spp==pick & YearsAgo >= 1e4 & YearsAgo <= 2e6)
  dat$Ne <- minmax(log(dat$Ne))
  yvar <- dat$Ne * space + i
  xvar <- log(dat$YearsAgo, base=10)
  if (i == 1) {
    # plot(yvar ~ xvar, type='S', log='x', ylim=c(1, Ntip(phy_genomes)), yaxt='n', xaxt='n', bty='n', col=pal[grp[pick]], xaxs='i', xlim=c(4, 6))
    plot(yvar ~ xvar, type='S', log='x', ylim=c(1, Ntip(phy_genomes)), yaxt='n', xaxt='n', bty='n', col=tipcols[i], xaxs='i', xlim=c(4, 6))

    # plot(yvar ~ xvar, type='S', log='x', ylim=c(1, Ntip(phy)), yaxt='n', xaxt='n', bty='n', lty=c(1,3)[grp[pick]], xaxs='i', xlim=c(4, 6))
    abline(h=i, lwd=0.5, col='gray')
  } else {
    # lines(yvar ~ xvar, type='S', col=pal[grp[pick]])
    lines(yvar ~ xvar, type='S', col=tipcols[i])
    # lines(yvar ~ xvar, type='S', lty=c(1,3)[grp[pick]])
    abline(h=i, lwd=0.5, col='gray')
  }
}
axis(side=1, line=-0.5, at=c(4,5,6), labels=c('10 Ky', '100 Ky', '1 My'), cex.axis=.5, tck=-0.01, padj=-4)
par(mar=c(2,0,0,0))
plot(0, type='n', xlim=c(0, 1), ylim=c(1, 31), bty='n', xaxt='n', yaxt='n')
text(x=0, y=1:Ntip(phy_mapped), labels=tipnames, adj=-1, font=3, pos=4, cex=0.6)
dev.off()

ggplot(df, aes(x=YearsAgo, y=Ne, color=spp)) + geom_line() + scale_y_log10(lim=c(1e3, 1e7)) + scale_x_log10(limits=c(17200, 368000))

#-------------------------------------------------------------------------------
# Statistical analysis of clustering
#-------------------------------------------------------------------------------

ranges <- df %>% group_by(spp) %>% summarize(tmin=min(YearsAgo), tmax=max(YearsAgo)) %>% summarize(max(tmin), min(tmax)) %>% as.numeric

l <- split(df, df$spp)

N <- 50
xout <- seq(ranges[1], ranges[2], length=N)
A <- lapply(seq_along(l), function(i) {approx(x=l[[i]]$YearsAgo, y=log(l[[i]]$Ne), xout=xout)$y})
A <- array(unlist(A), dim = c(N, 1,  length(A)))
dimnames(A)[[3]] <- names(l)
A <- A[,,phy$tip]

fish <- fish[phy$tip]

# Fish-eating behavior
fit.psmc.fish <- procD.pgls(t(A[,phy$tip]) ~ fish[phy$tip], phy = phy, iter=999)
summary(fit.psmc.fish)  # P = 0.25
mean.fish <- apply(A[, names(fish)[fish==1]],1,mean)
se.fish <- apply(A[, names(fish)[fish==1]],1,se)
mean.nofish <- apply(A[, names(fish)[fish==0]],1,mean)
se.nofish <- apply(A[, names(fish)[fish==0]],1,se)

# Island-dwelling
xvar <- setNames(ifelse(phy$tip %in% c("corMad","corVin","ceyFal","ceyWeb","ceyMel","ceyCya","ceyArg","ceyMar","ceyLep", "ceyCol","ceyNig","actHom","actLin","todAlb","todDio","todFar","todCin","todRuf"), "island", "cont"), phy$tip.label)
fit.psmc.island <- procD.pgls(t(A) ~ xvar, phy = phy, iter=999)
summary(fit.psmc.island)  # P = 0.009

mean.isl <- apply(A[, names(xvar)[xvar=="island"]],1,mean)
se.isl <- apply(A[, names(xvar)[xvar=="island"]],1,se)
mean.cont <- apply(A[, names(xvar)[xvar=="cont"]],1,mean)
se.cont <- apply(A[, names(xvar)[xvar=="cont"]],1,se)

pdf("figs/ne_islands.pdf", width=5, height=4)
par(mfrow=c(1, 1), mar=c(3,3,1,1), mgp=c(1.5,.5,0), mex=0.75, ps=8)
plot(mean.isl ~ xout, type='l', col=2, xlab="Mya", ylab="log Ne", ylim=c(8.5, 11))
lines(mean.isl + se.isl ~ xout, col=2, lty=2)
lines(mean.isl - se.isl ~ xout, col=2, lty=2)
lines(mean.cont ~ xout)
lines(mean.cont + se.cont ~ xout, lty=2)
lines(mean.cont - se.cont ~ xout, lty=2)
legend("topright", lwd=1, col=1:2, legend=c('Continent', 'Island'), bty='n')
dev.off()
```

## Assessing genome-wide convergence.

We tested the hypothesis that plunge-diving behavior results in genome-wide convergence following Sackton et al. (10). Briefly, we determined the total number of similar (i.e. changes in specific amino acids having evolved either convergently or through parallel evolution) and divergent amino acids for each pair of lineages. In calculating convergent and divergent sites, we opted to treat lineages as the whole path from an ancestor to tip rather than single branches (59).

```r
setwd("~/uce-alcedinidae")

source("genome-functions.R") # on my local machine :)

files <- list.files("paml/M0_results", full=T)
files <- files[!grepl("old", files)]
rstfiles <- paste0(files, "/rst")
mlcfiles <- paste0(files, "/M0_mlc")

# res = codeml result read in with read.codeml
res <- read.codeml(rstfiles[2], mlcfiles[2])
res@phylo$edge

# fish edges - 22, 38, 44, 51
# fish ALL edges - combn(c(22, 38, 44, 51, 23:32, 39:40, 45:48, 52:55), m=2)
# 23:32
# 39:40
# 45:48
# 52:55

tmp <- pbmclapply(seq_along(files), mc.cores=24, mc.preschedule=FALSE, function(i) {
	tryCatch(getChanges_edges(res = read.codeml(rstfiles[i], mlcfiles[i])), error = function(e) NULL)
})
names(tmp) <- basename(files)
saveRDS(tmp, file="convergence_analysis_edges.rds")

res <- plyr::ldply(tmp, .id="ID")
res$edge <- rep(1:nrow(tmp[[2]]))

df <- data.frame(conv = tapply(res$conv, res$edge, sum),
				 div = tapply(res$div, res$edge, sum),
				 par = tapply(res$par, res$edge, sum))

# df$fish <- ifelse(idx[1, ] %in% c(22,38,44,51) & idx[2, ] %in% c(22,38,44,51), "fishy", "not")
df$fish <- ifelse(idx[1, ] %in% c(22, 38, 44, 51, 23:32, 39:40, 45:48, 52:55) & idx[2,] %in% c(22, 38, 44, 51, 23:32, 39:40, 45:48, 52:55), "fishy", "not")

pdf("convtest.pdf")
ggplot(df, aes(x=div, y=conv+par, color=fish)) +
	geom_point() +
	geom_abline(slope=1, intercept=0, lty=2) +
	theme_bw()
dev.off()

# Takes a little while... (a few hours I think)
bigres <- pbmclapply(seq_along(files), mc.cores = 24, mc.preschedule = FALSE, function(i) {
	tryCatch(getChanges_tips(res = read.codeml(rstfiles[i], mlcfiles[i])), error = function(e) NULL)
})
names(bigres) <- basename(files)
saveRDS(bigres, file="paml/convergence_analysis.rds")
# bigres[["galGal_rna-XM_025142692.1_R0"]]  # Is this Tas1r3??


# for baseml stuff
library(treeio)
files <- list.files("~/uce-alcedinidae/cnee_redo/baseml/res_free_anc", full=T)
rstfiles <- paste0(files, "/rst")
# takes ~15 minutes
res_baseml <- pbmclapply(seq_along(files), mc.cores=48, mc.preschedule=FALSE, function(x) {
	tryCatch(getChanges_tips_baseml(res = read.paml_rst(rstfiles[x])), error = function(e) NULL)
})
names(res_baseml) <- basename(files)
saveRDS(res_baseml, file="~/uce-alcedinidae/cnee_redo/baseml/convergence_baseml.rds")
```

We then used phylogenetic linear mixed models implemented in MCMCglmm (60) to test whether the number of convergent changes was explained by foraging behavior, island dwelling, divergence time, number of divergent amino acid substitutions, or the interaction between divergent substitutions and divergence time. Given potential inflation of effect sizes due to non-independent samples (e.g., single lineage compared to multiple other lineages), we treated each lineage in a pairwise comparison as a random effect in the model, after shuffling lineage pairs to ensure even occurrence of each lineage in a comparison (61). We also accounted for phylogenetic non-independence by treating the most recent common ancestor (MRCA) of each pair as a random effect, with an expected covariance structure determined by the phylogenetic tree calculated with the ginverse function in MCMCglmm. We ran the model for 106 generations and discarded the first 25% as burnin. Phylogenetic signal was calculated as the posterior mean of the node random effect divided by the sum of all random effects scaled by dividing by π2/3 (61).

```r
#------------------------------------------------------------------------------
# Coding regions
#------------------------------------------------------------------------------

bigres <- readRDS("output/convergence_analysis.rds")

tr <- ape::read.tree("data/paml/kingTree_abbrev.phy")

fish <- setNames(rep(0, Ntip(tr)), tr$tip)
fish[c('megTor','megAlc','megMax','chlAma','chlInd','cerRud','ceyCya','ceyArg','alcSem','alcQua','alcAtt','corVin','corCri','corLeu')] <- 1
table(fish)

island <- ifelse(tr$tip %in% c("corMad","corVin","ceyFal","ceyWeb","ceyMel","ceyCya","ceyArg","ceyMar","ceyLep", "ceyCol","ceyNig","actHom","actLin","todAlb","todDio","todFar","todCin","todRuf"), 1, 0)
names(island) <- tr$tip

res <- plyr::ldply(bigres, .id="ID")
res$pair <- paste0(res$spp1, "-", res$spp2)
res$fish1 <- fish[res$spp1]
res$fish2 <- fish[res$spp2]
res$island1 <- island[res$spp1]
res$island2 <- island[res$spp2]

# Both lineages plunge-divers
res$newfish <- ifelse(res$fish1==1 & res$fish2==1, 1, 0)

# Both lineages island-dwelling
res$newisland <- ifelse(res$island1==1 & res$island2==1, 1, 0)

# Don't want to compare these because they all have the same diving phenotype
d1 <- combn(c('megMax','megTor','megAlc','chlAma','chlInd','cerRud'), m = 2)
d1 <- apply(d1,2,sort)
d1 <- apply(d1, 2, paste0, collapse="-")
d2 <- combn(c('corLeu','corVin','corCri'), m=2)
d2 <- apply(d2,2,sort)
d2 <- apply(d2, 2, paste0, collapse="-")
d3 <- c('ceyArg-ceyCya')
d4 <- combn(c('alcSem','alcQua','alcAtt'), m=2)
d4 <- apply(d4,2,sort)
d4 <- apply(d4, 2, paste0, collapse="-")
res$newfish[res$pair %in% c(d1,d2,d3,d4)] <- 2
table(res$newfish)

# Don't want to compare these because they all have the same island phenotype
d1 <- combn(c('ceyMel','ceyMar','ceyCya','ceyArg'), m = 2)
d1 <- apply(d1, 2,sort)
d1 <- apply(d1, 2, paste0, collapse="-")
d2 <- c('actHom-actLin')
res$newisland[res$pair %in% c(d1,d2)] <- 2
table(res$newisland)

dat_changes <- data.frame(conv = tapply(res$conv, res$pair, sum),
				 		  div = tapply(res$div, res$pair, sum),
				 		  par = tapply(res$par, res$pair, sum))

dat_changes$fish <- res$newfish[match(rownames(dat_changes), res$pair)]
dat_changes$island <- res$newisland[match(rownames(dat_changes), res$pair)]

dat_changes$converg_fish <- ifelse(dat_changes$fish==2, 0, 1)
dat_changes$converg_island <- ifelse(dat_changes$island==2, 0, 1)

dat_changes$fish <- factor(ifelse(dat_changes$fish==2, 1, dat_changes$fish))
dat_changes$island <- factor(ifelse(dat_changes$island==2, 1, dat_changes$island))

phydist <- ape::cophenetic.phylo(tr)

dat_changes[c('sp1','sp2')] <- do.call(rbind, strsplit(rownames(dat_changes), "-"))

# Getting phylogenetic distances from cophenetic result sorting by species name
dat_changes$phydist <- sapply(1:nrow(dat_changes), function(x) {
	phydist[dat_changes$sp1[x], dat_changes$sp2[x]]
})

ngenes <- length(bigres)

dat_changes$convprop <- dat_changes$conv/ngenes
dat_changes$divprop <- dat_changes$div/ngenes
dat_changes$relconv <- dat_changes$conv/dat_changes$div

#------------------------------------------------------------------------------
# Regulatory regions
#------------------------------------------------------------------------------

bigres_baseml <- readRDS("output/convergence_baseml.rds")

res_baseml <- plyr::ldply(bigres_baseml, .id="ID")
res_baseml$pair <- paste0(res_baseml$spp1, "-", res_baseml$spp2)
res_baseml$fish1 <- fish[res_baseml$spp1]
res_baseml$fish2 <- fish[res_baseml$spp2]
res_baseml$island1 <- island[res_baseml$spp1]
res_baseml$island2 <- island[res_baseml$spp2]

# Both lineages plunge-divers
res_baseml$newfish <- ifelse(res_baseml$fish1==1 & res_baseml$fish2==1, 1, 0)

# Both lineages island-dwelling
res_baseml$newisland <- ifelse(res_baseml$island1==1 & res_baseml$island2==1, 1, 0)

# Don't want to compare these because they all have the same diving phenotype
d1 <- combn(c('megMax','megTor','megAlc','chlAma','chlInd','cerRud'), m = 2)
d1 <- apply(d1,2,sort)
d1 <- apply(d1, 2, paste0, collapse="-")
d2 <- combn(c('corLeu','corVin','corCri'), m=2)
d2 <- apply(d2,2,sort)
d2 <- apply(d2, 2, paste0, collapse="-")
d3 <- c('ceyArg-ceyCya')
d4 <- combn(c('alcSem','alcQua','alcAtt'), m=2)
d4 <- apply(d4,2,sort)
d4 <- apply(d4, 2, paste0, collapse="-")
res_baseml$newfish[res_baseml$pair %in% c(d1,d2,d3,d4)] <- 2

# Don't want to compare these because they all have the same island phenotype
d1 <- combn(c('ceyMel','ceyMar','ceyCya','ceyArg'), m = 2)
d1 <- apply(d1, 2,sort)
d1 <- apply(d1, 2, paste0, collapse="-")
d2 <- c('actHom-actLin')
res_baseml$newisland[res_baseml$pair %in% c(d1,d2)] <- 2

dat_changes_baseml <- data.frame(conv = tapply(res_baseml$conv, res_baseml$pair, sum),
								 div = tapply(res_baseml$div, res_baseml$pair, sum),
								 par = tapply(res_baseml$par, res_baseml$pair, sum))

dat_changes_baseml$fish <- res_baseml$newfish[match(rownames(dat_changes_baseml), res_baseml$pair)]
dat_changes_baseml$island <- res_baseml$newisland[match(rownames(dat_changes_baseml), res_baseml$pair)]
dat_changes_baseml$converg_fish <- ifelse(dat_changes_baseml$fish==2, 0, 1)
dat_changes_baseml$converg_island <- ifelse(dat_changes_baseml$island==2, 0, 1)
dat_changes_baseml$fish <- factor(ifelse(dat_changes_baseml$fish==2, 1, dat_changes_baseml$fish))
dat_changes_baseml$island <- factor(ifelse(dat_changes_baseml$island==2, 1, dat_changes_baseml$island))
dat_changes_baseml[c('sp1','sp2')] <- do.call(rbind, strsplit(rownames(dat_changes_baseml), "-"))
# Getting phylogenetic distances from cophenetic result sorting by species name
dat_changes_baseml$phydist <- sapply(1:nrow(dat_changes_baseml), function(x) {
	phydist[dat_changes_baseml$sp1[x], dat_changes_baseml$sp2[x]]
})
ngenes <- length(bigres_baseml)
dat_changes_baseml$convprop <- dat_changes_baseml$conv/ngenes
dat_changes_baseml$divprop <- dat_changes_baseml$div/ngenes
dat_changes_baseml$relconv <- dat_changes_baseml$conv/dat_changes_baseml$div
rm(bigres_baseml,res_baseml)

save(dat_changes, dat_changes_baseml, file="output/convergence_results.rda")


##

load(file="output/convergence_results.rda")

df <- rbind(data.frame(dat_changes[names(dat_changes_baseml)], type="coding"), data.frame(dat_changes_baseml, type="regulatory"))
head(df)

df$fish[df$converg_fish==0] <- 2
df$island[df$converg_island==0] <- 2
df$group <- factor(paste0(df$fish, df$island))
table(df$group)

df$fish_bg <- plyr::revalue(df$fish, c('2'='0'))
df$island_bg <- plyr::revalue(df$island, c('2'='0'))

df$plotgroup <- factor(paste0(df$fish_bg, df$island_bg))
levels(dat_changes$plotgroup) <- c('background', 'island-dwelling', 'continental\nplunge-diver', 'island\nplunge-diver')

head(df)
ggplot(data=df, aes(divprop, convprop, color=fish_bg)) + geom_point() + stat_smooth(method="lm") + facet_wrap(~type, scales="free") + scale_color_manual(values=c('gray','blue'))

lm1 <- lm(convprop ~ divprop * phydist + fish_bg + island_bg, data = subset(df, type=="coding"))
summary(lm1)
lm2 <- lm(convprop ~ divprop * phydist + fish_bg + island_bg, data = subset(df, type=="regulatory"))
summary(lm2)

#------------------------------------------------------------------------------
# Statistical analyses
#------------------------------------------------------------------------------

# Fish stats
lm1 <- lm(convprop ~ divprop * phydist + fish + island, data = dat_changes)
summary(lm1) # model P < 0.001, fish P < 0.01
# r2full <- summary(lm1)$r.squared
# r2red <- summary(update(lm1, ~.-fish))$r.squared
# (r2full-r2red) / (1-r2full)

library(effectsize)
aov.table <- car::Anova(lm1, type=3)
effectsize::eta_squared(aov.table)  # 0.06 for fish-eating (= medium effect)

par(mfrow=c(1,2))
visreg::visreg2d(lm1, 'divprop', 'phydist')
visreg::visreg(lm1, 'divprop', by='fish', overlay=T)
visreg::visreg(lm1, 'fish')

# Island stats
lm2 <- lm(convprop ~ divprop * phydist + island, data=dat_changes)
summary(lm2) # model P < 0.001, fish P < 0.01

par(mfrow=c(1,2))
visreg::visreg2d(lm2, 'div', 'phydist')
visreg::visreg(lm2, 'island')

#-------------------------------------------------------------------------------
# MCMCglmm approach
#-------------------------------------------------------------------------------

# Adding outgroup for MCMCglmm
tr$root.edge <- 10
tr2 <- bind.tip(tree=tr, tip.label="outgroup", edge.length = 10+max(branching.times(tr)), position=10)

# Name nodes
tr2$node.label <- as.character((Ntip(tr2)+1) : (Ntip(tr2) + Nnode(tr2)))

dat_changes$node <- sapply(1:nrow(dat_changes), function(x) as.character(getMRCA(tr2, c(dat_changes$sp1[x], dat_changes$sp2[x]))))

# Convert expected covariance
ainv <- inverseA(tr2)$Ainv

# Prior for MCMClgmm
prior <- list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),R=list(V=1,nu=0.02))

## next, sort species names between the "sp1" and "sp2" random effects (as in Tobias et al. 2014), so that each species appears more or less the same number of times in "species1" and "species2"
##NB this takes several hours to run; can be saved and reloaded 
# lineages <- as.matrix(data.frame(df$sp1,df$sp2))
# lineages.sorted <- random.effect.sorting(lineages, counter.max=1e6)
# saveRDS(lineages.sorted, file="output/lineage_sorting.rds")
lineages.sorted <- readRDS(file="output/lineage_sorting.rds")
lineages.sorted$asymmetry
##note lineages.sorted$asymmetry shows the balance for each species between
# appearing first and second in random effects (==0 if species occurs the same
# amoutn of time as 'sp1' and 'sp2' in random effects)--this is about as good
# as we can do, and is mostly evenly spread

##lineages.sorted$sorted.data provides the balanced order for each species

dat_changes$sp1.sorted <- lineages.sorted$sorted.data[,1]
dat_changes$sp2.sorted <- lineages.sorted$sorted.data[,2]

# Needs to set levels to filter out species not in focal set
levels(dat_changes$fish) <- c('0','1','2')
levels(dat_changes$island) <- c('0','1','2')
dat_changes$fish[dat_changes$converg_fish==0] <- 2
dat_changes$island[dat_changes$converg_island==0] <- 2
dat_changes$group <- factor(paste0(dat_changes$fish, dat_changes$island))
table(dat_changes$group)
# 341 "background" species pairs (not plunge-diving, not islands)
# 27 pairs islands only
# 6 pairs islands but not convergent
# 67 pairs plunge-diving only
# 2 pairs BOTH
# 21 pairs plunge-diving but NOT convergent
# 1 pair BOTH but non-convergent

sum(dat_changes$group %in% c('00','01','10','11')) # 437 species pairs in target set
dim(dat_changes) # 465 total

dat_changes$fish_bg <- plyr::revalue(dat_changes$fish, c('2'='0'))
dat_changes$island_bg <- plyr::revalue(dat_changes$island, c('2'='0'))

# Fit model with both plunge-diving and island-dwelling predictors (uncomment to run)
set.seed(1980)
# fit_mcmc_bm <- MCMCglmm(convprop ~ (divprop * phydist) + fish_bg + island_bg, data = dat_changes, random = ~sp1.sorted+sp2.sorted+node, ginverse=list(node=ainv), nitt=1e6, burnin=0.25 * 1e6, prior=prior)
# summary(fit_mcmc_bm)
# saveRDS(fit_mcmc_bm, file="output/mcmc_fit.rds")

prior0 <- list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)),R=list(V=1,nu=0.02))

fit_mcmc_wn <- MCMCglmm(convprop ~ divprop + fish_bg + island_bg, data = dat_changes, random = ~sp1.sorted + sp2.sorted + node, nitt=1e5, burnin=0.25*1e5, prior=prior, ginverse=list(node=ainv))
summary(fit_mcmc_wn)
summary(fit_mcmc_bm)

dat_changes$clade <- NA
# dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 34)] <- "cor"
# dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 35)] <- "kin"
# dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 51)] <- "alc"
dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 59)] <- "alc1"
dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 52)] <- "alc2"
# dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 46)] <- "cer"
dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 49)] <- "cer1"
dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 47)] <- "cer2"
# dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 37)] <- "hal"
dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 42)] <- "hal1"
dat_changes$clade[dat_changes$node %in% getDescendants(tr2, 38)] <- "hal2"

ggplot(subset(dat_changes, !is.na(clade)), aes(x=divprop, y=convprop, color=clade)) + geom_point() + stat_smooth(method="lm", se=F)

summary(lm(convprop ~ (divprop * phydist) + fish + island, data = dat_changes))

fit1 <- MCMCglmm(convprop ~ (divprop * phydist) + fish, data = subset(dat_changes, converg_fish==1), random = ~sp1.sorted+sp2.sorted+node, ginverse=list(node=ainv), nitt=1e5, burnin=0.25 * 1e5, prior=prior)

summary(fit1)


fit_mcmc_bm <- readRDS(file="output/mcmc_fit.rds")
summary(fit_mcmc_bm)

out <- summary(fit3)$solutions

stats <- data.frame(estimate=round(out[,1],3),
	"(95% CI)" = paste0("(", round(out[,2],3), ", ", round(out[,3],3), ")"),
	pMCMC = ifelse(out[,"pMCMC"] < 0.001, 0.001, round(out[,"pMCMC"], 3)))

# Output for table
stats

# Phylogenetic signal
mean(fit3$VCV[,'node']/(rowSums(fit3$VCV)+(pi^2/3)))
coda::HPDinterval(fit3$VCV[,'node']/(rowSums(fit3$VCV)+(pi^2/3)))

pairs(df[c('convprop','divprop','phydist','fish')])

# Plots

# ggplot(data=dat_changes, aes(x=divprop, y=convprop)) +
# 	geom_point(aes(color=group)) +
# 	scale_color_manual(values=pal2) +
# 	stat_smooth(aes(linetype=fish), method="lm", se=F, col='black') +
# 	labs(x="Number of divergent sites per gene", y="Number of convergent sites per gene") +
# 	scale_linetype_manual(values = c(2, 1)) +
# 	theme_minimal() +
# 	theme(legend.position="none", text=element_text(size=8))
# ggsave("figs/gw_convergence.pdf", width=3.5, height=2.8)
# ggsave("figs/gw_convergence.jpg", width=5, height=4)

names(pal2) <- c('00','01','10','11')

# names(testpal) <- c('00','01','10','11')

bigpal <- c(pal2, pal2[-1])
names(bigpal)[c(5:7)] <- c('02','20','22')
bigpal_fill <- bigpal
bigpal_fill[5:7] <- "white"

ggplot(data=dat_changes, aes(x=phydist, y=convprop)) +
	# stat_smooth(aes(color=grp), method="lm", se=F) + 
	stat_smooth(data = subset(dat_changes, converg_fish==1 & converg_island==1), aes(linetype=fish), method="lm", se=F, col='black', lwd=0.5) +
	geom_point(data = subset(dat_changes, group %in% c('00','01','02')), aes(color=group, fill=group), pch=21, cex=1.5) +
	geom_point(data = subset(dat_changes, group %in% c('10','11','20','22')), aes(color=group, fill=group), pch=21, cex=1.5) +
	scale_color_manual(values = bigpal) +
	scale_fill_manual(values = bigpal_fill) +
	labs(x="Divergence time (My)", y="Number of convergent sites per gene") +
	scale_linetype_manual(values = c(2, 1)) +
	theme_minimal() +
	# theme_bw() +
	# theme(legend.position="none", text=element_text(size=8), panel.grid=element_blank())
	theme(legend.position="none", text=element_text(size=8))
# ggsave("figs/gw_convergence.jpg", width=5, height=4)
# ggsave("figs/gw_convergence.pdf", width=5, height=4)


pdf("figs/gw_convergence.pdf", width=5, height=4)
par(mgp=c(1.5,.5,0),mar=c(3,3,1,1),mex=.75,ps=8)
plot(convprop ~ phydist, data = arrange(dat_changes, group), col=bigpal[as.character(group)], bg=bigpal_fill[as.character(group)], pch=21, xlab="Divergence time (My)", ylab="Convergent sites per gene", bty='l')
abline(lm(convprop~phydist,subset(dat_changes, converg_fish==1 & converg_island==1 & fish==0)), lty=2)
abline(lm(convprop~phydist,subset(dat_changes, converg_fish==1 & converg_island==1 & fish==1)))
legend("bottomright", bty="n", legend=c('background', 'island-dwelling', 'continental plunge-diver', 'island plunge-diver'), col=pal2, pch=16, cex=1)
dev.off()


# ggplot(data=dat_changes, aes(x=phydist, y=convprop)) +
# 	geom_point(data=subset(dat_changes, fish==0), aes(color=group)) +
# 	geom_point(data=subset(dat_changes, fish==1), aes(color=group)) +
# 	scale_color_manual(values=pal2) +
# 	stat_smooth(aes(linetype=fish), method="lm", se=F, col='black') +
# 	labs(x="Divergence time (My)", y="Number of convergent sites per gene") +
# 	scale_linetype_manual(values = c(2, 1)) +
# 	theme_minimal() +
# 	# theme_bw() +
# 	# theme(legend.position="none", text=element_text(size=8), panel.grid=element_blank())
# 	theme(legend.position="none", text=element_text(size=8))
# ggsave("figs/gw_convergence.pdf", width=3.5, height=2.8)


# ggplot(data=dat_changes, aes(x=divprop, y=convprop, color=factor(island))) +
# 	geom_point() +
# 	stat_smooth(method="lm", se=F) +
# 	labs(x="Number of divergent sites per gene", y="Number of convergent sites per gene") +
# 	scale_color_manual(values=  c('0'='gray', '1'='blue')) +
# 	theme_classic() +
# 	theme(legend.position="none")
# ggsave("figs/gw_convergence_islands.pdf", width=5, height=4)
```

alternative approaches

```r
tmp <- subset(df, type=="regulatory")

# Convergent changes matrix
convdist_reg <- matrix(NA, nrow=31, ncol=31)
rownames(convdist_reg) <- rownames(phydist)
colnames(convdist_reg) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	convdist_reg[id1, id2] <- tmp$conv[i]
	convdist_reg[id2, id1] <- tmp$conv[i]
}
diag(convdist_reg) <- 0

# Divergent changes matrix
divdist_reg <- matrix(NA, nrow=31, ncol=31)
rownames(divdist_reg) <- rownames(phydist)
colnames(divdist_reg) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	divdist_reg[id1, id2] <- tmp$div[i]
	divdist_reg[id2, id1] <- tmp$div[i]
}
diag(divdist_reg) <- 0

# Fish matrix
fishdist_reg <- matrix(NA, nrow=31, ncol=31)
rownames(fishdist_reg) <- rownames(phydist)
colnames(fishdist_reg) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	fishdist_reg[id1, id2] <- as.numeric(tmp$fish_bg[i])-1
	fishdist_reg[id2, id1] <- as.numeric(tmp$fish_bg[i])-1
}
diag(fishdist_reg) <- 0

# Island matrix
islanddist_reg <- matrix(NA, nrow=31, ncol=31)
rownames(islanddist_reg) <- rownames(phydist)
colnames(islanddist_reg) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	islanddist_reg[id1, id2] <- as.numeric(tmp$island_bg[i])-1
	islanddist_reg[id2, id1] <- as.numeric(tmp$island_bg[i])-1
}
diag(islanddist_reg) <- 0


tmp <- subset(df, type=="coding")

# Convergent changes matrix
convdist_cod <- matrix(NA, nrow=31, ncol=31)
rownames(convdist_cod) <- rownames(phydist)
colnames(convdist_cod) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	convdist_cod[id1, id2] <- tmp$conv[i]
	convdist_cod[id2, id1] <- tmp$conv[i]
}
diag(convdist_cod) <- 0

# Divergent changes matrix
divdist_cod <- matrix(NA, nrow=31, ncol=31)
rownames(divdist_cod) <- rownames(phydist)
colnames(divdist_cod) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	divdist_cod[id1, id2] <- tmp$div[i]
	divdist_cod[id2, id1] <- tmp$div[i]
}
diag(divdist_cod) <- 0

# Fish matrix
fishdist_cod <- matrix(NA, nrow=31, ncol=31)
rownames(fishdist_cod) <- rownames(phydist)
colnames(fishdist_cod) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	fishdist_cod[id1, id2] <- as.numeric(tmp$fish_bg[i])-1
	fishdist_cod[id2, id1] <- as.numeric(tmp$fish_bg[i])-1
}
diag(fishdist_cod) <- 0

# Island matrix
islanddist_cod <- matrix(NA, nrow=31, ncol=31)
rownames(islanddist_cod) <- rownames(phydist)
colnames(islanddist_cod) <- colnames(phydist)
for (i in 1:nrow(tmp)) {
	id1 <- tmp$sp1[i]
	id2 <- tmp$sp2[i]
	islanddist_cod[id1, id2] <- as.numeric(tmp$island_bg[i])-1
	islanddist_cod[id2, id1] <- as.numeric(tmp$island_bg[i])-1
}
diag(islanddist_cod) <- 0


# Multiple regression on distance matrix approach (Lichstein 2007 Plant Ecol.)
library(ecodist)

MRM(as.dist(convdist_reg) ~ as.dist(divdist_reg) + as.dist(phydist) + as.dist(fishdist) + as.dist(islanddist), mrank=T)
# fish (506.6, P = 0.002), islands (685.5, P = 0.002) significant, div (0.0060, P=0.285) and phy NS (20.0, P=0.177)
# R^2 = 0.39
# this means, on average, plunge-diving species have 507 more convergent substitutions (total?)

MRM(as.dist(convdist_cod) ~ as.dist(divdist_cod) + as.dist(phydist) + as.dist(fishdist) + as.dist(islanddist), mrank=T)
# fish (6676.1, P = 0.001), islands (5890.7, P = 0.003), div (0.15, p=0.001), and phy significant (181.0,p=0.003)
# R^2 = 0.72
```

## Detecting signatures of positive selection

### CODEML positive selection analyses

```r
setwd("~/uce-alcedinidae/paml")

# setwd("~/genomeR")
# system("git pull https://celiason:eC6Na2hEb1Yi@github.com/celiason/genomeR.git")
devtools::load_all('~/genomeR')

#-------------------------------------------------------------------------------
# GeMoMa genes without early stop codons (N = 11938 genes)
#-------------------------------------------------------------------------------
library(pbmcapply)

fafiles <- list.files("~/uce-alcedinidae/output/cds_FINAL", full=TRUE)
length(fafiles)  # 11938 genes

# Fit M0 models (took 2-3d)
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M0", phy="/home/FM/celiason/uce-alcedinidae/paml/kingtree_nolabels.phy")
})

## NB These next models depend on results from M0, so wait until the above finishes

# Fit M1a models (<1h)
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M1a", fixedbl=TRUE)
})

# Fit M2a models
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M2a", fixedbl=TRUE)
})

# Fit M7 models
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M7", fixedbl=TRUE, outpath="M7_results")
})

# Fit M8 models
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M8", fixedbl=TRUE, outpath="M8_results")
})

# Fit branch-site Anull (all same) - BS1 model (?)
# runCodeml(fasta=fafiles[2], model="MAnull", fixedbl=TRUE, outpath="MAnull_results")
pbmclapply(fafiles, mc.cores=24, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="MAnull", fixedbl=TRUE, outpath="MAnull_results")
})

# Fit branch-site A (plunge-divers different) - BS2 model (?)
pbmclapply(fafiles, mc.cores=30, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="MA", fixedbl=TRUE, outpath="MA_results")
})

#-------------------------------------------------------------------------------
# Compare models using ML
#-------------------------------------------------------------------------------

setwd("~/uce-alcedinidae/paml")

library(pbmcapply)
library(stringr)

list.files("M1a_results")

# logLiks of M1 model
files <- paste0("M1a_results/", list.files("M1a_results", pattern="_"), "/M1a_mlc")
x1 <- parse_codeml(files=files, ncores=24)

# logLiks of M2 model
files <- paste0("M2a_results/", list.files("M2a_results", pattern="_"), "/M2a_mlc")
x2 <- parse_codeml(files=files, ncores=24)

# logLiks of M7 model
files <- paste0("M7_results/", list.files("M7_results", pattern="_"), "/mlc")
x3 <- parse_codeml(files=files, ncores=24)

# logLiks of M8 model
files <- paste0("M8_results/", list.files("M8_results", pattern="_"), "/mlc")
x4 <- parse_codeml(files=files, ncores=24)

# logLiks of BS1 model (=MA_null)
files <- paste0("MAnull_results/", list.files("MAnull_results", pattern="_"), "/mlc")
x5 <- parse_codeml(files=files, ncores=24)

# logLiks of BS2 model (=MA)
files <- paste0("MA_results/", list.files("MA_results", pattern="_"), "/mlc")
x6 <- parse_codeml(files=files, ncores=24)

allres <- data.frame(file=x1$file, ll_M1=x1$lnL, ll_M2=x2$lnL, ll_M7=x3$lnL, ll_M8=x4$lnL, ll_BS1=x5$lnL, ll_BS2=x6$lnL)

# Calculate P values
allres$P_M2 <- 1 - pchisq(2 * (allres$ll_M2 - allres$ll_M1), df=2)
allres$P_M8 <- 1 - pchisq(2 * (allres$ll_M8 - allres$ll_M7), df=2)
allres$P_BS2 <- 1 - pchisq(2 * (allres$ll_BS2 - allres$ll_BS1), df=1)

# Adjust P values with false discovery rate approach
allres$FDR_M2 <- p.adjust(allres$P_M2, method="fdr")
allres$FDR_M8 <- p.adjust(allres$P_M8, method="fdr")
allres$FDR_BS2 <- p.adjust(allres$P_BS2, method="fdr")

sum(allres$FDR_M2 < .05, na.rm=T)  # 3311 M2-M1
sum(allres$FDR_M8 < .05, na.rm=T)  # 8201 M8-M7
sum(allres$FDR_BS2 < .05, na.rm=T)  # 1446 (now 1427..?) BS2-BS1

tmp <- read.csv("~/uce-alcedinidae/output/gemoma_complete_overlapping.csv")

allres$complete <- ifelse(allres$file %in% gsub(":", "_", tmp$x), "complete", "fragment")

sum(allres$complete=="complete" & allres$FDR_M8 < .05, na.rm=T)  # 3155
sum(allres$complete=="complete" & allres$FDR_BS2 < .05, na.rm=T)  # 396

#-------------------------------------------------------------------------------
# Load functional annotations
#-------------------------------------------------------------------------------

mmseq <- read.table("~/uce-alcedinidae/annotations/gemoma/round2/uniprot.renamed.mmseqs")
mmseq2 <- read.table("~/uce-alcedinidae/annotations/gemoma/round2/chicken.mmseqs")
mmseq2$V1 <- gsub(":", "_", mmseq2$V1)

# Save unique uniprot IDs
# mylist <- unique(c(as.character(mmseq$V2),as.character(mmseq2$V2)))
# cat(mylist, sep="\n", file="uniprotIDs.txt")

# Lookup tables for converting gene IDs
map <- read.delim("uniprot.map")
m <- setNames(map$To[match(mmseq$V2, map$From)], mmseq$V1)
lookup <- read.table("~/uce-alcedinidae/annotations/gemoma/round2/filtered_predictions.map")
lookup$V1 <- gsub(":", "_", lookup$V1)

#-------------------------------------------------------------------------------
# Add gene data to output files
#-------------------------------------------------------------------------------

setdiff(allres$file, lookup$V1)

allres$geneID <- lookup$V2[match(allres$file, lookup$V1)]
allres$uni_gal <- mmseq2$V2[match(allres$file, mmseq2$V1)]
allres$gene <- toupper(m[as.character(allres$geneID)])
allres$gene_galGal <- as.character(map$To[match(allres$uni_gal, map$From)])
allres$finalgene <- toupper(ifelse(is.na(allres$gene), allres$gene_galGal, allres$gene))

sum(is.na(allres$finalgene))  # 1435 out of 11938 genes (12%) are NAs


#-------------------------------------------------------------------------------
# Output results
#-------------------------------------------------------------------------------

# write.csv(res, file="pos_selection_results_repeatmasked_20221209.csv")
write.csv(allres, file="~/uce-alcedinidae/output/allres_paml.csv")

# Get proportion of genes 1) annotated/not vs. 2) selected/not:
table(selec=allres$FDR_M2<.05, annot=!is.na(allres$finalgene))  # 8201 M2-M1
table(selec=allres$P_M7M8<.05, annot=!is.na(allres$finalgene))  # 8201 M8-M7
table(selec=allres$P_MAMA0<.05, annot=!is.na(allres$finalgene))  # 8201 MA-MAnull

#-------------------------------------------------------------------------------
# Export genes under pos selection for enrichment analyses
#-------------------------------------------------------------------------------

# Load positive selection results + annotations
# res <- read.csv(file="pos_selection_results_repeatmasked_20210504.csv")
# res <- read.csv(file="pos_selection_results_repeatmasked_20221209.csv")
allres <- read.csv(file="~/uce-alcedinidae/output/allres_paml.csv", row=1)
head(allres)

# outputting genes that don't have names- for manual
# cat(as.character(na.omit(res[is.na(res$finalgene) & res$P_M1M2<0.05, "file"])), sep="\n")

# cat(as.character(na.omit(res[is.na(res$finalgene) & res$P_MAMA0 < 0.05, "file"])), sep="\n")

#-------------------------------------------------------------------------------
# Output gene lists
#-------------------------------------------------------------------------------

# All genes
cat(sort(as.character(na.omit(allres$finalgene[res$FDR_M2<.05]))), file="genelist_M1M2.txt", sep="\n")
cat(sort(as.character(na.omit(allres$finalgene[res$FDR_M8<.05]))), file="genelist_M7M8.txt", sep="\n")
cat(sort(as.character(na.omit(allres$finalgene[res$FDR_BS2<.05]))), file="genelist_MAMA0.txt", sep="\n")

# Complete genes only
cat(sort(as.character(na.omit(allres$finalgene[allres$FDR_M2<.05&allres$complete=="complete"]))), file="genelist_M1M2_complete.txt", sep="\n")
cat(sort(as.character(na.omit(allres$finalgene[allres$FDR_M8<.05&allres$complete=="complete"]))), file="genelist_M7M8_complete.txt", sep="\n")
cat(sort(as.character(na.omit(allres$finalgene[allres$FDR_BS2<.05&allres$complete=="complete"]))), file="genelist_MAMA0_complete.txt", sep="\n")

# cat(unique(df$ipr), sep="\n", file="genelist_iprscan.txt")
cat(unique(as.character(allres$finalgene)), file="gene_background.txt", sep="\n")
# cat(na.omit(as.character(df$uni_gal)), sep="\n", file="unigal.txt")

table(is.na(allres$FDR_M2))[2]  # 1095 didn't run in PAML- mostly due to too much missing data (Ns) in alignments

#-------------------------------------------------------------------------------
# Get ancestral states
#-------------------------------------------------------------------------------

anc <- readLines("~/uce-alcedinidae/paml/M0_results/aquChr_transcript:ENSACCT00020000054_R0/rst")
phy <- read.tree("~/uce-alcedinidae/paml/kingtree_nolabels.phy")
f <- list.files("paml/M0_results",recursive=T,pattern="rst",full=TRUE)[1]

start <- grep("List of extant", anc)[1]+4
end <- grep("Overall accuracy of", anc)[1]-3
newanc <- t(do.call(cbind, str_match_all(anc[start:end], "[A-Z]{3,3}")))
rownames(newanc) <- c(sort(phy$tip), (Ntip(phy)+1) : (Ntip(phy)+Nnode(phy)))
newanc <- t(apply(newanc, 1, function(x) {
	translate(s2c(paste0(x, collapse="")))	
}))
```

Subset of taste receptor gene annotations from exonerate

```r
fafiles <- list.files("~/uce-alcedinidae/output/cds_sensory", full=T)

#-------------------------------------------------------------------------------
# Fit models
#-------------------------------------------------------------------------------

# M0 models
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M0", phy="/home/FM/celiason/uce-alcedinidae/paml/kingtree_nolabels.phy", outpath="sensory_M0")
})

# M1a models
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M1a", fixedbl=TRUE, outpath="sensory_M1a", M0path="sensory_M0", force=TRUE)
})

# M2a models
pbmclapply(fafiles, mc.cores=48, mc.preschedule=FALSE, function(f) {
	runCodeml(fasta=f, model="M2a", fixedbl=TRUE, outpath="sensory_M2a", M0path="sensory_M0")
})

#-------------------------------------------------------------------------------
# Compare models w/ML
#-------------------------------------------------------------------------------
resm1a <- list.files("sensory_M1a", pattern="mlc", recursive=TRUE, full=TRUE)
resm1a <- resm1a[!grepl("bad", resm1a)]  # remove ones in "bad" folder
resm2a <- list.files("sensory_M2a", pattern="mlc", recursive=TRUE, full=TRUE)
resm2a <- resm2a[!grepl("bad", resm2a)]  # remove ones in "bad" folder
ll.1a <- pbapply::pblapply(resm1a, function(x) {
	try(raw <- readLines(x))
	res <- stringr::str_match(raw[grep("lnL", raw)], "lnL.*?np.*?(\\d+)\\)\\:\\s+([\\.\\-\\d]+)")[,2:3]
	try(setNames(res, c('np', 'lnL')))
})
bad <- grep("Error", ll.1a)
ll.1a[bad] <- c(NA, NA)
x1 <- cbind(file=basename(dirname(resm1a)), do.call(rbind, ll.1a))
x1 <- as.data.frame(x1)
ll.2a <- pbmclapply(resm2a, mc.cores=24, function(x) {
	try(raw <- readLines(x))
	res <- str_match(raw[grep("lnL", raw)], "lnL.*?np.*?(\\d+)\\)\\:\\s+([\\.\\-\\d]+)")[,2:3]
	try(setNames(res, c('np', 'lnL')))
})
bad <- grep("Error", ll.2a)
ll.2a[bad] <- c(NA, NA)
x2 <- cbind(file=basename(dirname(resm2a)), do.call(rbind, ll.2a))
x2 <- as.data.frame(x2)

x1$lnL <- as.numeric(as.character(x1$lnL))
x2$lnL <- as.numeric(as.character(x2$lnL))

pvals <- 1-pchisq(2*(x2$lnL-x1$lnL), df=2)
pvals <- p.adjust(pvals, method="fdr")
table(pvals<.05)

#-------------------------------------------------------------------------------
# Calculate log likelihoods
#-------------------------------------------------------------------------------
res <- data.frame(file=x1$file, M1a_lnL = x1$lnL, M2a_lnL = x2$lnL, df=2, l=-2*(x1$lnL-x2$lnL), p=pvals)
res$file <- toupper(res$file)
subset(res, file %in% c("TAS2R4", "TAS1R1", "SCNN1A"))

#-------------------------------------------------------------------------------
# Output for RERconverge (on phoebe)
#-------------------------------------------------------------------------------
setwd("~/uce-alcedinidae/paml")

# On server-
files <- list.files("sensory_M0", pattern="M0_mlc", recursive=TRUE, full=TRUE)

# Takes ~6 mins for this to run-
trees <- sapply(files, function(f) {
	raw <- readLines(f)
	raw[grep("tree length", raw)[1]+4]
})
trees2 <- trees[!is.na(trees)]
trees2 <- trees2[!grepl("nan", trees2)]

# Table output for RERconverge
nms <- gsub("sensory_M0/", "", dirname(names(trees2)))

res <- cbind(nms[nms %in% c('Calhm1','Pkd2l1','Scnn1a','Scnn1b','Scnn1g','Tas1r1','Tas1r3','Tas2r4','Tas2r40','Trpm5')],
			 trees2[nms %in% c('Calhm1','Pkd2l1','Scnn1a','Scnn1b','Scnn1g','Tas1r1','Tas1r3','Tas2r4','Tas2r40','Trpm5')])

# nms <- str_extract(basename(names(raw)), "^.*?(?=\\.)")
write.table(res, file="taste_M0_genetrees.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=F)
```

### HYPHY aBS-REL analyses

```sh
# testing out with 1 run
# id=aquChr_transcript_ENSACCT00020000032_R5
# aln=~/uce-alcedinidae/output/cds_FINAL/$id.fa
# grep "((((((actHom" ~/uce-alcedinidae/paml/M0_results/$id/M0_mlc > tmp/$id.phy
# hyphy absrel --alignment $aln --code Universal --tree tmp/$id.phy --branches All --output hyphy_redo/$id.json
```

<!-- run all in parallel within R framework- -->

```r
setwd("~/uce-alcedinidae")

library(parallel)
library(pbapply)
library(pbmcapply)
library(jsonlite)
library(ape)

dir.create("hyphy/aBSREL")

files <- list.files("output/cds_FINAL", pattern=".fa")

ids <- gsub(".fa", "", basename(files))

mclapply(ids, mc.cores=24, function(x) {
	system(paste0("grep '((((((actHom' /home/FM/celiason/uce-alcedinidae/paml/M0_results/", x, "/M0_mlc > tmp/", x, ".phy"))
	run <- paste0("hyphy CPU=2 absrel --alignment /home/FM/celiason/uce-alcedinidae/output/cds_FINAL/", x, ".fa --code Universal --tree tmp/", x, ".phy --branches All --output hyphy/aBSREL/", x, ".json")
	system(run)
})

# Parse results
files <- list.files("hyphy/aBSREL", full=T)
res <- pbmclapply(files, mc.cores=12, function(x) try(read_json(x)))
pvals <- pbmclapply(seq_along(res), mc.cores=12, function(x) {
	try(sapply(res[[x]][["branch attributes"]][[1]], "[[", "Corrected P-value"))
})
names(pvals) <- gsub(".json", "", basename(files))
# Extract data
tmp <- do.call(rbind, pvals)
tmp <- apply(tmp, 2, as.numeric)
tmp <- as.data.frame(cbind(file=names(pvals), tmp))
write.csv(tmp, file="~/uce-alcedinidae/output/allres_absrel.csv")

# Parse results (uncorrected P-values)
files <- list.files("hyphy/aBSREL", full=T)
res <- pbmclapply(files, mc.cores=12, function(x) try(read_json(x)))
pvals <- pbmclapply(seq_along(res), mc.cores=12, function(x) {
	try(sapply(res[[x]][["branch attributes"]][[1]], "[[", "Uncorrected P-value"))
})
names(pvals) <- gsub(".json", "", basename(files))
# Extract data
tmp <- do.call(rbind, pvals)
tmp <- apply(tmp, 2, as.numeric)
tmp <- as.data.frame(cbind(file=names(pvals), tmp))
write.csv(tmp, file="~/uce-alcedinidae/output/allres_absrel_uncorrected.csv")
```

<!-- on local machine -->

```r
library(viridis)

tr <- ape::read.tree(text="(((((actHom:0.017832,actLin:0.019047)Node5:0.058249,((todChl:4e-06,todCin:4e-06)Node9:0.017672,todWin:0.018488)Node8:0.009789000000000001)Node4:0.011322,(((halLeu:0.029631,halSmy:0.038268)Node15:0.015679,halChe:0.031444)Node14:0.036666,(halBad:0.033486,halMal:0.042028)Node19:4e-06)Node13:0.009166000000000001)Node3:0.014156,((cerRud:0.035055,(chlInd:0.045467,chlAma:0.02306)Node25:0.024884)Node23:0.001019,((megAlc:0.012937,megTor:0.005506)Node29:0.008935999999999999,megMax:0.019787)Node28:0.029843)Node22:0.026609)Node2:0.003129,(((alcAtt:0.008836,(alcQua:0.019544,alcSem:0.050167)Node37:0.012772)Node35:0.031666,((((ceyArg:0.012803,ceyCya:0.034852)Node43:0.002585,ceyMar:4e-06)Node42:0.045643,ceyMel:0.055676)Node41:0.00151,ceyEri:0.080299)Node40:4e-06)Node34:4e-06,((((corCri:0.114857,corVin:0.007898000000000001)Node52:4e-06,corLeu:0.070631)Node51:4e-06,corMad:0.005392)Node50:0.006482,(ispLec:0.021555,ispPic:0.054634)Node57:0.0067)Node49:0.0483)Node33:0.152572,momMom:0.007462);")

# save(tr, tmp, file="hyphy/absrel_results.rda")
# load("hyphy/absrel_results.rda")

tmp <- read.csv(file="output/allres_merged.csv", row=1)
tmp <- tmp[, c(1, 22:80)]
names(tmp) <- gsub("P_aBSREL_", "", names(tmp))
rownames(tmp) <- tmp[,1]
tmp <- tmp[,-1]
head(tmp)

tr2 <- read.tree("data/paml/kingTree_abbrev.phy")
tr2 <- makeNodeLabel(tr2)
tr2.pruned <- drop.tip(tr2, "momMom")


# par(mfrow=c(1,2))
# plot(tr, show.node.label=T, cex=.5, no.margin=T)
# plot(tr2, show.node.label=T, cex=.5)

m <- read.csv("data/absrel_match_nodes.csv")
m[,1] <- paste0("Node", m[,1])
m[,2] <- paste0("Node", m[,2])

nms <- c(tr2.pruned$tip.label, tr2.pruned$node.label)
idx <- tr2.pruned$edge[,2]

edgenames <- nms[idx]
# edgenames[1] <- "Node1"
edgenames[match(m[,2], edgenames)] <- m[,1]

# Calculate number of positively selected genes per branch
num_pos <- apply(tmp, 2, function(x) sum(x<.05, na.rm=T))
num_pos_scaled <- minmax(num_pos[edgenames], na.rm=T)

# Setup color of branches
# make_colors <- colorRamp(c('blue','lightgray','red'))
make_colors <- colorRamp(rocket(10))
edgecols <- rgb(make_colors(num_pos_scaled), maxColorValue=255)

# edgecols <- c("gray", edgecols)
names(edgecols) <- edgenames

# Plot tree for figure
pdf("figs/absrel_tree.pdf", width=5.5, height=8)
par(mar=c(0,0,0,0))
plot(tr2.pruned, edge.col=edgecols, edge.width=3, label.offset=.5)
tiplabels(pch=22, bg=ifelse(tr2$tip.label %in% plungeFg, "black", "white"))
legend_image <- as.raster(matrix(rev(rocket(10)), ncol=1))
rasterImage(legend_image, 0, 26, 1, 30)
rect(0, 26, 1, 30)
lbsq <- seq.int(0, 4, l=5)  # seq. for labels
axis(4, at=lbsq+26, pos=1, labels=F, col=0, col.ticks=1, tck=-.01)  ## axis ticks
mtext(round(seq.int(10, 348, l=5)), side=2, line=-3.5, at=lbsq+26, las=2, cex=.6)
dev.off()

tmp2 <- ifelse(tmp<.05, 1, 0)
tmp2 <- tmp2[complete.cases(tmp2), ]
tmp2 <- tmp2[, !grepl("Node", colnames(tmp2))]
dmat <- dist(t(tmp2))
# apropos("graph_from")
# graph_from_adjacency_matrix(as.matrix(dmat))
plot(hclust(dmat))
obj = cophylo(tr2, as.phylo(hclust(dmat)))
plot(obj)

# taking too long so stopped it..
# pca <- prcomp(t(tmp2))
```

### HYPHY RELAX analyses

<!-- We're also going to validate the BS2 models with RELAX -->

```r
setwd("~/Storage/celiason_Phoebe")

library(parallel)
library(jsonlite)

fafiles <- list.files("cds_FINAL", pattern=".fa$", full=TRUE)
length(fafiles)  # 11938

annot <- read.csv("~/uce-alcedinidae/paml/pos_selection_results_repeatmasked_20221209_absrel.csv", row=1)

picks <- annot$file[annot$P_MAMA0<.05]
picks <- picks[!is.na(picks)]

# grep("19941", picks, value=T) # mapt to check..

# sanity check
all(file.exists(paste0("cds_FINAL/", picks, ".fa")))

# dir.create("temptrees")
dir.create("RELAX_trees")
dir.create("RELAX_output")

#-------------------------------------------------------------------------------
# Fit RELAX models in parallel
#-------------------------------------------------------------------------------
mclapply(picks, mc.cores=10, function(x) {
	system(paste0("grep '((((((actHom' ./paml_M0_results/", x, "/M0_mlc > RELAX_trees/", x, ".phy"))
	# load tree with branch lengths from paml M0 model fit
	edges <- ape::read.tree(paste0("RELAX_trees/", x, ".phy"))$edge.length
	# load foreground labeled tree
	tr <- ape::read.tree("mapt.phy")
	# replace edges with PAML-inferred edges
	tr$edge.length <- edges
	# re-write tree for RELAX
	ape::write.tree(tr, file=paste0("RELAX_trees/", x, ".phy"))
	# run model
	run <- paste0("hyphy CPU=4 relax --alignment cds_FINAL/", x, ".fa --code Universal --tree RELAX_trees/", x, ".phy --test Foreground --output RELAX_output/", x, ".json")
	system(run)
})

#-------------------------------------------------------------------------------
# Add results to positive selection output
#-------------------------------------------------------------------------------

files <- paste0("RELAX_output/", picks, ".json")
out <- lapply(files, function(x) try(read_json(x)))
pvals <- sapply(out, function(x) try(x[["test results"]][["p-value"]]))
pvals <- as.numeric(pvals)
names(pvals) <- picks
# padj <- p.adjust(pvals, method="fdr")
kvals <- sapply(out, function(x) try(x[["test results"]][["relaxation or intensification parameter"]]))
kvals <- as.numeric(kvals)
names(kvals) <- picks

df <- data.frame(file=names(kvals), k=kvals, p=pvals, row.names=NULL)
write.csv(df, "~/uce-alcedinidae/output/allres_relax.csv")
```

on personal computer

```r
meta <- read.csv(file="~/Downloads/pos_selection_results_repeatmasked_20221209_relax.csv", row=1)
head(meta)
hist(meta$Krelax, breaks=50)
plot(meta$Krelax[meta$Prelax < .15], log='y')

plot(Krelax ~ Prelax, meta, log='y')
plot(Prelax ~ P_MAMA0, meta, log='xy')
abline(v=.15)
abline(h=.15)
range(meta$P_MAMA0, na.rm=T)

table(bs2=meta$P_MAMA0 < .15, rel=meta$Prelax < .15)

out <- meta$finalgene[meta$P_MAMA0 < .15 & meta$Prelax < .15 & meta$Krelax > 1]
out <- sort(unique(out[!is.na(out)]))
out
clipr::write_clip(out)  # writing genes to analyze in gprofile2, string networks
```

## Testing for convergent amino acid substitutions

To test whether plunge-diving kingfishers had similar protein sequences due to convergent evolution, we used the program PCOC (33). This model tests for signatures of convergence in all residues of a protein. We set the cutoff at 0.8, or the posterior probability that a given residue is evolving convergently in plunge-diving lineages. Genes that showed support for convergence were retained and used in downstream enrichment analyses.

```sh
# All genes
cd ~/uce-alcedinidae/PCOC
CMD_PCOC_DOCKER="docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD -e CWD=$PWD carinerey/pcoc"
# echo $CMD_PCOC_DOCKER

# find output/cds_gemoma/*.fa | head -n1
# faTrans output/cds_gemoma/aquChr_transcript:ENSACCT00020000054_R0.fa test.fa

#-------------------------------------------------------------------------------
# A. Prepare tree for detection analysis: set putative convergent leaves/nodes
#-------------------------------------------------------------------------------

# i. The 1st step is to "number" the tree: each node in the input tree is labeled and given a number.
$CMD_PCOC_DOCKER pcoc_num_tree.py -t kingtree_nolabels.phy -o num_tree.pdf ##takes a few seconds; replace file name of tree.nw with your input tree

# ii. In the 2nd step, you use the numbered tree (output of the 1st step) to write a string corresponding to the convergent scenario you are interested in.

# iii. (Optional) The 3rd step is to visualize the final tree with colored labeled branches (to ascertain that no branch has been forgotten during 2nd step)
# this takes a few seconds
scenario="29,23,19,22,20,21,28,26,24,25,27/35,31,34,32,33/38,36,37/50,48,46,47,49"

$CMD_PCOC_DOCKER pcoc_num_tree.py -t kingtree_nolabels.phy -o num_tree.pdf -m $scenario 


# A. Translate into amino acids
for f in `find ../output/cds_FINAL`
do
	outfile=`basename $f .fa`.aa.fa; faTrans $f aa/$outfile
done

# B. Run detection analysis
for f in `find aa/*.fa`
do
	$CMD_PCOC_DOCKER pcoc_det.py -t kingtree_nolabels.phy -aa $f -o output_pcoc_det -m $scenario -f 0.8	
done

#------------------------------------------------------------------------------
# Just sensory genes (N=122)
#------------------------------------------------------------------------------

cd ~/uce-alcedinidae/PCOC/sensory

CMD_PCOC_DOCKER="docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v $PWD:$PWD -e CWD=$PWD carinerey/pcoc"

# A. Prepare tree for detection analysis: set putative convergent leaves/nodes
$CMD_PCOC_DOCKER pcoc_num_tree.py -t kingtree_nolabels.phy -o num_tree.pdf ##takes a few seconds; replace file name of tree.nw with your input tree
scenario="29,23,19,22,20,21,28,26,24,25,27/35,31,34,32,33/38,36,37/50,48,46,47,49"
$CMD_PCOC_DOCKER pcoc_num_tree.py -t kingtree_nolabels.phy -o num_tree.pdf -m $scenario 

# B. Translate sequences
for f in `find ../../output/cds_sensory/*.fa`
do
	outfile=`basename $f .fa`.aa.fa; faTrans $f aa/$outfile
done

# C. Run detection analysis
for f in `find aa/*.fa`
do
	$CMD_PCOC_DOCKER pcoc_det.py -t kingtree_nolabels.phy -aa $f -o output_pcoc_det -m $scenario -f 0.8	
done

# Pkd2l1 -  2/849 <1%
# Calhm1 -  3/321 <1%
# Scnn1a -  1/387 <1%
# Scnn1b -    N/A
# Scnn1g -  1/651 <1%
# Tas1r1 - 10/812 1.2%
# Tas1r3 - 14/835 1.6%
# Tas2r40 - 1/283 <1%
# Tas2r4  - 1/311 <1%
# Tas2r7  - 0/306 0%
# Trpm5   - 9/1169 <1%

# samtools tview -p 242:3360555 alignments/alcAtt-to-todChl.bam genomes/todChl/todChl.fasta

#-------------------------------------------------------------------------------
# 4. Graphical output-
#-------------------------------------------------------------------------------
$CMD_PCOC_DOCKER pcoc_det.py -t kingtree_nolabels.phy -aa aa/Tas1r1.aa.fa -o output_pcoc_det -m $scenario -f_pcoc 0.8 --plot --reorder
$CMD_PCOC_DOCKER pcoc_det.py -t kingtree_nolabels.phy -aa aa/Tas1r3.aa.fa -o output_pcoc_det -m $scenario -f_pcoc 0.8 --plot --reorder
```

<!-- Get genes convergently evolved -->

```r
setwd("~/uce-alcedinidae")

# Load set of genes with some sites convergent (i.e., has "filtered" in the name)
files <- list.files("PCOC/output_pcoc_det/", recursive=TRUE, pattern="filtered", full=TRUE)
length(files)  # 6536 genes

# TAS1R3 (=galGal_rna-XM_025142692.1_R0)
read.delim(grep("galGal_rna-XM_025142692.1_R0", files, value=TRUE))

# MAPT gene
read.delim(grep("aquChr_transcript:ENSACCT00020019941_R0", files, value=TRUE))

nms <- gsub(".aa.filtered.*?$", "", basename(files))
nms <- gsub(":", "_", nms)

# Append to pos selection results CSV-
# res <- read.csv(file="~/uce-alcedinidae/pos_selection_results_repeatmasked_20221209_absrel_relax.csv", row=1)

df <- data.frame(file=nms, PCOC="conv")

# res$P_PCOC <- "not_conv"
# res$P_PCOC[match(nms, res$file)] <- "conv"

write.csv(df, file="~/uce-alcedinidae/output/allres_pcoc.csv")

# We ran on all 11938 genes, so OK to make others "not_conv"
# allfiles <- list.files("../PCOC/output_pcoc_det/", pattern="aa.results.tsv", recursive=TRUE, full=TRUE)
# table(res$P_PCOC)

# Output gene list
# cat(sort(unique(res$finalgene[res$P_PCOC=="conv"])), sep="\n", file="../PCOC/pcoc_filtered_20221209.txt")
```

## Create final data set for gene enrichment analyses

```r
allres <- read.csv("~/uce-alcedinidae/output/allres_paml.csv", row=1)

absrel <- read.csv("~/uce-alcedinidae/output/allres_absrel.csv", row=1)
colnames(absrel)[-1] <- paste0("P_aBSREL_", colnames(absrel)[-1])
anysig <- apply(absrel[,-1], 1, function(x) any(x<.05))
absrel$SigAny_absrel <- ifelse(anysig, "yes", "no")
table(absrel$SigAny_absrel)  # 2210 significant

relax <- read.csv("~/uce-alcedinidae/output/allres_relax.csv", row=1)
names(relax) <- plyr::revalue(names(relax), c('k'='k_relax', 'p'='P_relax'))

pcoc <- read.csv("~/uce-alcedinidae/output/allres_pcoc.csv", row=1)

df <- dplyr::left_join(allres, relax)
df <- dplyr::left_join(df, absrel)
df <- dplyr::left_join(df, pcoc)
df$PCOC <- factor(df$PCOC, levels=c("not_conv", "conv"))
df$PCOC[is.na(df$PCOC)] <- "not_conv"
table(df$PCOC)

# df <- dplyr::left_join(df, absrel)

# res$Prelax <- NA
# res$Prelax[match(names(padj), res$file)] <- padj
# res$Krelax <- NA
# res$Krelax[match(names(kvals), res$file)] <- kvals
# head(res)

# table(res$Prelax<.05,res$P_aBSREL<.05)
# res <- res[, names(res)!="X"]
# write.csv(res, file="~/uce-alcedinidae/pos_selection_results_repeatmasked_20221209_absrel_relax.csv")

# absrel
# res <- read.csv("paml/pos_selection_results_repeatmasked_20221209.csv", row.names=1)

# res$P_aBSREL <- NA
# res$P_aBSREL[match(names(which(anysig)), res$file)] <- 0
# res$P_aBSREL[match(names(which(!anysig)), res$file)] <- 1

write.csv(df, "~/uce-alcedinidae/output/allres_merged.csv")
```

## Gene enrichment

```r
d <- data.frame(node=1:(Nedge(tree)+1), color="gray")
d$color[phytools::getDescendants(tree, getMRCA(tree, c('chlInd','megAlc')))] <- "blue"
d$color[phytools::getDescendants(tree, getMRCA(tree, c('alcAtt','alcSem')))] <- "blue"
d$color[phytools::getDescendants(tree, getMRCA(tree, c('ceyArg','ceyCya')))] <- "blue"
d$color[phytools::getDescendants(tree, getMRCA(tree, c('corCri','corLeu')))] <- "blue"
p <- ggtree(tree) %<+% d + aes(color=I(color)) + geom_tiplab()

# first plot some alignments next to the tree to show convergent sites:

# MAPT gene
x <- readAAStringSet("~/Downloads/aquChr_transcript:ENSACCT00020019941_R0.aa.fa")
data = tidy_msa(x, start=202, end=220)
p + geom_facet(geom = geom_msa, data = data, panel = 'msa',
               font = "helvetical", color = "Chemistry_AA", use_dot=TRUE, posHighligthed=c(207,209,215)) +
	# scale_alpha_identity() +
    xlim_tree(50) +
    theme(legend.position="none", strip.background=element_blank(), strip.text=element_blank())
ggsave("figs/mapt_with_tree.pdf", width=7, height=5)

# TAS1R3 gene
x <- readAAStringSet("data/pcoc_analyses/aln/TAS1R3_AA.fa")
xsub <- AAStringSet(sapply(x, "[", c(115:117,257:259,484:486,498:500,577:579,738:740)))
data = tidy_msa(xsub)
p + geom_facet(geom = geom_msa, data = data, panel = 'msa',
               font = "helvetical", color = "Chemistry_AA", use_dot=TRUE, posHighligthed=seq(2,17,by=3)) +
    xlim_tree(50) +
    theme(legend.position="none", strip.background=element_blank(), strip.text=element_blank())
ggsave("figs/tas1r3_with_tree.pdf", width=7, height=5)

# Load sets of positively selected/convergent genes
sets.list <- list(M2=readLines("output/genelist_M1M2_curated.txt"),
			MA=readLines("output/genelist_MAMA0_curated.txt"),
			M8=readLines("output/genelist_M7M8.txt"),
			PCOC=readLines("output/pcoc_filtered_20210805.txt"))
sets.ids <- unique(unlist(sets.list))


res_old <- read.csv("output/pos_selection_results_repeatmasked_20221209_absrel_relax_pcoc.csv", row=1)
head(res_old)

picks <- setdiff(res_old$file[which(res_old$P_MAMA0<.05)], allres$file[which(allres$FDR_BS2<.05)])

subset(allres, file %in% picks)

allres <- read.csv("output/allres_merged.csv", row=1)
head(allres)

sets.list <- list(
		M2 = sort(unique(na.omit(allres$geneID[allres$FDR_M2 < 0.05]))),
		BS2 = sort(unique(na.omit(allres$geneID[allres$FDR_BS2 < 0.05]))),
		M8 = sort(unique(na.omit(allres$geneID[allres$FDR_M8 < 0.05]))),
		PCOC = sort(unique(na.omit(allres$geneID[allres$PCOC == "conv"]))),
		# RELAX = sort(unique(na.omit(allres$geneID[allres$Prelax<.05]))),
		"aBS-REL" = sort(unique(na.omit(allres$geneID[allres$SigAny_absrel == "yes"])))
		)
sets.ids <- unique(unlist(sets.list))

dat <- sapply(sets.list, function(x) table(factor(x, levels = sets.ids)))
dat <- as.data.frame(dat)
dim(dat)
# colnames(dat) <- c('M2', 'MA', 'M8', 'PCOC')

sum(apply(dat, 1, function(x) all(x==1)))  # 585 overlapping genes
apply(dat, 2, sum)

plot(dat_paml$k_relax)


# pick absrel genes selected in at least 2 convergent lineages
clades <- list(c('corLeu','corVin','corCri'), c('ceyCya','ceyArg'), c('alcSem','alcQua','alcAtt'), c('megMax','megTor','megAlc','chlAma','chlInd','cerRud'))
subdat <- lapply(clades, function(x) {
	picks <- grep(paste0(x, collapse="|"), names(dat_paml))
	dat_paml[, picks]	
})
subdat_any <- sapply(subdat, function(x) which(apply(x, 1, function(x) any(x<.05))))
idx <- combn(1:length(subdat), m=2)
subdat_res <- sapply(1:ncol(idx), function(x) {
	Reduce(intersect, subdat_any[idx[, x]])
})
adap_conv_plunge <- sort(unique(dat_paml$finalgene[unlist(unique(subdat_res))]))  # used below-->

clades <- list('corMad', 'corVin', c('ceyEri','ceyMel','ceyMar','ceyCya'), 'todCin', c('actLin','actHom'))
subdat <- lapply(clades, function(x) {
	picks <- grep(paste0(x, collapse="|"), names(dat_paml))
	dat_paml[, picks, drop=F]
})
subdat_any <- sapply(subdat, function(x) which(apply(x, 1, function(x) any(x<0.05))))
idx <- combn(1:length(subdat), m=2)
subdat_res <- sapply(1:ncol(idx), function(x) {
	Reduce(intersect, subdat_any[idx[, x]])
})
adap_conv_island <- sort(unique(dat_paml$finalgene[unlist(unique(subdat_res))]))  # used below-->


sets_list_genes <- list(
	M2 = sort(unique(dat_paml$finalgene[which(dat_paml$FDR_M2<.05)])),
	# M2 = sort(unique((dat_paml$finalgene[dat_paml$FDR_M2 < 0.05]))),
	BS2 = sort(unique(dat_paml$finalgene[which(dat_paml$FDR_BS2<.05)])),
	# BS2 = sort(unique((dat_paml$finalgene[dat_paml$FDR_BS2 < 0.05]))),
	M8 = sort(unique(dat_paml$finalgene[which(dat_paml$FDR_M8<.05)])),
	# M8 = sort(unique(na.omit(dat_paml$finalgene[dat_paml$FDR_M8 < 0.05]))),
	PCOC = sort(unique(dat_paml$finalgene[which(dat_paml$PCOC == "conv")])),
	# PCOC = sort(unique(na.omit(dat_paml$finalgene[dat_paml$PCOC == "conv"]))),
	RELAX = sort(unique(dat_paml$finalgene[which(dat_paml$k_relax > 0.05 & dat_paml$P_relax<.05)])),
	# "aBSREL" = sort(unique(na.omit(dat_paml$finalgene[dat_paml$SigAny_absrel == "yes"])))
	aBSREL_any = sort(unique(dat_paml$finalgene[which(dat_paml$SigAny_absrel == "yes")])),
	aBSREL = adap_conv_plunge,
	aBSREL_isl = adap_conv_island
)

sets_ids_genes <- unique(unlist(sets_list_genes))

dat_genes <- sapply(sets_list_genes, function(x) table(factor(x, levels = sets_ids_genes)))
dat_genes <- as.data.frame(dat_genes)
head(dat_genes)

venn(dat_genes[c('M2','M8')])

names(dat_genes)

genelist_m2_notconv <- sort(unique(rownames(subset(dat_genes, BS2==0&M2==1&RELAX==0&PCOC==0))))
clipr::write_clip(genelist_m2_notconv)

plot(venneuler(dat_genes[c('BS2','PCOC','aBSREL','RELAX','M2')]))

df <- data.frame(aBSREL=dat_genes$aBSREL, M2=dat_genes$M2, PCOC=dat_genes$PCOC, BS2=ifelse(dat_genes$BS2==1&dat_genes$RELAX==1, 1, 0))
upset(df, nsets=ncol(df), order.by="freq")

pdf("figs/euler_m2.pdf", width=2, height=2)
par(mar=rep(0,4), ps=8)
venneuler(dat_genes[c('M2','M8')])
plot(venneuler(c(M2=25, `M2&M8`=2976, `M8`=4666), cex=1))
dev.off()



genelist_m2 <- names(which(apply(dat_genes[,c('M2','M8')],1,sum)==2))
clipr::write_clip(genelist_m2)

genelist <- names(which(apply(dat_genes[,c('M2','M8','aBSREL')],1,sum)==3))

genelist_adap_conv <- names(which(apply(dat_genes[,c('RELAX','PCOC','BS2','aBSREL')],1,sum)==4))
clipr::write_clip(genelist)

genelist_isl <- names(which(apply(dat_genes[,c('RELAX','PCOC','BS2','aBSREL_isl')],1,sum)==4))



picks <- grep(paste0(plungeFg, collapse="|"), names(dat_paml))
genes <- lapply(picks, function(x) {
	sort(unique(dat_paml$finalgene[which(dat_paml[, x] < 0.05)]))
})
names(genes) <- gsub("P_aBSREL_", "", names(dat_paml)[picks])
df <- matrix("", nrow=max(sapply(genes, length)), ncol=length(genes))
for (i in 1:length(genes)) {
	df[1:length(genes[[i]]), i] <- genes[[i]]
}
colnames(df) <- names(genes)
write.csv(df, file="output/mult_plunge_genes.csv", row.names=F)


genelist <- names(which(apply(dat_genes[,c('RELAX','PCOC','BS2')],1,sum)==3))
clipr::write_clip(genelist)

x1=list(intersect(sets_list_genes$BS2, sets_list_genes$RELAX))
library(VennDiagram)
library(gplots)

venn(c(RELAX=x1, sets_list_genes[c('PCOC','aBSREL')]))

limma::vennDiagram(c(RELAX=x1, sets_list_genes[c('PCOC','aBSREL')]))
library(venneuler)
dat_genes$BS2_filt <- ifelse(dat_genes$BS2==1 & dat_genes$RELAX==1, 1, 0)
df <- dat_genes[c('BS2_filt','aBSREL','PCOC')]
df <- df[!apply(df, 1, sum)==0, ]

sum(df[,1]==1 & df[,2]==1)
sum(df[,1]==1 & df[,2]==1)
venn(df)
plot(venneuler(df))

pdf("figs/euler.pdf", width=2, height=2)
par(mar=rep(0,4), ps=8)
plot(venneuler(c(BS2=25, PCOC=5554, `BS2&PCOC`=160, aBSREL=40, `aBSREL&PCOC`=284, `BS2&aBSREL`=1, `PCOC&aBSREL&BS2`=93)), cex=1)
dev.off()

genelist <- names(which(apply(dat_genes[,c('M2','M8','aBSREL','RELAX','PCOC','BS2')],1,sum)==6))
clipr::write_clip(genelist)

# KIF6 is involved in ATP binding and microtubule movement (cite honeycreeper paper)
# AGT is a dietary enzyme(!) that has been linked to diet in other birds (Wang et al. 2020)
# CENPJ role in brain size!? - a microcephaly gene in primates (montgomery et al. 2011)

#-------------------------------------------------------------------------------
# Upset plot showing overlapping sets
#-------------------------------------------------------------------------------
library(ComplexUpset)
library(ggplot2)

# conflicted::conflict_prefer("upset", "ComplexUpset")
conflicted::conflict_prefer("upset", "UpSetR")

# pdf("figs/upset_plot_new.pdf", width=6.5, height=4)
picks=c('PCOC','aBSREL','RELAX','BS2')
upset(dat_genes[picks], intersect=picks, width_ratio=0.05, height_ratio=.2,
	    set_sizes=
        upset_set_size()
        + theme(axis.text.x=element_text(angle=0, size=6))
    )
ggsave("figs/upset_plot_new.pdf", width=6.5*2, height=4*2)
# dev.off()

upset(dat_genes, intersect=c('M2','M8'), width_ratio=0.1, height_ratio=.2,
	    set_sizes=
        upset_set_size()
        + theme(axis.text.x=element_text(angle=45, size=6))
    )


picks <- c('PCOC','BS2','RELAX','aBSREL')
df = dat_genes[, picks]
df <- df[!apply(df, 1, sum)==0, ]

pdf("figs/upset_adap_conv.pdf", width=6.5, height=3)
upset(df, nsets=6, nintersects=20, order.by="freq", decreasing=T, text.scale=.85)
dev.off()

upset(df, intersect=names(df), width_ratio=0.1, height_ratio=.2,
	    set_sizes=
        upset_set_size()
        + theme(axis.text.x=element_text(angle=45, size=6))
    ) + theme(text=element_text(size=9))
ggsave("figs/upset_adap_conv.pdf", width=6.5, height=3)
table(bs2=allres$FDR_BS2<.05, relax=allres$P_relax<.05)  # 308 genes verified


#-------------------------------------------------------------------------------
# Convert gene names to GO terms
#-------------------------------------------------------------------------------
conv <- gconvert(rownames(dat), organism="ggallus", target="ENTREZGENE_ACC")
# clipr::write_clip(unique(conv$target))

#-------------------------------------------------------------------------------
# Run enrichment analyses of sets
#-------------------------------------------------------------------------------
# Order: M2 M8 MA PCOC aBS-REL (so, e.g., 01000 indicates M8 model, 11000 indicates M2 & M8 overlap, etc.)

gp_11000 <- gost(query = rownames(dat_genes)[dat_genes$PCOC==0 & dat_genes$BS2==0 & dat_genes$M2==1 & dat_genes$M8==1 & dat_genes$"aBS-REL"==0], organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")

# all genes
gp_11111 <- gost(query = rownames(dat_genes)[dat_genes$PCOC==1 & dat_genes$BS2==1 & dat_genes$M2==1 & dat_genes$M8==1 & dat_genes$aBSREL==1], organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")
pdf(file="figs/network_11111.pdf", width=8, height=4.5)
set.seed(1982)
emapplot(gost_to_enrich(gp_11111), showCategory=30, layout="dh", line_scale=.1) + scale_color_gradient(low='blue', high='lightgray', name="P value")
dev.off()

# just M8 model
gp_0100 <- gost(query = rownames(dat_genes)[dat_genes$M8==1], organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")


gp_0101 <- gost(query = rownames(dat_genes)[dat_genes$PCOC==1 & dat_genes$M8==1],
				organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")
gp_1100 <- gost(query = rownames(dat_genes)[dat_genes$M2==1 & dat_genes$M8==1 & dat_genes$MA==0 & dat_genes$PCOC==0], 
				organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")
gp_1101 <- gost(query = rownames(dat_genes)[dat_genes$PCOC==1 & dat_genes$M8==1 & dat_genes$M2==1],
				organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")
gp_0011 <- gost(query = rownames(dat_genes)[dat_genes$M2==0 & dat_genes$M8==0 & dat_genes$PCOC==1 & dat_genes$MA==1], 
				organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")
gp_1111 <- gost(query = rownames(dat_genes)[dat_genes$PCOC==1 & dat_genes$MA==1 & dat_genes$M2==1 & dat_genes$M8==1], 
				organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")

gp_0010 <- gost(query = rownames(dat_genes)[dat_genes$PCOC==0 & dat_genes$MA==1 & dat_genes$M2==0 & dat_genes$M8==0], 
				organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")
gp_0001 <- gost(query = rownames(dat_genes)[dat_genes$PCOC==1 & dat_genes$MA==0 & dat_genes$M2==0 & dat_genes$M8==0], 
				organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP")

#-------------------------------------------------------------------------------
# Make plots
#-------------------------------------------------------------------------------
pdf(file="figs/network_m2m8.pdf", width=7, height=6)
emapplot(gost_to_enrich(gp_1100), showCategory=10)
dev.off()


pdf(file="figs/network_11000.pdf", width=7, height=6)
emapplot(gost_to_enrich(gp_11000), showCategory=10, layout="kk")
dev.off()


pdf(file="figs/network_m8.pdf", width=7, height=6)
emapplot(gost_to_enrich(gp_0100), showCategory=10)
dev.off()

# PCOC + M8
pdf(file="figs/network_pcoc+m8.pdf", width=7, height=6)
emapplot(gost_to_enrich(gp_0101), showCategory=10)
dev.off()

# PCOC + M8 + M2
pdf(file="figs/network_pcoc+m8+m2.pdf", width=7, height=6)
emapplot(gost_to_enrich(gp_1101))
dev.off()

# PCOC + Ma
pdf(file="figs/network_pcoc+ma.pdf", width=7, height=6)
emapplot(gost_to_enrich(gp_0011), showCategory=10)
dev.off()

# All 4
pdf(file="figs/network_pcoc+ma+m8+m2.pdf", width=7, height=6)
emapplot(gost_to_enrich(gp_1111))
dev.off()

dat[c('TAS1R3', 'TAS2R40','TAS2R7','SCNN1B','SCNN1G','CALHM1','PKD2L1'),]
```

## Studying tongue shape and taste receptors

<!-- prepare data sets for statistical analysis -->

```r
setwd("~/Dropbox/Projects/king_tongues")

# setwd("C:/Users/celiason/Dropbox/Kingfisher taste receptors")
# setwd("/Users/Laur/Dropbox/Kingfisher")

dat <- readxl::read_xlsx("data/All.tongues.xlsx")
length(unique(dat$Label)) #42 good

labels <- split(dat, dat[,2])
names(labels) <- str_extract(names(labels), "[A-Z]+\\_\\d+")
names(labels) <- str_replace(names(labels), "BBM", "BPBM")
names(labels) <- str_replace(names(labels), "344953394", "344953")
names(labels) <- str_replace(names(labels), "472379", "427379")

#set resolution
reso <- read.csv("data/stacksizes.csv")
res <- reso$reso[match(names(labels), reso$X)]

tonguewidths <- numeric()
tonguelengths <- numeric()
XYZ <- list()

## going from 2,12 to 1 and then adding 1 to 22
## preferred??? this is what Lauren's stats were based upon
for (i in 1:length(labels)) {
  X <- labels[[i]]$X
  Y <- labels[[i]]$Y
  Z <- res[i] * labels[[i]]$Slice
  xyz <- cbind(X,Y,Z)
  dc <- dist(xyz[c(2,1),])   # hypoteneuse
  da <- dist(xyz[c(2,12), ])/2
  db <- ifelse(dc>da, sqrt(dc^2 - da^2), sqrt(da^2 - dc^2))
  tonguewidths[i] <- dist(xyz[c(2, 12), ])
  tonguelengths[i] <- (db + dist(xyz[c(1,22), ]))
  XYZ[[i]] <- xyz
}

#assigning label names to species name
shape <- tonguelengths/tonguewidths
Species <- readxl::read_xlsx("data/kingfisher_CT_scans_050817.xlsx")
names(shape) <- names(labels)
names(shape) <- Species$spp_andersen2017[match(names(shape), Species$ID)]
names(shape) <- gsub(" ", "_", names(shape))

data.frame(species=names(diet), biosample=names(labels), fish=diet, width.mm=round(tonguewidths,3), length.mm=round(tonguelengths,3), aspect=round(shape,2)) %>% write.csv("output/tongue_data_for_supp.csv", row.names=FALSE)
```

<!-- analyses -->

```r
setwd("~/Dropbox/Projects/King_genome")

dat <- read.csv("doc/ms_comparative/5-sciadv/SI/dataS1-tongue_data_for_supp.csv")

t.test(aspect ~ fish, data = dat)  # t = -2.45, pval = 0.024


# load tree
list.files("data")
tre <- read.tree("data/kingTree.BEAST.mean.age.newick copy.tre")
tre$tip.label <- gsub("Ceyx_melanurus_melanurus", "Ceyx_melanurus", tre$tip, fixed=TRUE)
plot(tre)

# add body mass
#diet data for 0=nonfish and 1=fish
DietData <- readxl::read_xlsx("data/kingdata_v2_edited.xlsx")
DietData$spp <- gsub("Ceyx melanurus melanurus", "Ceyx melanurus", DietData$sppAndersen, fixed=TRUE)
DietData$spp <- gsub(" ", "_", DietData$spp)

bmass <- setNames(DietData$weight_g, DietData$spp)
bmass <- bmass[!is.na(bmass)]

dat$logmass <- log(bmass[dat$species])

# setup trait matrix
# trait <- cbind(shape, length=tonguelengths, width=tonguewidths,diet,BodMass)

rownames(dat) <- dat$species
dat <- dat[, -1]

# drop tips without data
tre2 <- drop.tip(tre, which(!tre$tip %in% rownames(dat)))

# reorder trait matrix based on tip labels
dat <- dat[tre2$tip, ]

plm_tshape <- phylolm(aspect ~ fish + logmass, data = dat, phy = tre2, 
                 model = "lambda",
                 boot = 0, full.matrix = TRUE)
summary(plm_tshape) # diet 0.28 +/- 0.13 P = 0.046
                    # mass 0.12 +/- 0.09 P = 0.19

plm_tlength <- phylolm(length.mm ~ fish + logmass, data = dat,
	phy = tre2, model = "lambda", boot = 0, full.matrix = TRUE, lower.bound=1e-10)
summary(plm_tlength) # diet 0.22 +/- 0.52, P = 0.68
					 # mass 2.60 +/- 0.38 P < 0.001

plm_twidth <- phylolm(width.mm ~ fish + logmass, data = dat,
	phy = tre2, model = "lambda", boot = 0, full.matrix = TRUE, lower.bound=1e-10)
summary(plm_twidth) # diet -0.12 +/- 0.26, P = 0.64
					# mass 1.41 +/- 0.14 P < 0.011

ggplot(data = dat, aes(x=factor(fish), y=aspect, color=factor(fish))) + geom_boxplot() + theme_classic() + scale_color_manual(values=c('0'='gray', '1'='blue')) + labs(x="Plunge-diving", y="Tongue aspect ratio") + theme(legend.position="none")
ggsave("figs/tongue_boxplot.pdf", width=3.25, height=3)
```

#-------------------------------------------------------------------------------
# Correlate rates of gene evolution with tongue shape, etc.
#-------------------------------------------------------------------------------

```r
# devtools::install_github("nclark-lab/RERconverge")

# Time tree
timetree <- read.tree("data/paml/kingTree_abbrev.phy")

# Gene trees (~7 minutes for 10945 trees)
# trees_all <- readTrees("output/M0_genetrees.txt")
# saveRDS(trees_all, file="output/M0_genetrees_RER.rds")
trees_all <- readRDS(file="output/M0_genetrees_RER.rds")

trees_taste <- readTrees("output/taste_M0_genetrees.txt", masterTree=timetree)

# Phenotypic data
trait <- read.csv("output/tongue.csv", row=1)
rownames(trait) <- plyr::revalue(rownames(trait), c("Halcyon_senegaloides" = "hal_Sng"))
rownames(trait) <- sppAbbrev(rownames(trait))

# Overlapping species
# idx <- intersect(rownames(trait), trees_all[[1]][[1]]$tip.label)
idx <- intersect(rownames(trait), trees_taste[[1]][[1]]$tip.label)
length(idx)  # N = 23

# Phylogenetic residuals for tongue width
td <- treeplyr::make.treedata(timetree, trait)
# contMap(td$phy, log(td[["width"]]))
# contMap(td$phy, log(td[["area"]]))
# plot(area ~ BodMass, data=trait, log='x')

pres <- phyl.resid(td$phy, x=log(td[["BodMass"]]), Y=log(td[["width"]]))
pres <- phyl.resid(td$phy, x=log(td[["BodMass"]]), Y=log(td[["area"]]))
pres <- phyl.resid(td$phy, x=log(td[["BodMass"]]), Y=log(td[["length"]]))


coords_hyoid <- readRDS("~/Dropbox/Projects/king_tongues/output/coords_hyoid.rds")
curvesH <- as.matrix(read.csv("~/Dropbox/Projects/king_tongues/data/curves2.csv", header=F))
gpa_hyoid <- geomorph::gpagen(A = coords_hyoid, curves = curvesH, ProcD=FALSE, print.progress = T)
dat <- geomorph::two.d.array(gpa_hyoid$coords)
Species <- readxl::read_xlsx("~/Dropbox/Projects/king_tongues/data/kingfisher_CT_scans_050817.xlsx")
spp <- Species$spp_andersen2017[match(rownames(dat), Species$ID)]
spp <- gsub(" ", "_", spp)
rownames(dat) <- sppAbbrev(spp)
pca <- prcomp(dat)
# plot(pca)
# lines(vegan::bstick(pca))

# Setup
cnee_trees <- readRDS(file="output/baseml_free_ASHCE_trees_RER.rds")
gene_trees <- readRDS(file="output/M0_genetrees_RER.rds")


rer_cnee <- getAllResiduals(cnee_trees, plot=TRUE, transform="sqrt", scale=TRUE, weighted=TRUE, useSpecies=idx)  # fast
rer_gene <- getAllResiduals(gene_trees, plot=TRUE, transform="sqrt", scale=TRUE, weighted=TRUE, useSpecies=idx)  # fast

#-------------------------------------------------------------------------------
# Permulations for tongue shape
#-------------------------------------------------------------------------------

# Genes
res_tongue_gene <- correlateWithContinuousPhenotype(RERmat=rer_gene, charP=charpaths)
permCC_tongue_gene <- getPermsContinuous(numperms=500, traitvec=td[["shape"]], RERmat=rer_gene, trees=gene_trees, mastertree=keep.tip(gene_trees$masterTree, td$phy$tip.label), calculateenrich=FALSE)
permpvalCC_tongue_gene <- permpvalcor_fast(res_tongue_gene, permCC_tongue_gene)
save(res_tongue_gene, permCC_tongue_gene, permpvalCC_tongue_gene, file="output/perm_rer_tongue_gene.rda")

# CNEEs
res_tongue_cnee <- correlateWithContinuousPhenotype(RERmat=rer_cnee, charP=charpaths)
permCC_tongue_cnee <- getPermsContinuous(numperms=500, traitvec=td[["shape"]], RERmat=rer_cnee, trees=cnee_trees, mastertree=keep.tip(cnee_trees$masterTree, td$phy$tip.label), calculateenrich=FALSE)
permpvalCC_tongue_cnee <- permpvalcor_fast(res_tongue_cnee, permCC_tongue_cnee)
save(res_tongue_cnee, permCC_tongue_cnee, permpvalCC_tongue_cnee, file="output/perm_rer_tongue_cnee.rda")

# Stats/output
picks <- names(which(p.adjust(pvals_genes_tongue, method="fdr") < .15))
picks <- gsub(":", "_", picks)
clipr::write_clip(unique(na.omit(annot$finalgene[match(picks, annot$file)])))

# Just taste receptors
# RERconverge - tongue shape (aspect ratio)
rer_tongue_taste <- getAllResiduals(trees_taste, plot=TRUE, transform="sqrt", scale=TRUE, weighted=TRUE, useSpecies=idx)  # fast
# rer <- getAllResiduals(trees_all, plot=TRUE, transform="sqrt", scale=TRUE, weighted=TRUE, useSpecies=idx)  # fast
charpaths <- char2Paths(tip.vals=td[["shape"]], treesObj=trees_taste, metric="diff")
# charpaths <- char2Paths(tip.vals = td[["shape"]], treesObj = trees_all, altMasterTree = mtree)

par(mfrow=c(1,2))
hist(res_tongue_gene$P, main="Genes tongue")
hist(res_tongue_cnee$P, main="CNEs tongue")


res_tongue_taste <- correlateWithContinuousPhenotype(RERmat = rer_tongue_taste, charP = charpaths)

permCC_tongue_taste <- getPermsContinuous(numperms=1000, traitvec=td[["shape"]], RERmat=rer_tongue_taste, trees=trees_taste, mastertree=keep.tip(trees_taste$masterTree, td$phy$tip.label), calculateenrich=FALSE)
permpvalCC_tongue_taste <- permpvalcor(res_tongue_taste, permCC_tongue_taste)

p.adjust(permpvalCC_tongue_taste[c('Calhm1','Trpm5')], method="fdr")
p.adjust(permpvalCC_tongue_taste[c('Tas1r1','Tas1r3','Scnn1a','Scnn1b','Scnn1g')], method="fdr")

padj <- p.adjust(permpvalCC_tongue_taste, method="fdr")



# Plot
pdf("figs/taste_rer_sig.pdf" ,width=6.5, height=3.25)
par(oma=c(0,2,0,0), mfrow=c(1,2), mar=c(3,1,1,1), mgp=c(1.5,.5,0),ps=8, mex=.75)
x <- charpaths
y <- rer_tongue_taste['Scnn1g', ]
plot(x, y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Tongue shape change", 
     ylab="Molecular Rate", main="",
     pch=19, cex=1, xlim=c(-.25,.3))# ,col=ifelse(names(y) %in% plungeFg, "blue", "gray"))
text(x,y, labels=names(y), pos=4)
abline(lm(y ~ x), col='black', lwd=2)
x=charpaths
y=rer['Trpm5', ]
plot(x, y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Tongue shape change", 
     ylab="Molecular Rate", main="",
     pch=19, cex=1, xlim=c(-.25,.3)) #     col = ifelse(names(y) %in% plungeFg, "blue", "gray"))
text(x,y, labels=names(y), pos=4)
abline(lm(y ~ x), col='black', lwd=2)
# mtext(side=1, line=.5, outer=T, text="Tongue shape change")
mtext(side=2, line=1, outer=T, text="Molecular rate")
dev.off()

#-------------------------------------------------------------------------------
# hyoid shape
#-------------------------------------------------------------------------------
charpaths <- char2Paths(tip.vals = pca$x[, 1], treesObj = trees_taste, altMasterTree = mtree)
charpaths <- char2Paths(tip.vals = pca$x[, 2], treesObj = trees_taste, altMasterTree = mtree)
res <- correlateWithContinuousPhenotype(RERmat = rer, charP = charpaths)

res <- getAllCor(rer, charpaths, method="p", min.pos=0, winsorizeRER=3)
# res

perms <- getPermsContinuous(numperms=100, traitvec=td[["shape"]], RERmat=rer, trees=trees_all, mastertree=mtree, calculateenrich=FALSE)
# perms
corpermpvals=permpvalcor(res, perms)
plot(corpermpvals)
hist(res$P)
hist(corpermpvals, add=T, col='red')

# Figure output
jpg("figs/rer_tas1r3_tongue.jpg", width=5, height=5)
x=charpaths
y=rer['Tas1r3', ]
plot(x, y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Weight Change", 
     ylab="Evolutionary Rate", main="Gene TTN Pearson Correlation",
     pch=19, cex=1, xlim=c(-.25,.3))
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
dev.off()

# x=charpaths
# y=rer['Tas2r7', ]
plot(x, y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Tongue shape change", 
     ylab="Molecular Rate", main="Gene TTN Pearson Correlation",
     pch=19, cex=1)
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
for (i in 2:nrow(rer)) {
	x=charpaths
	y=rer[i, ]
	abline(lm(y~x), col='red',lwd=3)
}

# RERconverge - plunge-diving behavior
fg = c(rownames(trait)[which(trait[,'diet']==1)], 'megAlc','ceyMel')
mtree2 = foreground2Tree(fg, trees_all, clade="terminal")
plot(mtree2)
tp <- tree2Paths(mtree2, trees_all)
rer2 <- getAllResiduals(trees_all, plot=TRUE, transform="sqrt", scale=TRUE, weighted=TRUE)
res2 <- correlateWithBinaryPhenotype(rer2, tp, min.sp=3)
table(res2$p.adj < 0.05) # nothing significant
# Figure output
p <- plotRers(rer2, "Tas1r3", phenv=tp)
p + scale_color_manual(values = c("gray", "blue"))
ggsave("figs/rer_tas1r3_dive.jpg", width=6.5, height=5)

# Beak width and BMP4
mtree$node.label <- NULL
tr <- readTrees("data/BMP4_M0_genetree.txt")
beak <- read.csv("/Users/chadeliason/Dropbox (The Field Museum)/Projects/Coraciiform morphometrics/data/Linear-dataset/morph-data.csv")
beak <- beak %>% group_by(Species) %>% summarize_if(is.numeric, mean)
beak$spp <- sppAbbrev(beak$Species)
trait <- beak$width.nares.1[match(tr$trees[[1]]$tip, beak$spp)]
names(trait) <- tr$trees[[1]]$tip
rer <- getAllResiduals(tr, transform="sqrt", plot=TRUE, scale=F, weighted=TRUE)
charpaths <- char2Paths(tip.vals=trait, treesObj=tr)
res <- correlateWithContinuousPhenotype(rer, charpaths, winsorizeRER = 3, winsorizetrait = 3)
res[1, ]

# Brain shape and MAPT...

```

## Permulation tests of RERconverge results

```r
setwd("~/Dropbox/Projects/king_genome")

library(stringr)
library(RERconverge)

source("R/genome-sensory-ms-functions.R")

# Load trees from PAML (gene_trees) and BASEML (cnee_trees)
gene_trees <- readRDS(file="output/M0_genetrees_RER.rds")
cnee_trees <- readRDS(file="output/baseml_free_ASHCE_trees_RER.rds")

# Timetree
mtree <- read.tree("data/paml/kingTree_abbrev.phy")

# Define sister clades
sisters_plunge <- list(clade1=c('corVin','corCri'), clade2=c('clade1','corLeu'),
					   clade3=c('ceyArg','ceyCya'), clade4=c('alcSem','alcQua'),
					   clade5=c('clade4','alcAt'), clade6=c('megTor','megAlc'),
					   clade7=c('clade5','megMax'), clade8=c('chlInd','chlAma'),
					   clade9=c('clade8','cerRud'), clade10=c('clade7','clade8'))

sisters_island <- list(clade1=c('ceyArg','ceyCya'), clade2=c('clade1','ceyMar'),
					   clade3=c('clade2', 'ceyMel'), clade4=c('actLin','actHom'))

# Plotting trees to check

plungeFgTree <- foreground2TreeClades(plungeFg, sisters_plunge, gene_trees, plotTree=F)

pdf("figs/plunge_tree.pdf", width=4, height=5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(plungeFgTree, hlspecies = which(plungeFgTree$edge.length==1),
	hlcols="blue", main="Kingfishers trait tree")
dev.off()

islandFgTree <- foreground2TreeClades(islandFg, sisters_island, gene_trees, plotTree=F)

pdf("figs/island_tree.pdf", width=4, height=5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(islandFgTree, hlspecies=which(islandFgTree$edge.length==1),
	hlcols="red", main="Kingfishers trait tree")
dev.off()

# Calculating paths from the foreground tree
pathvec_plunge <- foreground2Paths(foreground = plungeFg, treesObj = gene_trees, clade="all")
pathvec_island <- foreground2Paths(foreground = islandFg, treesObj = gene_trees, clade="all")

# Calculate RERs
kingRERw_gene <- getAllResiduals(gene_trees, transform="sqrt", weighted=T, scale=T)
kingRERw_cnee <- getAllResiduals(cnee_trees, transform="sqrt", weighted=T, scale=T)

# Calculate correlations
res_plunge_gene <- correlateWithBinaryPhenotype(kingRERw_gene, pathvec_plunge, min.sp=2, min.pos=2, weighted="auto")
res_island_gene <- correlateWithBinaryPhenotype(kingRERw_gene, pathvec_island, min.sp=2, min.pos=2, weighted="auto")
res_plunge_cnee <- correlateWithBinaryPhenotype(kingRERw_cnee, pathvec_plunge, min.sp=2, min.pos=2, weighted="auto")
res_island_cnee <- correlateWithBinaryPhenotype(kingRERw_cnee, pathvec_island, min.sp=2, min.pos=2, weighted="auto")

# Define the root species
root_sp = "momMom"

# Extract master trees from genes/CNEEs
masterTree_gene = gene_trees$masterTree
masterTree_cnee = cnee_trees$masterTree

#-------------------------------------------------------------------------------
# Perform 500 permulations (UNCOMMENT to run- takes ~2h for CNEE runs)
#-------------------------------------------------------------------------------

# permCC_plunge_gene <- getPermsBinary(500, plungeFg, sisters_plunge, root_sp, kingRERw_gene, gene_trees, masterTree_gene, permmode="cc", ncores=4)
# permCC_island_gene <- getPermsBinary(500, islandFg, sisters_list=NULL, root_sp, kingRERw_gene, gene_trees, masterTree_gene, permmode="cc", ncores=4)
# permCC_island_cnee <- getPermsBinary(500, islandFg, sisters_list=sisters_island, root_sp, kingRERw_cnee, cnee_trees, masterTree_cnee, permmode="cc", ncores=4)
# permCC_plunge_cnee <- getPermsBinary(500, plungeFg, sisters_list=sisters_plunge, root_sp, kingRERw_cnee, cnee_trees, masterTree_cnee, permmode="cc", ncores=4)

#-------------------------------------------------------------------------------
# Calculate P values
#-------------------------------------------------------------------------------
permpvalCC_plunge_gene <- permpvalcor_fast(res_plunge_gene, permCC_plunge_gene)
permpvalCC_island_gene <- permpvalcor_fast(res_island_gene, permCC_island_gene)
permpvalCC_island_cnee <- permpvalcor_fast(res_island_cnee, permCC_island_cnee)
permpvalCC_plunge_cnee <- permpvalcor_fast(res_plunge_cnee, permCC_plunge_cnee)

#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------
save(res_island_gene, permCC_island_gene, permpvalCC_island_gene, file="output/perm_rer_island_gene.rda")
save(res_plunge_gene, permCC_plunge_gene, permpvalCC_plunge_gene, file="output/perm_rer_plunge_gene.rda")
save(res_island_cnee, permCC_island_cnee, permpvalCC_island_cnee, file="output/perm_rer_island_cnee.rda")
save(res_plunge_cnee, permCC_plunge_cnee, permpvalCC_plunge_cnee, file="output/perm_rer_plunge_cnee.rda")

#-------------------------------------------------------------------------------
# Generate figure
#-------------------------------------------------------------------------------

pdf("figs/rer_histo_plunge.pdf", width=3.5, height=3)
par(oma = c(0,0,0,0), mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.5,.5,0), mex=.75, ps=8)
plot(density(res_plunge_cnee$Rho), main="", xlab="Correlation")
lines(density(res_plunge_gene$Rho), col=2)
legend("topleft", bty="n", legend=c("genes", "CNEEs"), col=2:1, lwd=1)
abline(v=0, lty=2)
dev.off()

pdf("figs/rer_histo_island.pdf", width=3.5, height=3)
par(oma = c(0,0,0,0), mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.5,.5,0), mex=.75, ps=8)
plot(density(res_island_cnee$Rho), main="", xlab="Correlation")
lines(density(res_island_gene$Rho), col=2)
legend("topleft", bty="n", legend=c("genes", "CNEEs"), col=2:1, lwd=1)
abline(v=0, lty=2)
dev.off()

#-------------------------------------------------------------------------------
## Stats
#-------------------------------------------------------------------------------
# load(file="output/perm_rer_island_gene.rda")
load(file="output/perm_rer_plunge_gene.rda")
# load(file="output/perm_rer_island_cnee.rda")
load(file="output/perm_rer_plunge_cnee.rda")

padj_plunge_cnee <- p.adjust(permpvalCC_plunge_cnee, method="fdr")
# padj_island_cnee <- p.adjust(permpvalCC_island_cnee, method="fdr")
padj_plunge_gene <- p.adjust(permpvalCC_plunge_gene, method="fdr")
# padj_island_gene <- p.adjust(permpvalCC_island_gene, method="fdr")
names(padj_plunge_gene) <- gsub(":", "_", names(padj_plunge_gene))
names(padj_plunge_cnee) <- gsub(":", "_", names(padj_plunge_cnee))
# names(padj_island_gene) <- gsub(":", "_", names(padj_island_gene))
# names(padj_island_cnee) <- gsub(":", "_", names(padj_island_cnee))

table(padj_plunge_cnee < .05)  # 37
# table(padj_island_cnee < .05)  # 44
table(padj_plunge_gene < .05)  # 21
# table(padj_island_gene < .05)  # 43

# Find accelerations/deccelerations
table(res_plunge_gene[which(padj_plunge_gene < 0.05), "Rho"] >= 0)
# table(res_island_gene[which(padj_island_gene < 0.05), "Rho"] >= 0)
table(res_plunge_cnee[which(padj_plunge_cnee < 0.05), "Rho"] >= 0)
# table(res_island_cnee[which(padj_island_cnee < 0.05), "Rho"] >= 0)

id1=which(padj_plunge_cnee<.05)
id2=which(padj_island_cnee<.05)
intersect(id1,id2) # none in common

id3=which(padj_plunge_gene<.05)
id4=which(padj_island_gene<.05)
intersect(id3,id4) # none in common

annot <- read.csv("output/pos_selection_results_repeatmasked_20221209_absrel_relax_pcoc.csv", row=1)
annot$fileshort <- gsub("_R[0-9]$", "", annot$file)

# Significantly accelerated genes
picks <- names(which(padj_plunge_gene < .15 & res_plunge_gene$Rho > 0))
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

# Significantly decelerated genes
picks <- names(which(padj_plunge_gene < .15 & res_plunge_gene$Rho <= 0))
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

# Significant plunge-diving genes
picks <- names(which(padj_plunge_gene < .15))
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

# Significant island-dwelling genes
picks <- names(which(padj_island_gene < .15))
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

# Output to find closest genes
picks <- names(c(which(padj_plunge_cnee < .15), which(padj_island_cnee < .15)))
out <- do.call(rbind, strsplit(picks, "_"))[, 1:3]
head(out)
write.table(out, row.names=F, col.names=F, quote=F, file="output/rer_cnee_picks.bed", sep="\t")
# use this to find closest in phoebe...
```


```sh
cd ~/uce-alcedinidae

# Get CNEEs closest to annotation genes (GeMoMa)
wc -l annotations/gemoma/round2/filtered_predictions.mRNA.bed  # 16540

# Filter out CNEEs overlapping transcripts in Todiramphus chloris
bedtools intersect -v -a cnee_redo/todChl_cnee.bed -b annotations/gemoma/round2/filtered_predictions.mRNA.bed > cnee_redo/todChl_ASHCE_nonoverlap.bed
wc -l cnee_redo/todChl_ASHCE_nonoverlap.bed # 194298 non-overlapping CNEEs
awk '{if(($3-$2)>=50) print}' cnee_redo/todChl_ASHCE_nonoverlap.bed > cnee_redo/todChl_ASHCE_nonoverlap_50bp.bed
wc -l cnee_redo/todChl_ASHCE_nonoverlap_50bp.bed  # 45613

# Find CNEEs closest to each gene, ignore overlap (-io)
cat cnee_redo/rer_cnee_picks.bed

bedtools sort -i cnee_redo/rer_cnee_picks.bed |\
bedtools closest -io -d -a - -b annotations/gemoma/round2/filtered_predictions.mRNA.bed |\
bedtools sort > cnee_redo/todChl_ASHCE_nonoverlap_closest_RER.bed
```

Now more plots of RER results

```r
bed <- read.delim("output/todChl_ASHCE_nonoverlap_closest_RER.bed", head=F)
bed$V7 <- gsub(":", "_", bed$V7)
bed$id <- paste0(bed$V1, "_", bed$V2, "_", bed$V3)

# Plunge-diving CNEEs (for STRING)
idx <- names(which(padj_plunge_cnee < .15))
picks <- bed$V7[match(gsub("_trimmed", "", idx), bed$id)]
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

idx <- names(which(padj_plunge_cnee < .15 & res_plunge_cnee$Rho > 0))
picks <- bed$V7[match(gsub("_trimmed", "", idx), bed$id)]
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

idx <- names(which(padj_plunge_cnee < .15 & res_plunge_cnee$Rho < 0))
picks <- bed$V7[match(gsub("_trimmed", "", idx), bed$id)]
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

# Island-dwelling CNEEs (for STRING)
idx <- names(which(padj_island_cnee < .15))
picks <- bed$V7[match(gsub("_trimmed", "", idx), bed$id)]
cat(sort(unique(na.omit(annot$finalgene[match(picks, annot$file)]))), sep="\n")

# Histograms
hist(res_plunge_gene$P)
hist(res_plunge_cnee$P)
hist(res_island_gene$P)
hist(res_island_cnee$P)

# Plots for figure (plunge-diving genes)
idx <- which(padj_plunge_gene < .15 & res_plunge_gene$Rho > 0)
idx <- which(padj_plunge_gene < .15 & res_plunge_gene$Rho < 0)
res_plunge_gene[idx, ]
hi="aquChr_transcript:ENSACCT00020012913_R1"
lo="aquChr_transcript:ENSACCT00020009335_R0"
pdf("figs/rer_plunge_accel_gene_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(gene_trees[[1]][[hi]], hlspecies=which(plungeFgTree$edge.length==1), hlcols="blue", main="", outgroup="momMom")
dev.off()
pdf("figs/rer_plunge_decel_gene_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(gene_trees[[1]][[lo]], hlspecies=which(plungeFgTree$edge.length==1), hlcols="blue", main="", outgroup="momMom")
dev.off()

annot$finalgene[match(hi, annot$file)] # CSF2RA gene
annot$finalgene[match(gsub(":", "_", lo), annot$file)] # HELZ gene


# Plots for figure (plunge-diving CNEEs)
idx <- which(padj_plunge_cnee < .15 & res_plunge_cnee$Rho < 0)
dplyr::arrange(res_plunge_cnee[idx, ], Rho)
lo="129_11935471_11935626_trimmed"
idx <- which(padj_plunge_cnee < .15 & res_plunge_cnee$Rho > 0)
dplyr::arrange(res_plunge_cnee[idx, ], Rho)
hi="193_4562700_4562763_trimmed"
pdf("figs/rer_plunge_accel_cnee_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(cnee_trees[[1]][[hi]], hlspecies=which(plungeFgTree$edge.length==1), hlcols="blue", main="", outgroup="momMom")
dev.off()
pdf("figs/rer_plunge_decel_cnee_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(cnee_trees[[1]][[lo]], hlspecies=which(plungeFgTree$edge.length==1), hlcols="blue", main="", outgroup="momMom")
dev.off()

idx <- gsub(":", "_", bed$V7[which(bed$id=="129_11935471_11935626")]) # accel plunge-diving
annot$finalgene[which(annot$file==idx)]  # SCOC
idx <- bed$V7[which(bed$id=="193_4562700_4562763")] # decel plunge-diving
annot$finalgene[which(annot$file==idx)]  # ADGRB1


# GO plots

# idx <- which(padj_plunge_gene < .15)
# picks <- sort(unique(na.omit(annot$finalgene[match(names(idx), annot$file)])))
# g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
# pdf("figs/emap_plunge_gene.pdf", width=7, height=6)
# emapplot(gost_to_enrich(g), showCategory=30, layout="kk")
# dev.off()


# idx <- which(padj_island_gene < .15)
# picks <- sort(unique(na.omit(annot$finalgene[match(names(idx), annot$file)])))
# g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
# pdf("figs/emap_island_gene.pdf", width=7, height=6)
# emapplot(gost_to_enrich(g), showCategory=10, layout="kk")
# dev.off()

my_scale_size <- scale_size(limits=c(1, 7), range=c(4,10))

idx <- which(padj_plunge_gene < .15 & res_plunge_gene$Rho > 0)
picks <- sort(unique(na.omit(annot$finalgene[match(names(idx), annot$file)])))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_plunge_gene_pos.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=2, high=2) + my_scale_size
dev.off()
pdf("figs/heatplot_plunge_gene_pos.pdf", width=7, height=6)
heatplot(gost_to_enrich(g), showCategory=10)
dev.off()

idx <- which(padj_plunge_gene < .15 & res_plunge_gene$Rho < 0)
picks <- sort(unique(na.omit(annot$finalgene[match(names(idx), annot$file)])))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_plunge_gene_neg.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=2, high=2) + my_scale_size
dev.off()
pdf("figs/heatplot_plunge_gene_neg.pdf", width=7, height=6)
heatplot(gost_to_enrich(g), showCategory=10)
dev.off()

# idx <- which(padj_plunge_gene < .15 & res_plunge_gene$Rho > 0)
# picks <- sort(unique(na.omit(annot$finalgene[match(names(idx), annot$file)])))
# g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
# pdf("figs/emap_plunge_gene_pos.pdf", width=7, height=6)
# emapplot(gost_to_enrich(g), showCategory=10, layout="kk")
# dev.off()

# heatplot(gost_to_enrich(g), showCategory=30, symbol="dot")



idx <- which(padj_island_gene < .15 & res_island_gene$Rho < 0)
picks <- sort(unique(na.omit(annot$finalgene[match(gsub(":", "_", names(idx)), annot$file)])))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_island_gene_neg.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=2, high=2) + my_scale_size
dev.off()

idx <- which(padj_island_gene < .15 & res_island_gene$Rho > 0)
picks <- sort(unique(na.omit(annot$finalgene[match(gsub(":", "_", names(idx)), annot$file)])))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_island_gene_pos.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=2, high=2) + my_scale_size
dev.off()

# CNEE GO plots
idx <- which(padj_plunge_cnee < .15 & res_plunge_cnee$Rho > 0)
picks <- bed$V7[match(gsub("_trimmed", "", names(idx)), bed$id)]
picks <- gsub("_R[0-9]$", "", picks)
picks <- sort(na.omit(annot$finalgene[match(picks, annot$fileshort)]))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_plunge_cnee_pos.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=1, high=1) + my_scale_size
dev.off()
pdf("figs/heatplot_plunge_cnee_pos.pdf", width=7, height=6)
heatplot(gost_to_enrich(g), showCategory=10)
dev.off()

idx <- which(padj_plunge_cnee < .15 & res_plunge_cnee$Rho < 0)
picks <- bed$V7[match(gsub("_trimmed", "", names(idx)), bed$id)]
picks <- gsub("_R[0-9]$", "", picks)
picks <- sort(na.omit(annot$finalgene[match(picks, annot$fileshort)]))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_plunge_cnee_neg.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=1, high=1) + my_scale_size
dev.off()
pdf("figs/heatplot_plunge_cnee_neg.pdf", width=7, height=6)
heatplot(gost_to_enrich(g), showCategory=10)
dev.off()


idx <- which(padj_island_cnee < .15 & res_island_cnee$Rho > 0)
picks <- bed$V7[match(gsub("_trimmed", "", names(idx)), bed$id)]
picks <- gsub("_R[0-9]$", "", picks)
picks <- sort(na.omit(annot$finalgene[match(picks, annot$fileshort)]))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_island_cnee_pos.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=1, high=1) + my_scale_size
dev.off()

idx <- which(padj_island_cnee < .15 & res_island_cnee$Rho < 0)
picks <- bed$V7[match(gsub("_trimmed", "", names(idx)), bed$id)]
picks <- gsub("_R[0-9]$", "", picks)
picks <- sort(na.omit(annot$finalgene[match(picks, annot$fileshort)]))
g = gost(query=picks, organism="hsapiens", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=F, user_threshold=0.15)
pdf("figs/emap_island_cnee_neg.pdf", width=7, height=6)
emapplot(gost_to_enrich(g), showCategory=10, layout="kk") + scale_color_gradient(low=1, high=1) + my_scale_size
dev.off()

# barplot(gost_to_enrich(g), showCategory=10) 
# heatplot(gost_to_enrich(g), showCategory=30)


# Plots for figure (island genes)
idx <- which(padj_island_gene < .15)
dplyr::arrange(res_island_gene[idx, ], Rho)
lo="aquChr_transcript_ENSACCT00020013663_R0"
idx <- which(padj_island_gene < .15 & res_island_gene$Rho > 0)
dplyr::arrange(res_island_gene[idx, ], Rho)
hi="aquChr_transcript_ENSACCT00020003476_R0"
pdf("figs/rer_island_accel_gene_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(gene_trees[[1]][[hi]], hlspecies=which(islandFgTree$edge.length==1), hlcols="red", main="", outgroup="momMom")
dev.off()
pdf("figs/rer_island_decel_gene_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(gene_trees[[1]][[lo]], hlspecies=which(islandFgTree$edge.length==1), hlcols="red", main="", outgroup="momMom")
dev.off()

annot$finalgene[match(hi, annot$file)] # TLN2 gene
annot$finalgene[match(lo, annot$file)] # PIMREG gene

# Plots for figure (island CNEEs)
idx <- which(padj_island_cnee < .15 & res_island_cnee$Rho < 0)
dplyr::arrange(res_island_cnee[idx, ], Rho)
lo="74_17372822_17372941_trimmed"  # 2nd lowest, bc first wasn't an annotated gene
idx <- which(padj_island_cnee < .15 & res_island_cnee$Rho > 0)
dplyr::arrange(res_island_cnee[idx, ], Rho)
hi="71_57179131_57179251_trimmed"
pdf("figs/rer_island_accel_cnee_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(cnee_trees[[1]][[hi]], hlspecies=which(islandFgTree$edge.length==1), hlcols="red", main="", outgroup="momMom")
dev.off()
pdf("figs/rer_island_decel_cnee_tree.pdf", width=4.3, height=6.5)
par(mar=c(0,0,0,0))
plotTreeHighlightBranches(cnee_trees[[1]][[lo]], hlspecies=which(islandFgTree$edge.length==1), hlcols="red", main="", outgroup="momMom")
dev.off()

idx = bed$V7[which(bed$id==gsub("_trimmed", "", hi))] # accel plunge-diving
annot$finalgene[which(annot$file==idx)]  # LEG1 gene
idx=bed$V7[which(bed$id==gsub("_trimmed", "", lo2))] # accel plunge-diving
annot$finalgene[which(annot$file==idx)]  # CLDN10 gene

# Output for STRING, gprofiler, etc.
picks <- names(which(p.adjust(permpvalCC_island, method="fdr") < .15))
picks <- rownames(res_island)[which(res_island$P<.05)]
genes <- unique(annot$finalgene[match(picks, annot$file)])
genes <- genes[!is.na(genes)]
clipr::write_clip(genes)

# par(ask=T)
# for (i in names(padj[padj<.05])) {
# 	tr=gene_trees$trees[[i]]
# 	plot(tr, tip.col=ifelse(tr$tip %in% plungeFg, "blue", "black"), main=i)
# }

# [ ] Now I need to add this to the upset plot...and do an enrichment analysis

pathvec_tongueshape <- char2Paths(tip.vals=td[["shape"]], treesObj=gene_trees)
rer2 <- getAllResiduals(gene_trees, transform="sqrt", weighted=T, scale=T, useSpecies=names(td[["shape"]]))
res_tongue <- correlateWithContinuousPhenotype(RERmat=rer2, charP=pathvec_tongueshape, min.sp=2)
plot(res_tongue$p.adj)

#-------------------------------------------------------------------------------
# Brain shape
#-------------------------------------------------------------------------------

brainpc1 <- setNames(brain$brainPC1, brain$abbrev)
# brainpc2 <- setNames(brain$brainPC2, brain$abbrev)
# brainpc3 <- setNames(brain$brainPC3, brain$abbrev)
# brainpc4 <- setNames(brain$brainPC4, brain$abbrev)
# brainsize <- setNames(brain$ln.brainsize, brain$abbrev)

pathvec_brainpc1 <- char2Paths(tip.vals=brainpc1, treesObj=gene_trees)
# pathvec_brainpc2 <- char2Paths(tip.vals=brainpc2, treesObj=gene_trees)
# pathvec_brainpc3 <- char2Paths(tip.vals=brainpc3, treesObj=gene_trees)
# pathvec_brainpc4 <- char2Paths(tip.vals=brainpc4, treesObj=gene_trees)
# pathvec_brain_size <- char2Paths(tip.vals=brainsize, treesObj=gene_trees)

# Genes
rer3 <- getAllResiduals(gene_trees, transform="sqrt", weighted=T, scale=T, useSpecies=names(brainpc1))
res_brainPC1_gene <- correlateWithContinuousPhenotype(RERmat=rer3, charP=pathvec_brainpc1, min.sp=2)
# res_brain2 <- correlateWithContinuousPhenotype(RERmat=rer3, charP=pathvec_brainpc2, min.sp=2)
# res_brain3 <- correlateWithContinuousPhenotype(RERmat=rer3, charP=pathvec_brainpc3, min.sp=2)
# res_brain4 <- correlateWithContinuousPhenotype(RERmat=rer3, charP=pathvec_brainpc4, min.sp=2)
# res_brain_size <- correlateWithContinuousPhenotype(RERmat=rer3, charP=pathvec_brain_size, min.sp=2)
masterTree_gene = gene_trees$masterTree
mtree2 = drop.tip(masterTree_gene, setdiff(masterTree_gene$tip.label, names(brainpc1)))
trait <- setNames(brain$brainPC1, brain$abbrev)
permCC_brainPC1_gene <- getPermsContinuous(numperms = 500, traitvec = trait, RERmat = rer3, trees = gene_trees, mastertree = mtree2, calculateenrich = FALSE)
permpvalCC_brainPC1_gene <- permpvalcor_fast(res_brainPC1_gene, permCC_brainPC1_gene)
save(res_brainPC1_gene, permCC_brainPC1_gene, permpvalCC_brainPC1_gene,
	 file="output/perm_rer_brainPC1_gene.rda")

# CNEs
rer_brainPC1_cnee <- getAllResiduals(cnee_trees, transform="sqrt", weighted=T, scale=T, useSpecies=names(brainpc1))
res_brainPC1_cnee <- correlateWithContinuousPhenotype(RERmat=rer_brainPC1_cnee, charP=pathvec_brainpc1, min.sp=2)
permCC_brainPC1_cnees <- getPermsContinuous(numperms=500, traitvec=trait, RERmat=rer_brainPC1_cnee, trees=cnee_trees, mastertree=cnee_trees$masterTree, calculateenrich=FALSE)
permpvalCC_brainPC1_cnee <- permpvalcor_fast(res_brainPC1_cnee, permCC_brainPC1_cnees)
save(res_brainPC1_cnee, rer_brainPC1_cnee, permCC_brainPC1_cnees, permpvalCC_brainPC1_cnee,
	 file="output/perm_rer_brainPC1_cnee.rda")

# Plots
par(mfrow=c(2,2))
hist(res_brainPC1_gene$P); abline(v=.05, lty=2)
hist(res_brainPC1_cnee$P); abline(v=.05, lty=2)
hist(res_gene_tongue$P); abline(v=.05, lty=2)
hist(res_cnee_tongue$P); abline(v=.05, lty=2)

# Analysis
table(p.adjust(permpvalCC_brainPC1_cnee, method="fdr") < 0.05)  # 64
table(p.adjust(permpvalCC_brainPC1_gene, method="fdr") < 0.05)  # 17

sum(p.adjust(permpvalCC_tongue_cnee, method="fdr") < 0.05)  # 75
sum(p.adjust(permpvalCC_tongue_gene, method="fdr") < 0.05)  # 39

nrow(res_brainPC1_cnee)
nrow(res_brainPC1_gene)

fisher.test(rbind(c(21, 10945 - 21), c(37, 39618 - 37)))  # plunge-diving
fisher.test(rbind(c(39, 10945 - 39), c(75, 39618 - 75)))  # tongue shape
fisher.test(rbind(c(17, 10945 - 17), c(64, 39618 - 64)))  # brain PC1

# compare genes brain shape vs behavior
fisher.test(rbind(c(17, 10945 - 17), c(21, 10945 - 21))) # or=0.81, p=0.63
# compare genes tongue shape vs behavior
fisher.test(rbind(c(39, 10945 - 39), c(21, 10945 - 21))) # or=1.86, p=0.027
# compare CNEs brain shape vs behavior
fisher.test(rbind(c(64, 39618 - 64), c(37, 39618 - 37))) # or=1.73, p=0.0093
# compare CNEs tongue shape vs behavior
fisher.test(rbind(c(75, 39618 - 75), c(37, 39618 - 37))) # or=2.03 P<.001

# fisher.test(rbind(c(21, 10945 - 21), c(37, 39618 - 37)))  # plunge-diving

# Output gene names for webgestalt, STRING protein analysis...
padj <- p.adjust(permpvalCC_brainPC1_gene, method="fdr")
picks <- sort(names(padj[padj<.05]))
clipr::write_clip(unique(na.omit(annot$finalgene[match(gsub(":", "_", picks), annot$file)])))
```

## Reshaping files, phastcons...

```sh
# script modified from taylor
cd ~/uce-alcedinidae
mkdir genomes_aligned

for sp in `ls genomes`

do
	if [[ $sp =~ "todChl" ]]; then continue; fi
	# echo $sp
	sed -r 's/^>[0-9]+/& '$sp'/g' genomes/${sp}/${sp}-to-todChl.rg.indelrealigner.consensus.fa > genomes_aligned/"$sp"_renamed.fa
done

cp genomes/todChl/todChl.scaffolds.full_mask.fa genomes_aligned/todChl_renamed.fa
sed -r -i 's/^>[0-9]+/& todChl/g' genomes_aligned/todChl_renamed.fa

cat genomes_aligned/*_renamed.fa > genomes_aligned/all.fasta


for loci in $(loci_list.txt)
do
	# loci=60
	# echo $loci
	# using line length of 80
	echo $loci | seqtk subseq -l 80 genomes_aligned/all.fasta - | sed -r 's/^>[0-9]+/>/g' > genomes_aligned/contig_${loci}.fasta
done

# get 4-fold degenerate sites

# get consensus (code from Sackton et al. 2019)
phastCons \
	--expected-length=45 \
	--target-coverage=0.3 \
	--estimate-rho ./mods/$SAMP.rho \
	--no-post-probs \
	--msa-format FASTA \
	--log ./logs/$SAMP.log \
	../chunks/$SAMP.ss ../../../neutMods/neut_ver${VERSION}_final.named.mod &> ./logs/$SAMP.out

```


## New figure with absrel results

```{r}
#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------

dat_paml <- read.csv("output/allres_merged.csv", row=1)

islandFg <- c("corMad","corVin","ceyFal","ceyWeb","ceyMel","ceyCya","ceyArg","ceyMar","ceyLep", "ceyCol","ceyNig","actHom","actLin","todAlb","todDio","todFar","todCin","todRuf")

plungeFg <- c('corLeu','corVin','corCri','ceyCya','ceyArg','alcSem','alcQua','alcAtt','megMax','megTor','megAlc','chlInd','chlAma','cerRud')

# Plot the convergent origins on the tree

mtree <- read.tree("data/paml/kingTree_abbrev.phy")

mtree_pruned <- drop.tip(mtree, "momMom")

pal <- RColorBrewer::brewer.pal(8, "Dark2")[c(8,4)]

ecols <- rep('lightgray', Nedge(mtree_pruned))
ecols[match(which(mtree_pruned$tip %in% 'corVin'), mtree_pruned$edge[, 2])] <- pal[2]
ecols[match(which(mtree_pruned$tip %in% c('corLeu','corCri')), mtree_pruned$edge[, 2])] <- pal[1]
ecols[match(which(mtree_pruned$tip %in% c('ceyCya','ceyArg')), mtree_pruned$edge[, 2])] <- pal[2]
ecols[match(which(mtree_pruned$tip %in% c('alcSem','alcQua','alcAtt')), mtree_pruned$edge[, 2])] <- pal[1]
ecols[match(which(mtree_pruned$tip %in% c('megMax','megTor','megAlc','chlAma','chlInd','cerRud')), mtree_pruned$edge[, 2])] <- pal[1]

nodes <- sapply(list(c('corLeu','corCri'), c('ceyCya','ceyArg'), c('alcSem','alcAtt'), c('cerRud','megMax')), getMRCA, phy=mtree_pruned)

pdf("figs/focal_tree.pdf", width=2.5, height=3.5)
plot(mtree_pruned, edge.col=ecols, edge.width=1, no.margin=TRUE, cex=.6, type='ph')
nodelabels(text=1:4, node=nodes, frame="circ", bg='white', cex=.65)
dev.off()


#------------------------------------------------------------------------------
# Create data set for enrichment tests (e.g., GO terms-by-branch)
#------------------------------------------------------------------------------

# Picking focal (plunge-diving) species-
picks <- c('corVin','ceyCya','ceyArg','corLeu','corCri', 'alcSem','alcQua','alcAtt','megMax','megTor','megAlc','chlAma','chlInd','cerRud')

tmp <- pblapply(picks, function(p) {
	pvals <- dat_paml[,grep(p, colnames(dat_paml))]
	genelist <- dat_paml$finalgene[which(pvals<.05)]
	genelist <- unique(genelist[!is.na(genelist)])
	g <- gprofiler2::gost(query = genelist, organism="ggallus", correction_method="fdr", evcodes=TRUE, sources="GO:BP", significant=FALSE)
	x=gost_to_enrich(g)
	n <- enrichplot:::update_n(x, showCategory=30)
	geneSets <- enrichplot:::extract_geneSets(x, n)
	d <- enrichplot:::list2df(geneSets)
	table(d$categoryID)
})
nms <- lapply(tmp, names)
unique_nms <- unique(unlist(nms))
dat_enrich <- as.data.frame(sapply(1:length(nms), function(x) {
	tmp[[x]][unique_nms]
}))
names(dat_enrich) <- picks
dat_enrich[is.na(dat_enrich)] <- 0
dat_enrich$term <- unique_nms
head(dat_enrich)

#-------------------------------------------------------------------------------
# PGLS analysis for number of positively selected genes
#-------------------------------------------------------------------------------

DietData <- readxl::read_xlsx("data/kingdata_v2_edited.xlsx")
DietData$spp <- gsub("Ceyx melanurus melanurus", "Ceyx melanurus", DietData$sppAndersen, fixed=TRUE)
DietData$spp <- gsub(" ", "_", DietData$spp)
DietData$spp_abbrev <- sppAbbrev(DietData$spp)

df3 = data.frame(numpos=apply(dat_paml[,50:80], 2, function(x) sum(x<.05,na.rm=T)))
df3$spp <- gsub("P_aBSREL_", "", rownames(df3))
rownames(df3) <- df3$spp
df3$lnmass <- log(DietData$weight_g[match(df3$spp, DietData$spp_abbrev)])
df3$dives <- ifelse(df3$spp %in% plungeFg, 1, 0)
df3$island <- ifelse(df3$spp %in% islandFg, 1, 0)
df3 <- subset(df3, spp!="momMom")
phy <- ape::drop.tip(phy_genomes, "momMom")
df3 <- df3[phy$tip, ]

lm1 <- phylostep(numpos ~ (dives + lnmass + island)^2, data=df3, phy=phy_genomes, model="lambda")
summary(lm1)

lm1_best <- phylolm(numpos ~ dives + lnmass * island, data=df3, phy=phy_genomes, model="lambda")
summary(lm1_best)

lm1_ols <- lm(numpos ~ dives + lnmass * island, data=df3)

df3$group <- NA
df3$group <- ifelse(df3$dives==0 & df3$island==0, 1, df3$group)
df3$group <- ifelse(df3$dives==0 & df3$island==1, 2, df3$group)
df3$group <- ifelse(df3$dives==1 & df3$island==0, 3, df3$group)
df3$group <- ifelse(df3$dives==1 & df3$island==1, 4, df3$group)
df3$group <- factor(df3$group)
df3$group <- factor(paste0(df3$dives, df3$island))

pars <- coef(lm1_best)
# pars <- coef(lm1_ols)

# Panel C
# ggplot(df3, aes(x=lnmass, y=numpos, fill=group, color=group)) +
# 	geom_point(pch=21, cex=1.5) +
# 	scale_fill_manual(values=pal2) +
# 	scale_color_manual(values=pal2) +
# 	theme_minimal() +
# 	geom_abline(slope=pars['lnmass'], intercept=pars[1]+pars['dives'], col=pal2[3]) +
# 	geom_abline(slope=pars['lnmass']+pars['lnmass:island'], intercept=pars[1]+pars['dives']+pars['island'], col=pal2[2]) +
# 	labs(x="Body mass (ln g)", y="Number of PSGs") +
# 	theme(legend.pos="none", text=element_text(size=8))
# ggsave("figs/psg_size.pdf", width=2.75, height=2.75)

pdf("figs/psg_size.pdf", width=3, height=2.75)
par(mgp=c(1.5,.5,0),mar=c(3,3,1,1),mex=.75,ps=8)
plot(numpos ~ lnmass, data = df3, col=pal2[group], bg=pal2[group], pch=21, xlab="Body mass (ln g)", ylab="Number of PSGs", bty='l')
abline(a=pars[1]+pars['dives'], b=pars['lnmass'], lty=2)
abline(a=pars[1], b=pars['lnmass'], lty=3)
abline(a=pars[1]+pars['dives']+pars['island'], b=pars['lnmass']+pars['lnmass:island'])
dev.off()

#-------------------------------------------------------------------------------
# overlap analysis with MCMCglmm
#-------------------------------------------------------------------------------

tmp <- as.matrix(dat_paml[, grep("P_aBSREL", names(dat_paml))])
colnames(tmp) <- gsub("P_aBSREL_", "", colnames(tmp))
idx <- combn(colnames(tmp), m=2)

dat_changes$absrel_overlap <- sapply(1:nrow(dat_changes), function(x) {
	id1 = dat_changes[x, 'sp1']
	id2 = dat_changes[x, 'sp2']
	sum(tmp[, id1] < 0.05 & tmp[, id2] < 0.05, na.rm=T)
})

ggplot(dat_changes, aes(x=phydist, y=absrel_overlap, color=group)) + geom_point() + scale_color_manual(values=pal2) + geom_smooth(method="lm", se=F)

dat_changes$plotgroup <- factor(paste0(dat_changes$fish_bg, dat_changes$island_bg))
levels(dat_changes$plotgroup) <- c('background', 'island-dwelling', 'continental\nplunge-diver', 'island\nplunge-diver')

datsum <- dat_changes %>% group_by(plotgroup) %>% summarize_if(is.numeric, mean)

ggplot(dat_changes, aes(x=absrel_overlap, y=plotgroup, color=plotgroup)) +
	geom_crossbar(data=datsum, aes(xmin = absrel_overlap, xmax = absrel_overlap, color=plotgroup),
                  size=.25, width = .5) +
	geom_jitter(height=.05, size=.5) +
	scale_fill_manual(values = setNames(pal2, levels(dat_changes$plotgroup))) +
	scale_color_manual(values = setNames(pal2, levels(dat_changes$plotgroup))) +
	theme_minimal() +
	# theme_bw() +
	labs(x="Pairwise overlap in PSGs", y="") +
	theme(legend.pos="none", panel.grid=element_blank(), text=element_text(size=8), axis.text.y=element_text(hjust=1, size=8))

ggsave("figs/dots_absrel_psgs.pdf", width=3, height=2.75)


prior0 <- list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)),R=list(V=1,nu=0.02))

dat_changes$grp <- factor(paste0(dat_changes$fish_bg, dat_changes$island_bg))

prior <- list(G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)),R=list(V=1,nu=0.02))

# set.seed(1980)
# fit_m_absrel_bm <- MCMCglmm(absrel_overlap ~ grp, data = dat_changes,
# 							random = ~sp1.sorted + sp2.sorted + node, nitt=1e5,
# 							burnin=0.25*1e5, prior=prior, ginverse=list(node=ainv))
# saveRDS(fit_m_absrel_bm, file="output/psg_overlap_mcmc.rds")

fit_m_absrel_bm <- readRDS(file="output/psg_overlap_mcmc.rds")

summary(fit_m_absrel_bm)

#-------------------------------------------------------------------------------
# significance testing for overlap
#-------------------------------------------------------------------------------

source("R/pMCMC.R")
pars = fit_m_absrel_bm$Sol
pars[,2] <- pars[,2]+pars[,1]
pars[,3] <- pars[,3]+pars[,1]
pars[,4] <- pars[,4]+pars[,1]
pairids <- combn(1:ncol(pars), m=2)
diffs = sapply(1:ncol(pairids), function(x) {
	id1=pairids[1,x]
	id2=pairids[2,x]
	pars[,id1]-pars[,id2]
})
pvals = p.adjust(pMCMC(as.data.frame(diffs)), method="fdr")
names(pvals) <- apply(combn(levels(dat_changes$grp), m=2), 2, paste0, collapse="-")
multcompView::multcompLetters(pvals)

tapply(dat_changes_sub$absrel_overlap, dat_changes_sub$group, mean)

# We found that the number of positively selected genes was significantly higher in small-bodied island kingfishers
# than larger island lineages. In continental species this body mass effect was not found.
# Plunge-diving species also had significantly more PSGs (P = xx) compared to non-plunge-diving lineages.

# "Uncorrected P values were reported for aBSREL; as we were testing the hypothesis that the short lived would share selective pressures, not all branches were individually of interest." (rockfish paper)

tr2 <- read.tree("data/paml/kingTree_abbrev.phy")
tr2 <- makeNodeLabel(tr2)
tr2.pruned <- drop.tip(tr2, "momMom")
m <- read.csv("data/absrel_match_nodes.csv")
m[,1] <- paste0("Node", m[,1])
m[,2] <- paste0("Node", m[,2])
nms <- c(tr2.pruned$tip.label, tr2.pruned$node.label)
idx <- tr2.pruned$edge[,2]
edgenames <- nms[idx]
# edgenames[1] <- "Node1"
edgenames[match(m[,2], edgenames)] <- m[,1]

# Calculate number of positively selected genes per branch
# num_pos <- apply(tmp, 2, function(x) sum(x<.05, na.rm=T))
# num_pos_scaled <- minmax(num_pos[edgenames], na.rm=T)

#-------------------------------------------------------------------------------
# Plot GO terms
#-------------------------------------------------------------------------------

df2 <- dat_enrich %>% gather(branch,genes,-term)

foc1 = intersect(islandFg, plungeFg)
foc2 = setdiff(plungeFg, islandFg)

df2$group <- ifelse(df2$branch %in% foc1, 'island', 'cont')
df2$lnmass <- setNames(brain$lnmass, brain$abbrev)[df2$branch]

df2$origin <- NA
df2$origin <- ifelse(df2$branch %in% c('corLeu','corVin','corCri'), '1', df2$origin)
df2$origin <- ifelse(df2$branch %in% c('ceyCya','ceyArg'), '2', df2$origin)
df2$origin <- ifelse(df2$branch %in% c('alcSem','alcQua','alcAtt'), '3', df2$origin)
df2$origin <- ifelse(df2$branch %in% c('megMax','megTor','megAlc','chlAma','chlInd','cerRud'), '4', df2$origin)

df2$term_abbrev = df2$term %>% 
  map_chr(~ str_split(.x, pattern=" ", simplify=TRUE) %>% 
            gsub("(.{5}).*", "\\1.", .) %>% 
            paste(., collapse=" "))

df2$Distribution <- factor(df2$group)
df2$Distribution <- plyr::revalue(df2$Distribution, c('cont'='Continental', 'island'='Insular'))

df2$PSGs <- df2$genes

conflicted::conflict_prefer("reorder", "stats")
conflicted::conflict_prefer("filter", "dplyr")

# Panel A

# color axis labels:
df2$color <- factor(ifelse(df2$term %in% unique(gp_11111$result[,"term_name"]), "black", "gray"))
df2$term_colored <- paste0("<span style=\"color: ", df2$color, "\">", df2$term_abbrev, "</span>")

df2 %>% filter(genes>=5 & !grepl("adapt. immun.|somat. diver.", df2$term_abbrev)) %>%
	ggplot(aes(x=reorder(branch, lnmass, mean), y=reorder(term_colored, PSGs, sum),
		size=PSGs, color=Distribution)) +
		geom_point() +
		scale_color_viridis() +
		scale_size_continuous(range = c(.5, 3)) +
		scale_color_manual(values = pal2[3:4], guide="none") +
		scale_x_discrete(position = "top") +
		theme_minimal() +
		facet_grid(~origin, scales="free_x", space="free_x") +
  		theme(strip.placement="outside",
  			  axis.text.x = element_text(angle=-45),
  			  axis.text.y = ggtext::element_markdown(),
  			  panel.spacing = unit(0, "lines"),
  			  text = element_text(size=7),
  			  legend.pos="bottom",
			  legend.margin=margin(0,0,0,0),
        	  panel.border = element_rect(fill = NA, color = "black", size=.25),
        	  panel.grid.major.x=element_blank()) +
  		labs(x = "Plunge-diving origin:", y = "GO term")

ggsave("figs/go_terms_islands.pdf", width=4, height=6.5)
```


## Liftoff

```sh
# y axis - plot pos selection estimates for each gene?

# conda create -n liftoff
# conda activate liftoff
conda install -c bioconda liftoff

liftoff

annot=~/uce-alcedinidae/annotations/gemoma/round2/final_annotation.gff
target=~/uce-alcedinidae/GCA_009819595.1_bMerNub1.pri_genomic.fna
ref=~/uce-alcedinidae/genomes/todChl/todChl.scaffolds.full_mask.fa

liftoff -g $annot $target $ref -p 24 -o todChl-to-merNub.liftoff.gff3
```

plot in R-

```{r}
library(ggtext)

gff <- read.delim("output/todChl-to-merNub.liftoff.gff3", header=F, comment.char="#")
gff <- gff[gff[,3]=="mRNA",]
gff$ID <- stringr::str_match(gff$V9, "ID=(.*?);")[,2]
gff$ID <- gsub(":", "_", gff$ID)

# 35=w, 36=z, 37=mt
lookup <- setNames(c(1:37), unique(gff$V1)[1:37])

gff$chrom <- plyr::revalue(gff$V1, lookup)

gff2 <- subset(gff, chrom %in% lookup)
gff2$chrom <- as.numeric(gff2$chrom)
# gff2$p <- runif(nrow(gff2))

gff2$genename <- dat_paml$finalgene[match(gff2$ID, dat_paml$file)]

gff2$P <- dat_paml$FDR_M2[match(gff2$ID, dat_paml$file)]

gff3 <- subset(gff2, !is.na(P))

# gff3$P <- ifelse(gff3$P==1, .999, gff3$P)
gff3$P <- ifelse(gff3$P==0, min(gff3$P[gff3$P>0]), gff3$P)


# par(bg="gray25")
# manhattan(gff3, chr="chrom", bp="V4", p="P", snp="genename", highlight = genelist_adap_conv, logp=TRUE, chrlabs=c(1:34, 'w','z'))

gff3$chrom <- plyr::revalue(as.character(gff3$chrom), c('35'='W', '36'='Z'))
gff3$chrom <- factor(gff3$chrom, levels=c(1:34, 'W', 'Z'))

data_cum <- gff3 %>% 
  group_by(chrom) %>% 
  summarise(max_bp = max(V4)) %>% 
  mutate(bp_add = dplyr::lag(cumsum(max_bp), default = 0)) %>% 
  select(chrom, bp_add)

gwas_data <- gff3 %>% 
  inner_join(data_cum, by = "chrom") %>% 
  mutate(bp_cum = V4 + bp_add)

axis_set <- gwas_data %>% 
  group_by(chrom) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  dplyr::filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)

# sig <- 5e-8
sig <- 0.05

gwas_data$foc <- factor(ifelse(gwas_data$genename %in% genelist_adap_conv & gwas_data$P < .05, 1, 0))

manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P), 
                                  color = factor(chrom))) +
  # geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(pch=16, size=0.5) +
  scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  # scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chrom)))) +
  scale_color_manual(values = rep(grey(c(.25,.75),alpha=0.25), unique(length(axis_set$chrom)))) +
  # scale_fill_manual(values = c('white', 'red')) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_point(data=subset(gwas_data, foc==1), pch=16, col='blue', size=1) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(size=6),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 4, vjust = 0.5)
  )
# manhplot
ggsave(manhplot, file="figs/manhat.pdf", width=4.5, height=1)
```
