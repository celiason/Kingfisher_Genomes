#' Function to align raw reads to reference genome
#' 
#' TODO: add description
#' 
#' @param ref path to reference genome (genome will be indexed if not already)
#' @param reads path to raw reads file (PE)
#' @param cores how many cores to use for BWA MEM
#' @param ram how much RAM (GB) to allot to BAM file output
#' 
#' @export
#' 
alignReads <- function(ref, reads, cores=48, ram=150, suffix=NULL, test=FALSE, force=FALSE) {
	# Paths to programs
	# bwa="/home/FM/celiason/bwa/bwa"
	# fastp="/home/FM/celiason/fastp"
	# Setup
	oldwd <- getwd()
	# id <- 
	id <- gsub(".fa(.*?)$", "", basename(ref))
	from <- substr(basename(reads), 1, 6)
	to <- substr(basename(ref), 1, 6)
	prefix <- paste0(from, "-to-", to)
	# Index reference genome
	setwd(paste0("genomes/", to))
	if (!file.exists(paste0(basename(ref), ".bwt"))) {
		system(paste0("bwa index ", basename(ref)))
	}
	setwd(oldwd)
	# Setup name	
	if (!is.null(suffix)) {
		prefix <- paste0(prefix, suffix)
	}
	if (!force & file.exists(paste0("alignments/", prefix, ".bam"))) {
		stop("BAM file already exists! Please add or change `suffix` argument.")
	}
	# Align raw reads to reference genome
	run <- paste0("fastp -i ", reads, " --interleaved_in --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --stdout -h alignments/", prefix, ".html | sed -E 's/^((@|\\+)SRR[^.]+\\.[^.]+)\\.(1|2)/\\1/' | bwa mem -p -t ", cores, " ", ref, " - | samtools view -Sb - | samtools sort -m ", ram, "G > alignments/", prefix, ".bam")
	# run
	if (test) {
		run
	} else {
		system(run)
	}
}
