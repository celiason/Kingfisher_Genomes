#' Function to process GFF files
#' 
#' @param gff sorted GFF file (e.g., output from MAKER), or CSV file
#' @param bamfile path to a BAM file (must be a bam BAI index file in same dir)
#' @param ref path to reference genome in fasta format
#' 
#' @export
#' 
#' Outputs: CDS from BAM files
#' 
#' @examples \dontrun{
#' 
#' # Example for one BAM file
#' bam2fasta(gff = "todChl_rnd3.all.maker.noseq.sorted.gff", bamid = "SRR9644504", ref="King_genome/kingfisher.fasta")
#' 
#' # Get CDS for a batch of BAM files
#' bamfiles <- list.files(pattern=".bam$")
#' bamids <- stringr::str_match(bamfiles, "(SRR\\d+)\\.bam")[, 2]
#' library(parallel)
#' mclapply(bamids, bam2fasta, gff="todChl_rnd3.all.maker.noseq.sorted.gff", ref="King_genome/kingfisher.fasta", mc.cores = 4)
#' bam2fasta(gff = "todChl_rnd3.all.maker.noseq.sorted.gff", bamid = "SRR9644504", ref="King_genome/kingfisher.fasta")
#' 
#' # Get reference CDS
#' bam2fasta(gff = "todChl_rnd3.all.maker.noseq.sorted.gff", ref="King_genome/kingfisher.fasta")
#' 
#' }
#' 
#' NB: run R on server with *capital* R
#' 
#' # TODO: parallelize this!
#' TODO: stop outputting calls.vcf files, etc. - figure out how to use STDIN

# It might be better to create a BED file and then use that to get sequences in Mpileup..
# files <- list.files("bam", pattern=".bam$", full=T)
# picks <- sapply(c("SRR10349430", "SRR10350663", "SRR10351892", "SRR10352149"), grep, files)
# gff="data/kingfisher_sensory_cds_gff3.csv"
# ref="reference_genome/kingfisher.fasta"
# bamfile=files[2]
# outpath="extracted_seqs/sensory"

# alias bcftools=bcftools-1.9/bcftools
# alias samtools=samtools-1.9/samtools

# TESTING ZONE
# bamfile="bam/SRR10349430.bam"
# can be 1 bam file, or list of them
# ref="reference_genome/kingfisher.fasta"
# outpath=NULL
# gff="data/kingfisher_sensory_cds_gff3.csv"

# bamfile="bam/kingfisher10x.bam"
# ref="reference_genome/kingfisher.fasta"
# region="132:829255" # exon 1 of tau
# gff="data/Mapt_gff.csv"
# outpath=NULL

bam2fasta <- function(gff, ref, bamfile=NULL, outpath=NULL) {

	require(seqinr)
	require(stringr)

	if (is.null(outpath)) {
		outpath <- getwd()
	}

	# outpath="extracted_seqs/sensory"

	if (!file.exists(outpath)) {
		system(paste0("mkdir ", outpath))
	}
	
	# Add slash to output directory for later:
	if (!grepl("\\/$", outpath)) {
		outpath <- paste0(outpath, "/")
	}

	fileext <- stringr::str_extract(gff, "[a-z]+$")

	if (fileext == "csv") {
		cat("Reading GFF from CSV file...\n")
		temp <- read.csv(gff)
	}

	if (fileext == "gff") {
		cat("Reading GFF file...\n")
		temp <- read.delim(gff, sep="\t", header=FALSE, stringsAsFactors=FALSE)
		names(temp) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
		temp <- subset(temp, feature=="CDS")
	}

	if (any(names(temp) == "gene")) {
		l <- split(temp, temp$gene)
	} else {
		l <- split(temp, stringr::str_match(temp$attribute, "ID=(.*?)(?=:)")[, 2])	
	}

	# Check to make sure all strands are in the same direction:
	junk <- sapply(1:length(l), function(x) {
		length(unique(l[[x]]$strand))
	})
	
	# TODO: check this, might need to figure out what to do if + AND - exons
	if (length(unique(junk))!=1) {
		stop("Not all genes have exons on the same strand")
	}

	# Checks
	# nm <- "maker-91627-augustus-gene-0.1-mRNA-1"
	# tmp <- temp[grep(nm, temp$attribute, fixed=TRUE), ]
	# sum(tmp$end - (tmp$start - as.numeric(tmp$frame)))
	# 492 bp .. cool! matches AA sequence length

	pb <- txtProgressBar(min = 0, max = length(l), initial = 0, style = 3)

	# for (i in 1:3) {
	for (i in 1:length(l)) {

		# i=1

		# Subset the GFF file
		gffsub <- l[[i]]
		nm <- names(l)[i]
		
		seqs <- list()

		if (!is.null(bamfile)) {
			bamid <- str_extract(basename(bamfile), ".*?(?=.bam)")
		}
		
		# Process individual genes
		# cat("\n###### Processing gene ", i, "/", length(l), " ###### \n\n", sep='')

		setTxtProgressBar(pb,i)

# j=2
		for (j in 1:nrow(gffsub)) {
		# j=11
			region <- paste0(gffsub$seqname[j], ":", gffsub$start[j], "-", gffsub$end[j])
		
			if (!is.null(bamfile)) {
				
				# Extract sequence from BAM file using given GFF coordinates
				# TODO check q/Q/d parameters
				# NB: 2>/dev/null silences output notifications and warnings!
				system(paste0("bcftools mpileup -q 30 -Q 30 -d 100 -r ", region, " -Ou -f ", ref, " ", bamfile, " 2>/dev/null | bcftools call -mv -Oz -o ", bamid, "_calls.vcf.gz 2>/dev/null"))

# zcat SRR10349430_calls.vcf.gz

				# Normalize indels (CHECK THIS)
				# -m +any - join biallelic sites into multiallelic records (+). SNPs and indels should be merged into a single record (any)
				# -Oz output as compressed VCF (-z)
				system(paste0("bcftools norm -f ", ref, " -m +any -Oz ", bamid, "_calls.vcf.gz -o ", bamid, "_calls.norm.vcf.gz 2>/dev/null"))
				
				# TODO? Filter adjacent indels within 5bp??
				
				# Index VCF file
				system(paste0("bcftools index -f ", bamid, "_calls.norm.vcf.gz"))

# zcat SRR10349430_calls.norm.vcf.gz

				# Extract consensus sequences
				# -H 1: haplotype select first allele, regardless of phasing
				myseq <- system(paste0("samtools faidx ", ref, " ", region, " | bcftools consensus -H 1 ", bamid, "_calls.norm.vcf.gz 2>/dev/null"), intern=TRUE)
				# myseq <- system(paste0("samtools faidx ", ref, " ", region, " | bcftools consensus ", bamid, "_calls.norm.vcf.gz 2>/dev/null"), intern=TRUE)				
			} else {
				
				# If no bam file, just extract CDS from reference genome
				myseq <- system(paste0("samtools faidx ", ref, " ", region), intern=TRUE)
			
			}
			
			seqs[[j]] <- paste0(myseq[-1], collapse="")
			# seqs[[j]] <- c(paste0(">", bamid, "|", nm, "|", region), myseq)
		}

		# Convert sequence fragments into single CDS
		
		# Check if need to do reverse complement-
		if (all(gffsub$strand=="-")) {
			seqs <- lapply(seqs, function(x) {
				c2s(rev(comp(s2c(x))))	
			})
		}
		rawc <- paste0(unlist(seqs), collapse="")
		
		# Reverse complement
		# if (all(gffsub$strand=="-")) {
		# 	rawc <- c2s(rev(comp(s2c(rawc))))
		# }

		# TODO: add translate option at some point?
		# c2s(translate(s2c(rawc)))
		# DONE do something about strand (-/+)
		# what if some -, some +??		
		# Works if all neg..
		# See added check above that stops the script if not all same strand

		# Convert to upper case
		rawc <- toupper(rawc)

		# Format as fasta (60 wide columns)
		starts <- seq(1, nchar(rawc), 60)
		stops <- seq(1, nchar(rawc), 60) + 59
		rawc <- substring(rawc, first = starts, last = stops)
		
		# Combine single fasta files to multifasta
		if (!is.null(bamfile)) {
			# Append gene name and output
			outfile <- paste0(outpath, bamid, ".cds.fa")
			# TODO check this works correctly
			# if (file.exists(outfile)) {
			# 	stop("File already exists, creating copy")
			# 	outfile <- paste0(outpath, bamid, ".copy.cds.fa")
			# }
			cat(paste0(">", bamid, "|", nm), rawc, file=outfile, append=TRUE, sep="\n")
			# Cleanup
			system(paste0("rm ", bamid, "_calls*"))
		} else {
			# TODO this might not be needed now? Since we input a GFF from the reference
			cat(paste0(">reference|", nm), rawc, append=TRUE, sep="\n", file = paste0(outpath, "reference.cds.fa"))
		}
	}
}
