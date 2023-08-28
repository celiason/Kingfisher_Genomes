#' Convert gene-wise fasta files to species-wise fastas
#' 
#' Function converts gene-wise fasta files to species-wise fastas
#' (e.g., multiple sp. per one gene file)
#' Relies on faSomeRecords
#' Download at http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
#' Put in /bin folder
#' Add execute permissions chmod a+x faSomeRecords.
#' 
#' @param files list of FASTA files
#' @param cores number of cores
#' @param outpath path to output FASTAs
#' @param suffix suffix to remove from file names when naming species in out fasta (e.g., ".masked.fa")
#' @param picks numeric vector of which sub sequences to output (e.g, c(1,3,50))
#' @param lookup named character vector (names = old sequence names, vector iteself = new sequence names)
#' @param namesep unused? e.g., "|" for ">todChl|Mc1r" format
#' 
genewise2sppwise <- function(files, cores=1, outpath=".", suffix=NULL, lookup=NULL, picks=NULL, contig=NULL, group=FALSE, namesep=NULL, force=FALSE) {
	require(parallel)
	require(stringr)
	require(pbapply)
	seqs <- pblapply(files, readLines)
	seqstarts <- lapply(seqs, grep, pattern=">")
	if (length(unique(sapply(seqstarts, length))) != 1) {
		stop("Multifastas do not all have the same number of genes.")
	}
	nseq <- length(seqstarts[[1]])
	nsamp <- length(seqs)
	nms <- gsub(">", "", seqs[[1]][seqstarts[[1]]])	
	if (!is.null(lookup)) {
		if (any(duplicated(lookup))) {
			warning("Some names are duplicated in `lookup`, creating unique IDs but might want to check.")
			oldnames <- names(lookup)
			lookup <- make.unique(as.character(lookup))
			names(lookup) <- oldnames
		}
		nms <- lookup[nms]
	}
	names(nms) <- NULL
	# Replace hyphens and dots with underscores to make names cleaner
	nms <- gsub("(-|\\.)", "_", nms)
	# Sub sequences
	if (is.null(picks)) {
		picks <- 1:nseq
	}
	# Parse sequences by gene and output as new files
	for (i in picks) {
		outfile <- paste0(outpath, "/", nms[i], ".fa")
		if (file.exists(outfile) & !force) {
			stop("Files exists, consider using `force=TRUE`.")
		} else if (force) {
			file.remove(outfile)
		}
		for (j in 1:nsamp) {
			start <- seqstarts[[j]][i]
			if (i == length(seqstarts[[j]])) {
				end <- length(seqs[[j]])
			} else {
				end <- seqstarts[[j]][i+1] - 1
			}
			if (!is.null(suffix)) {
				spp <- gsub(suffix, "", basename(files[j]))
			} else {
				spp <- basename(files[j])
			}
			bits <- seqs[[j]][start : end]
			bits <- gsub(">.*?$", paste0(">", spp), bits)
			cat(bits, file=outfile, append=TRUE, sep="\n")
		}
	}
}

#---
# Testing zone
# files <- list.files("uces", pattern="masked.fa$", full=T)[1:3]
# meta <- read.table("todChl_uce5k_sorted.bed")
# lookup <- setNames(meta$V4, paste0(meta$V1, ":", meta$V2, "-", meta$V3))
# lookup[1:3]
# dir.create("test_rename")
# genewise2sppwise(files = files, outpath = "test_rename", suffix = ".uce5k.masked.fa", lookup = lookup)
#---

