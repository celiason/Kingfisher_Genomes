#' Remove stop codons
#' 
#' @param seqs list of DNA sequences in seqinr(?) format (e.g., `seqinr::read.fasta`)
#' 
cutstops <- function(seqs, fastapath=NULL) {
	if (!is.null(fastapath)) {
		seqs <- multi2single(fastapath)
	}
	# Check if codons
	id <- sapply(seqs, function(x) {
		all((sapply(x, nchar) %% 3) == 0)
	})
	if (!all(id)) {
		stop("Sequences do not appear to be codon alignments, try aligning (e.g., with MAFFT)")
	}
	for (i in seq_along(seqs)) {
		for (j in seq_along(seqs[[i]])) {
			tmp <- substring(seqs[[i]][[j]], first=seq(1, nchar(seqs[[i]][[j]]), by=3), last=seq(1, nchar(seqs[[i]][[j]]), by=3) + 2)
			seqs[[i]][[j]] <- paste0(gsub("TAG|TAA|TGA", "---", tmp), collapse="")
		}
	}
	return(seqs)
}
