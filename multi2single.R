#' Convert multi-line to single-line fastas
#' 
#' @param fastapath path to FASTA files
#' @param outpath output path
#' 
#' @export
#' 
# fasta <- files[1]
multi2single <- function(fasta, fastapath=NULL, outpath=NULL) {
	# lots of files
	if (!is.null(fastapath)) {
		files <- list.files(fastapath, pattern=".fa$", full=TRUE)
		seqs <- lapply(files, readLines)
	# single fasta file
	} else {
		seqs <- list(readLines(fasta))
	}
	for (i in seq_along(seqs)) {
		# i=1
		firstlines <- grep(">", seqs[[i]])
		lastlines <- c(grep(">", seqs[[i]])[-1] - 1, length(seqs[[i]]))
		newseq <- lapply(seq_along(firstlines), function(x) {
			c(seqs[[i]][firstlines[x]], paste0(seqs[[i]][(firstlines[x]+1) : lastlines[x]], collapse=""))
		})
		newseq <- unlist(newseq)
		nms <- newseq[seq(1, length(newseq)-1, by=2)]
		seqs[[i]] <- setNames(as.list(newseq[seq(2, length(newseq), by=2)]), nms)
	}
	if (!is.null(fastapath)) {
		setNames(seqs, basename(files))
	} else if (length(seqs)==1) {
		return(seqs[[1]])
	}
}
