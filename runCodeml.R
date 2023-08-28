#' Function to run codeml within R
#' 
#' @param phy path to phylogenetic tree
#' @param model which model to use (see XX for details)
#' @param fasta path to FASTA codon alignment file
#' @param delete whether to force analysis (delete original) and mask early stop codons
#' @param outpath output path for results
#' @param silent whether to output run progress in console
#' @param M0path path for M0 results to use for branch lengths
#' @param droptest
#' @param fg vector with foreground species to drop
#' 
#' @author Chad M. Eliason
#' 
#' testing zone:
#' fasta=files[1]
#' model="M0"
#' phy="/home/FM/celiason/uce-alcedinidae/paml/kingtree_nolabels.phy"
#' 
runCodeml <- function(fasta, model=c("M0", "M1a", "M2a", "M7", "M8", "MA", "MAnull"), phy=NULL, fixedbl=FALSE, droptest=FALSE, outpath=NULL, silent=FALSE, M0path=NULL, fg=NULL) {

# TODO consider changing model names? (M1, M2, BS1, BS2)

	# Load required packages
	require(stringr)
	require(treeio)
	require(seqinr)
	require(ape)

	oldwd <- getwd()
	
	prefix <- gsub("\\..[^\\.]*$", "", basename(fasta))
	
	model <- match.arg(model)

	if (droptest & is.null(fg)) {
		stop("No foreground branches specified")
	}
	
	ctl <- readLines(paste0(model, ".ctl"))

	# For model M0 we are using phylogeny without any nodes labeled
	if (model == "M0") {
		ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = kingtree_nolabels.phy"))
		warning("Ignoring `phy` argument. Using kingtree_nolabels.phy in working directory.")
	}

	# Get path of M0 file if using that tree
	if (is.null(M0path)) {
		M0path <- paste0("M0_results/", prefix)
	} else {
		M0path <- paste0(M0path, "/", prefix)
	}
	
	# Setup outut path for files
	if (is.null(outpath)) {
		outpath <- paste0(model, "_results/", prefix)
	} else {
		outpath <- paste0(outpath, "/", prefix)
	}

	# For models M1a, M2a, and free-ratio we are using branch lengths estimated under M0 model:
	if (model %in% c("M1a", "M2a", "M7", "M8", "free", "MA", "MAnull")) {
		raw <- readLines(paste0(M0path, "/M0_mlc"))
		textphy <- raw[grep("tree length", raw)[1]+4]
	}

	if (fixedbl & is.na(textphy)) {
		warning("No branch lengths found in M0 folder...moving to next sequence")
		return(NA)
	}

	# Create output path, make sure not to overwrite
	if (!dir.exists(outpath)) {
		dir.create(outpath, recursive=TRUE)
	} else if (dir.exists(outpath) & !delete) {
		stop(paste0("Path ", outpath, " already exists, consider deleting?"))
	}
	
	# Move into folder where codeml results will be output
	setwd(outpath)

	# Load alignment and tree
	if (fixedbl) {		
		phy <- read.tree(text=textphy)
	} else {
		phy <- read.tree(phy)
	}
	
	# Load alignment
	aln <- read.fasta(fasta)

	# Replacing names in the control file
	# ctl <- str_replace(ctl, "treefile = .*", paste0("treefile = ", phy))
	# ctl <- str_replace(ctl, "seqfile = .*", paste0("seqfile = ", fasta))

	# Modify the M0 trees and add annotations for nodes (e.g., "#1")
	if (model %in% c("MA", "MAnull")) {		
		phy$node.label <- NULL
		phy <- treeio::label_branch_paml(phy, 44, "#1")
		phy <- treeio::label_branch_paml(phy, 56, "#1")
		phy <- treeio::label_branch_paml(phy, 59, "#1")
		phy <- treeio::label_branch_paml(phy, 51, "#1")
		phy$node.label <- gsub(" ", "", phy$node.label)
		phy$tip.label <- gsub(" ", "", phy$tip.label)
	}
	
	# Drop foreground species if doing drop test (see REF)
	if (model %in% c("M1a", "M2a") & droptest) {
		phy <- drop.tip(phy, fg)
		aln <- aln[!(names(aln) %in% fg)]
	}

	# Replace name in codeml control file
	ctl <- str_replace(ctl, "treefile = .*$", "treefile = ./tree.phy")
	ctl <- str_replace(ctl, "seqfile = .*", "seqfile = ./alignment.fa")

	# Keep species with < 50% Ns in sequence
	# aln <- aln[!grepl("n", aln)]
	# drops <- sapply(aln, function(x) sum(x=="n"))/length(aln[[1]]) > 0.5
	# ntips <- Ntip(phy)
	# aln <- aln[-drops]
	# aln <- lapply(aln, function(x) gsub("n", "-", x))
	# phy <- keep.tip(phy, names(aln))

	# Output tree and alignment
	treeio::write.tree(phy, file="tree.phy")
	write.fasta(aln, names=names(aln), file="alignment.fa")

	# Output ctl file in current directory
	cat(ctl, file="temp.ctl", sep="\n")

	# Run codeml
	# silent=FALSE
	system("codeml temp.ctl", ignore.stdout=silent)
	
	# Set back to original dir
	setwd(oldwd)
}

# M1 ll = -6320.2
# M2 ll = -6293.48

parseCodeml <- function(respath) {
	list.files(respath, pattern="mlc")
	res <- readLines("mlc")
	omega <- as.numeric(stringr::str_extract(grep("omega", res, value=T), "[\\d\\.]+"))
	setwd(oldwd)
	omega
}
