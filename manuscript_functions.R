# This function gets ancestral DNA sequences and calculates the number of convergent
# and divergent changes for each pair of species (from baseml model)
getChanges_tips_baseml <- function(res) {
	tr <- res@phylo	
	idx <- combn(tr$tip, m=2)
	changes <- sapply(1:ncol(idx), function(x) {
		#x=1
		nodeID <- getMRCA(tr, idx[, x])
		tips <- idx[, x]
		tipID <- match(idx[, x], tr$tip)
		anc <- as.character(res@anc_seq[[as.character(nodeID)]])
		seq1 <- as.character(res@tip_seq[[tips[1]]])
		seq2 <- as.character(res@tip_seq[[tips[2]]])
		np <- sum(seq1==seq2 & seq1==anc & seq2==anc, na.rm=T)  # parallel changes
		nc <- sum(seq1==seq2 & seq1!=anc & seq2!=anc, na.rm=T)  # convergent changes
		nd <- sum(seq1!=seq2, na.rm=T)  # divergent changes
		c(par = np, conv = nc, div = nd)
	})
	fishtips <- c('corLeu','corVin','corCri','ceyPus','ceyWeb','ceyAzu','ceyCya','ceyArg','alcQua','alcAtt','megMax','megTor','megAlc','chlAma','chlAme','chlAen','chlInd','halCor','halPil')
	bothfish <- ifelse(idx[1, ] %in% fishtips & idx[2, ] %in% fishtips, 1, 0)
	changes <- as.data.frame(t(changes))
	changes$fish <- factor(bothfish)
	changes$spp1 <- idx[1, ]
	changes$spp2 <- idx[2, ]
	changes$phydist <- as.numeric(as.dist(cophenetic(tr)))
	return(changes)
}

# This function gets ancestral AA sequences and calculates the number of convergent
# and divergent changes for each pair of species
getChanges_tips <- function(res) {
	tr <- res@phylo	
	idx <- combn(tr$tip, m=2)	
	changes <- sapply(1:ncol(idx), function(x) {
		nodeID <- getMRCA(tr, idx[, x])
		tips <- idx[, x]
		tipID <- match(idx[, x], tr$tip)
		anc <- as.character(res@anc_seq[[as.character(nodeID)]])
		seq1 <- as.character(res@tip_seq[[tips[1]]])
		seq2 <- as.character(res@tip_seq[[tips[2]]])
		aa1 <- translate(seq1)  # translate
		aa2 <- translate(seq2)
		aa3 <- translate(anc)
		aa1 <- ifelse(aa1=="X", NA, aa1)  # fix X amino acids
		aa2 <- ifelse(aa2=="X", NA, aa2)
		aa3 <- ifelse(aa3=="X", NA, aa3)
		np <- sum(aa1==aa2 & aa1==aa3 & aa2==aa3, na.rm=T)  # parallel changes
		nc <- sum(aa1==aa2 & aa1!=aa3 & aa2!=aa3, na.rm=T)  # convergent changes
		nd <- sum(aa1!=aa2, na.rm=T)  # divergent changes
		c(par = np, conv = nc, div = nd)
	})
	fishtips <- c('corLeu','corVin','corCri','ceyPus','ceyWeb','ceyAzu','ceyCya','ceyArg','alcQua','alcAtt','megMax','megTor','megAlc','chlAma','chlAme','chlAen','chlInd','halCor','halPil')
	bothfish <- ifelse(idx[1, ] %in% fishtips & idx[2, ] %in% fishtips, 1, 0)
	changes <- as.data.frame(t(changes))
	changes$fish <- factor(bothfish)
	changes$spp1 <- idx[1, ]
	changes$spp2 <- idx[2, ]
	changes$phydist <- as.numeric(as.dist(cophenetic(tr)))
	return(changes)
}

getChanges_edges <- function(res) {
	tr <- res@phylo
	idx <- combn(1:Nedge(tr), m=2)
	all_seq <- c(res@tip_seq, res@anc_seq)
	aa <- apply(sapply(all_seq, as.character), 2, translate)
	aa <- gsub("X", NA, aa)
	n1 <- tr$edge[idx[1, ], 1]
	n2 <- tr$edge[idx[1, ], 2]
	n3 <- tr$edge[idx[2, ], 1]
	n4 <- tr$edge[idx[2, ], 2]
	changes <- sapply(1:ncol(idx), function(x) {
		nc <- sum(aa[,n2[x]]==aa[,n4[x]] & aa[,n1[x]]!=aa[,n3[x]], na.rm=T)  # convergent changes
		nd <- sum(aa[,n2[x]]!=aa[,n4[x]] & aa[,n1[x]]==aa[,n3[x]], na.rm=T)  # divergent changes
		np <- sum(aa[,n1[x]]==aa[,n3[x]] & aa[,n2[x]]==aa[,n4[x]], na.rm=T)  # parallel changes
		c(par = np, conv = nc, div = nd)
		# c(conv = nc, div = nd)
	})
	return(t(changes))
}

random.effect.sorting<-function(data,counter.max=1e6,seed=Sys.time()){
	if(class(data)!="matrix"){data<-as.matrix(data)}
	check<-unique(c(levels(as.factor(data[,1])),levels(as.factor(data[,2])))) ##arrange the orders of sp1 and sp2 so that they are as evenly distributed as possible
	ord=2
	counter=0
	half=dim(data)[1]/2
	half.samp=sample(1:dim(data)[1],half)
	data[half.samp,]<-data[half.samp,c(2,1)]
	set.seed(seed)
	while(ord[1]>1 && counter < counter.max){
		counter=counter+1
		t1<-table(data[,1])
		t2<-table(data[,2])
		t1<-t1[match(check,names(t1))]
		names(t1)<-check
		t2<-t2[match(check,names(t2))]
		names(t2)<-check
		#int<-data.frame(check,as.vector(t1),as.vector(t2))
		#colnames(int)<-c("check","t1","t2")	
		t1[is.na(t1)] <- 0
		t2[is.na(t2)] <- 0
		diff<-t1-t2
		ord<-sort(abs(diff),decreasing=TRUE)
		for(j in 1:4){
			if(ord[j]>1){
				rowx<-sample(seq(1:(dim(data[which(data[,1]==names(ord[j])| data[,2]==names(ord[j])),])[1])),(ord[j]/2))
				data[which(data[,1]==names(ord[j])| data[,2]==names(ord[j])),][rowx,]<-data[which(data[,1]==names(ord[j])| data[,2]==names(ord[j])),][rowx,c(2,1)]
			}
		}	
	}	
	if(counter==counter.max){print("hit counter max")}
	return(list(sorted.data=data,asymmetry=ord))
}	

# Function to scale from min (0)-max (1)
minmax <- function(x, ...) {
  (x-min(x, ...)) / (max(x-min(x, ...), ...))
}

se <- function(x) {
	sd(x)/sqrt(length(x))
}

# x = scientific names
# x=gsub(" ", "_", meta$Accepted_sciname)
# output genome species abbreviates (e.g., "homSap" for Homo sapiens)
sppAbbrev <- function(x) {
	tmp <- str_split(x, " |_")
	paste0(tolower(substring(sapply(tmp, "[[", 1), 1, 3)),
		   toupper(substring(sapply(tmp, "[[", 2), 1, 1)),
		   substring(sapply(tmp, "[[", 2), 2, 3))
}

# Trait rate equal splits (see Cooney et al. 2019 Nat. Comm)
tres <- function(phy, node=FALSE) {
	require(ape)
	Ntip <- length(phy$tip.label)
	nnode <- Nnode(phy)
	rootnd <- Ntip + 1L
	if (node) {
        ids <- c(1:Ntip, (Ntip+1) : (Ntip + nnode))
    } else {
        ids <- c(1:Ntip)
    }
	# ids <- c(1:Ntip)
	es <- numeric(length(ids))
	for (k in 1:length(ids)) {
		# k=4
	    # get lineage
	    lin <- c(ids[k], phangorn::Ancestors(phy, ids[k]));
	    # drop root node
	    lin <- lin[-length(lin)]
	    if (length(lin) == 0) {
	        es[k] <- 0
	        next
	    }
	    # set indices
	    inds <- seq(1, length(lin), 1)
	    # get edge lengths
	    els <- phy$edge.length[match(lin, phy$edge[,2])]
		# downweight rootward edges
		w <- (1/ (2^(inds-1)))
	    # compute rates following Cooney et al. (2019) Nat. Comm
		es[k] <- sum(els * w / w)
	    # es[k] <- 1/sum(els * (1 / 2^(inds-1)))
	}
	if (node) {
		nms <- c(phy$tip.label, (Ntip(phy)+1) : (Ntip(phy) + Nnode(phy)) )
    	names(es) <- nms[ids]
    	es <- es[phy$edge[, 2]]
    } else {
    	names(es) <- phy$tip.label	
    }
	return(es)
}

# Plotting the result of PSMC analysis using R.
# By Shenglin Liu, Apr 3, 2016.

# https://datadryad.org/bitstream/handle/10255/dryad.126825/plotPsmc.r?sequence=1

##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation

psmc.resultOLD<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,"temp.psmc.result")
	theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
	N0<-theta0/4/mu/s
	
	write(RS,"temp.psmc.result")
	a<-read.table("temp.psmc.result")
	Generation<-as.numeric(2*N0*a[,3])
	Ne<-as.numeric(N0*a[,4])
	
	file.remove("temp.psmc.result")
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	data.frame(YearsAgo,Ne)
}


##-------Rescale the ith iteration result of PSMC, and make ready for plotting
# file: result file from PSMC
# i.iteration: the ith iteration
# mu: mutation rate
# s: bin size
# g: years per generation

# file="data/psmc/todChl_diploid.psmc"
# file="data/psmc/todChl_combined.psmc"
# mu=4.6e-9
# g=5
# i.iteration=30
# s = 100

psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	# X<-readLines(f)

	START<-grep("^RD", X)
	END<-grep("^//", X)

	idx <- seq(i.iteration+1, length(START), by=i.iteration+1)

	START <- START[idx]
	END <- END[idx]

	res <- lapply(1:length(START), function(i) {
		# i=1
		tmp <- X[START[i] : END[i]]
		TR<-grep("^TR",tmp,value=TRUE)
		RS<-grep("^RS",tmp,value=TRUE)
		theta0 <- as.numeric(stringr::str_match(TR, "TR\\t(.*?)\\t")[,2])
		N0 <- theta0/4/mu/s
		# write(RS,"temp.psmc.result")
		# a<-read.table("temp.psmc.result")
		# a[,3]
		Generation <- 2 * N0 * as.numeric(do.call(rbind, strsplit(RS, "\t"))[, 3])
		# Generation<-as.numeric(2*N0*a[,3])
		# Ne<-as.numeric(N0*a[,4])
		Ne <- N0 * as.numeric(do.call(rbind, strsplit(RS, "\t"))[, 4])
		# file.remove("temp.psmc.result")
		n.points<-length(Ne)
		YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
			Generation[n.points])*g
		Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
			Ne[n.points])
		data.frame(YearsAgo,Ne)
		# plot(data.frame(YearsAgo,Ne), log='y', type='s')
	})

	res
}



##-------Turn "ms" commandline into history, and make ready for plotting
# N0: diploid population size

history<-function(command,N0,g=1)
{
	X<-unlist(strsplit(command," "))
	X<-matrix(X,3,length(X)/3)
	Generation<-as.numeric(X[2,])*4*N0
	Generation<-c(10,Generation)
	Ne<-as.numeric(X[3,])*N0
	Ne<-c(N0,Ne)
	
	n.points<-length(Ne)
	Generation<-c(Generation,2*Generation[n.points])
	Ne<-c(Ne,Ne[n.points])
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*g
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	data.frame(YearsAgo,Ne)
}


##-------Plot population by population; real history can be added in ms commandline
# keywords: regular expression; beginning format of the psmc files' names
# labels: names of the output plots
# command: e.g., "-eN 0.01 0.1 -eN 0.06 1 -eN 0.2 0.5 -eN 1 1 -eN 2 2"

plotPsmc.popwise<-function(keywords, labels,
	command=NA, N0=25000,
	save.as="png", height=7, width=12,
	i.iteration=25, mu=1e-8, s=100, g=1,
	ylim=c(0,100000), xlim=c(200,500000),
	col=rep("red",length(keywords)), col.hist="grey")
{
	n<-length(keywords)
	files<-grep("\\.psmc$",dir(),value=TRUE)
	dev.new(height=height,width=width)
	for(i in 1:n)
	{
		subfiles<-grep(paste("^",keywords[i],sep=""),files,value=TRUE)
		n.sub<-length(subfiles)
		plot(1,1,
			ylim=ylim,xlim=xlim,
			log="x",type="n",
			ylab="Ne",xlab="Generations ago")
		for(i.sub in 1:n.sub)
		{
			lines(psmc.result(subfiles[i.sub],i.iteration,mu,s,g),
				type="l",col=col[i],lwd=1)
		}
		if(!is.na(command))
		{
			lines(history(command,N0,g),
				type="l",col=col.hist,lwd=2)
		}
		savePlot(filename=paste(labels[i],save.as,sep="."),type=save.as)
	}
	dev.off()
}


##-------Plot all populations together

plotPsmc.allPops<-function(keywords, label, legend.names,
	save.as="png", height=7, width=12,
	i.iteration=25, mu=3.7e-8, s=100, g=2,
	ylim=c(0,60000), xlim=c(200,500000),
	col=rainbow(length(keywords)))
{
	n<-length(keywords)
	files<-grep("\\.psmc$",dir(),value=TRUE)
	dev.new(height=height,width=width)
	plot(1,1,
		ylim=ylim,xlim=xlim,
		log="x",type="n",
		main=label[1],ylab="Ne",xlab="Years ago")
	for(i in 1:n)
	{
		subfiles<-grep(paste("^",keywords[i],sep=""),files,value=TRUE)
		n.sub<-length(subfiles)
		for(i.sub in 1:n.sub)
		{
			lines(psmc.result(subfiles[i.sub],i.iteration,mu,s,g),
				type="l",col=col[i],lwd=2)
		}
	}
	legend("topright",legend=legend.names,col=col,lty=1,lwd=2)
	savePlot(filename=paste(label[1],save.as,sep="."),type=save.as)
	dev.off()
}


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
	if (any(duplicated(nms))) {
		warning("Some names are duplicated, creating unique IDs but might want to check.")
		nms <- make.unique(as.character(nms))
	}
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
#
# Revise code in RERconverge to allow for a progress update and parallel processing

# simBinPhenoCC(trees, mastertree, root_sp, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=T)

# simBinPhenoCC(trees=tree_rep[[1]], cnee_trees$masterTree, "momMom", islandFg, pathvec=pathvec, sisters_list=sisters_island)

generatePermulatedBinPhen <- function (tree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode = "cc", ncores) {
	require(pbmcapply)
	require(pbapply)
	# require(ape)
    if (permmode == "cc") {
        tree_rep = lapply(1:numperms, RERconverge:::rep_tree, tree = trees)
        if (!is.null(ncores)) {
    		permulated.binphens <- pbmclapply(tree_rep, mc.cores=ncores, simBinPhenoCC, 
            	mastertree = trees$masterTree, root_sp = root_sp, 
            	fg_vec = islandFg, sisters_list = sisters_island, pathvec = pathvec, 
            	plotTreeBool = F)
        } else {
        	permulated.binphens <- pblapply(tree_rep, simBinPhenoCC, 
	            mastertree = trees$masterTree, root_sp = root_sp, 
	            fg_vec = fg_vec, sisters_list = sisters_list, pathvec = pathvec, 
	            plotTreeBool = F)
        }
    }
    else if (permmode == "ssm") {
        tree_rep = lapply(1:numperms, rep_tree, tree = tree)
        permulated.binphens = lapply(tree_rep, simBinPhenoSSM, 
            trees = trees, root_sp = root_sp, fg_vec = fg_vec, 
            sisters_list = sisters_list, pathvec = pathvec)
    }
    else {
        stop("Invalid binary permulation mode.")
    }
    output.list <- list()
    output.list[[1]] <- permulated.binphens
    return(output.list)
}
environment(generatePermulatedBinPhen) <- asNamespace('RERconverge')
assignInNamespace("generatePermulatedBinPhen", generatePermulatedBinPhen, ns = "RERconverge")


# numperms=500
# fg_vec=islandFg
# sisters_list=NULL
# root_sp="momMom"
# RERmat=kingRERw_cnee
# trees=cnee_trees
# mastertree=masterTree_cnee
# permmode="cc"
# ncores=4

# getDepthOrder <- RERconverge:::getDepthOrder

getPermsBinary <- function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees, mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, calculateenrich = F, annotlist = NULL, ncores=NULL) {
    require(pbmcapply)
    require(pbapply)
    pathvec = foreground2Paths(islandFg, cnee_trees, clade = "all", plotTree = F)
    col_labels = colnames(trees$paths)
    names(pathvec) = col_labels
    if (permmode == "cc") {
        print("Running CC permulation")
        print("Generating permulated trees")
        # fg_tree_depth_order not found
        permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 
            numperms, cnee_trees, root_sp, islandFg, sisters_list=NULL, pathvec, 
            permmode = "cc", ncores = ncores)
        permulated.fg = mapply(RERconverge:::getForegroundsFromBinaryTree, 
            permulated.binphens[[1]])
        permulated.fg.list = as.list(data.frame(permulated.fg))
        phenvec.table = mapply(foreground2Paths, permulated.fg.list, 
            MoreArgs = list(treesObj = trees, clade = "all"))
        phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[, 
            i])
        print("Calculating correlations")
        if (!is.null(ncores)) {
            corMatList = pbmclapply(phenvec.list, mc.cores=ncores, correlateWithBinaryPhenotype,
                RERmat = RERmat)
        } else {
            corMatList = pblapply(phenvec.list, correlateWithBinaryPhenotype,
                RERmat = RERmat)
        }
        permPvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
        rownames(permPvals) = rownames(RERmat)
        permRhovals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
        rownames(permRhovals) = rownames(RERmat)
        permStatvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
        rownames(permStatvals) = rownames(RERmat)
        for (i in 1:length(corMatList)) {
            permPvals[, i] = corMatList[[i]]$P
            permRhovals[, i] = corMatList[[i]]$Rho
            permStatvals[, i] = sign(corMatList[[i]]$Rho) * -log10(corMatList[[i]]$P)
        }
    }
    else if (permmode == "ssm") {
        print("Running SSM permulation")
        if (is.null(trees_list)) {
            trees_list = trees$trees
        }
        RERmat = RERmat[match(names(trees_list), rownames(RERmat)), 
            ]
        print("Generating permulated trees")
        permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list, 
            numperms, trees, root_sp, fg_vec, sisters_list, pathvec)
        df.list = lapply(trees_list, getSpeciesMembershipStats, 
            masterTree = masterTree, foregrounds = fg_vec)
        df.converted = data.frame(matrix(unlist(df.list), nrow = length(df.list), 
            byrow = T), stringsAsFactors = FALSE)
        attr = attributes(df.list[[1]])
        col_names = attr$names
        attr2 = attributes(df.list)
        row_names = attr2$names
        colnames(df.converted) = col_names
        rownames(df.converted) = row_names
        df.converted$num.fg = as.integer(df.converted$num.fg)
        df.converted$num.spec = as.integer(df.converted$num.spec)
        spec.members = df.converted$spec.members
        grouped.trees = groupTrees(spec.members)
        ind.unique.trees = grouped.trees$ind.unique.trees
        ind.unique.trees = unlist(ind.unique.trees)
        ind.tree.groups = grouped.trees$ind.tree.groups
        unique.trees = trees_list[ind.unique.trees]
        unique.map.list = mapply(matchAllNodesClades, unique.trees, 
            MoreArgs = list(treesObj = trees))
        unique.permulated.binphens = permulated.binphens[ind.unique.trees]
        unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens, 
            unique.map.list, trees)
        permulated.paths = vector("list", length = length(trees_list))
        for (j in 1:length(permulated.paths)) {
            permulated.paths[[j]] = vector("list", length = numperms)
        }
        for (i in 1:length(unique.permulated.paths)) {
            ind.unique.tree = ind.unique.trees[i]
            ind.tree.group = ind.tree.groups[[i]]
            unique.path = unique.permulated.paths[[i]]
            for (k in 1:length(ind.tree.group)) {
                permulated.paths[[ind.tree.group[k]]] = unique.path
            }
        }
        attributes(permulated.paths)$names = row_names
        print("Calculating correlations")
        RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i, 
            ])
        corMatList = mapply(calculateCorPermuted, permulated.paths, 
            RERmat.list)
        permPvals = extractCorResults(corMatList, numperms, mode = "P")
        rownames(permPvals) = names(trees_list)
        permRhovals = extractCorResults(corMatList, numperms, 
            mode = "Rho")
        rownames(permRhovals) = names(trees_list)
        permStatvals = sign(permRhovals) * -log10(permPvals)
        rownames(permStatvals) = names(trees_list)
    }
    else {
        stop("Invalid binary permulation mode.")
    }
    if (calculateenrich) {
        realFgtree = foreground2TreeClades(fg_vec, sisters_list, 
            trees, plotTree = F)
        realpaths = tree2PathsClades(realFgtree, trees)
        realresults = getAllCor(RERmat, realpaths, method = method, 
            min.pos = min.pos)
        realstat = sign(realresults$Rho) * -log10(realresults$P)
        names(realstat) = rownames(RERmat)
        realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, 
            outputGeneVals = F)
        groups = length(realenrich)
        c = 1
        while (c <= groups) {
            current = realenrich[[c]]
            realenrich[[c]] = current[order(rownames(current)), 
                ]
            c = c + 1
        }
        permenrichP = vector("list", length(realenrich))
        permenrichStat = vector("list", length(realenrich))
        c = 1
        while (c <= length(realenrich)) {
            newdf = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[c]])))
            rownames(newdf) = rownames(realenrich[[c]])
            permenrichP[[c]] = newdf
            permenrichStat[[c]] = newdf
            c = c + 1
        }
        counter = 1
        while (counter <= numperms) {
            stat = permStatvals[, counter]
            names(stat) = rownames(RERmat)
            enrich = fastwilcoxGMTall(na.omit(stat), annotlist, 
                outputGeneVals = F)
            groups = length(enrich)
            c = 1
            while (c <= groups) {
                current = enrich[[c]]
                enrich[[c]] = current[order(rownames(current)), 
                  ]
                enrich[[c]] = enrich[[c]][match(rownames(permenrichP[[c]]), 
                  rownames(enrich[[c]])), ]
                permenrichP[[c]][, counter] = enrich[[c]]$pval
                permenrichStat[[c]][, counter] = enrich[[c]]$stat
                c = c + 1
            }
            counter = counter + 1
        }
    }
    if (calculateenrich) {
        data = vector("list", 5)
        data[[1]] = permPvals
        data[[2]] = permRhovals
        data[[3]] = permStatvals
        data[[4]] = permenrichP
        data[[5]] = permenrichStat
        names(data) = c("corP", "corRho", "corStat", "enrichP", 
            "enrichStat")
    }
    else {
        data = vector("list", 3)
        data[[1]] = permPvals
        data[[2]] = permRhovals
        data[[3]] = permStatvals
        names(data) = c("corP", "corRho", "corStat")
    }
    data
}
environment(getPermsBinary) <- asNamespace('RERconverge')
assignInNamespace("getPermsBinary", getPermsBinary, ns = "RERconverge")

# files <- paste0("MAnull_results/", list.files("MAnull_results", pattern="_"), "/mlc")

# files = list of mlc files to parse
parse_codeml <- function(files, ncores=1) {
	require(pbmcapply)
	ll <- pbmclapply(files, mc.cores=24, function(x) {
		try(raw <- readLines(x))
		res <- stringr::str_match(raw[grep("lnL", raw)], "lnL.*?np.*?(\\d+)\\)\\:\\s+([\\.\\-\\d]+)")[,2:3]
		# Remove genes that had NaN errors in PAML model fitting
		tryCatch(setNames(res, c('np', 'lnL')), error=function(err) c(NA, NA))
	})
	# Create sub data frame
	res <- cbind(file = basename(dirname(files)), do.call(rbind, ll))
	res <- as.data.frame(res)
	res$lnL <- as.numeric(as.character(res$lnL))
	res$np <- as.numeric(as.character(res$np))
	return(res)
}


permpvalcor_fast <- function(realcor, permvals) {
	denom <- apply(permvals$corRho, 1, function(x) sum(!is.na(x)))
	numer <- t(apply(permvals$corRho, 1, function(x) abs(x) ))
	numer2 <- pbsapply(1:nrow(numer), function(x) {
		sum(numer[x, ] > abs(realcor$Rho[x]))
	})
	numer2 / denom # P values
}

