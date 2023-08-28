#' Function to get consensus sequences from BAM alignments
#' 
#' @param ref path to reference genome
#' @param bam path to BAM alignment file
#' @param call whether to call variants
#' @param index whether to index the VCF variant file output
#' @param cons whether to generate consensus sequence
#' @param suffix suffix of genome file (e.g., ".con.fa")
#' 
#' @export
#' 
getConsensus <- function(ref, bam, vcfpath=NULL, outpath=NULL, call=TRUE, index=TRUE, cons=TRUE, suffix=".con.fa", force=FALSE, filter="QUAL>30 && DP>5", onlyvar=TRUE) {

	nm <- stringr::str_extract(basename(bam), ".*?(?=\\.bam)")
	
	if (is.null(outpath)) {
		outpath <- paste0("genomes/", substr(nm, 1, 6))
	}
	if (!dir.exists(outpath)) {
		dir.create(outpath)
	}
	# call variants (filter out indels with -I - Taylor Hains suggested doing this to maintain reference coordinates)
	if (call) {
		if (!force & file.exists(paste0(outpath, "/", nm, ".calls.vcf.gz"))) {
			stop("VCF file already exists, use `force=TRUE` to overwrite")
		} else {
			if (onlyvar) {
				system(paste0("bcftools mpileup -I -Ou -f ", ref, " ", bam, " | bcftools call -mv -Oz -o ", vcfpath, "/", nm, ".calls.vcf.gz"))
			} else {
				system(paste0("bcftools mpileup -I -Ou -f ", ref, " ", bam, " | bcftools call -m -Oz -o ", vcfpath, "/", nm, ".calls.vcf.gz"))	
			}
		}		
	}
	# index
	if (index) {
		system(paste0("bcftools index ", vcfpath, "/", nm, ".calls.vcf.gz"))
	}
	# normalize indels
	# system(paste0("bcftools norm -f ", ref, " ", nm, ".calls.vcf.gz -Ob -o ", nm, ".calls.norm.bcf"))
	# filter adjacent indels within 5bp
	# system(paste0("bcftools filter --IndelGap 5 calls.norm.bcf -Ob -o calls.norm.flt-indels.bcf")
	if (cons) {
		if (!force & file.exists(paste0(outpath, "/", nm, suffix))) {
			stop("Consensus FASTA already exists, try a different suffix")
		} else {
			system(paste0("cat ", ref, " | bcftools consensus -i '", filter, "' ", vcfpath, "/", nm, ".calls.vcf.gz > ", outpath, "/", nm, suffix))
		}		
	}
	# This assumes we have already made the calls, normalized indels and filtered. There is another page which goes deeper and is devoted just to this, but in brief, the variant calling command in its simplest form is:
}
