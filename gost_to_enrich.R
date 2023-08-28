gost_to_enrich <- function(gp) {
	require(DOSE)
	require(enrichplot)
	# modify the g:Profiler data frame
	gp_mod = gp$result[,c("query", "source", "term_id", "term_name", "p_value", "query_size", "intersection_size", "term_size", "effective_domain_size", "intersection")]
	gp_mod$GeneRatio = paste0(gp_mod$intersection_size,  "/", gp_mod$query_size)
	gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
	names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust",
	                  "query_size", "Count", "term_size", "effective_domain_size",
	                  "geneID", "GeneRatio", "BgRatio")
	gp_mod$geneID = gsub(",", "/", gp_mod$geneID)
	row.names(gp_mod) = gp_mod$ID
	# Convert to enrichResult object for plotting
	gp_mod_enrich = new("enrichResult", result = gp_mod)
	gp_mod_enrich	
}
