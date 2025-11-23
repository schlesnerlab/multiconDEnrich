library(glmmSeq)
library(ggplot2)
library(tibble)
library(dplyr)
run_gsea_glmmseq <- function(fixed_eff,
                             coef,
                             glmmseq_obj,
                             cate = "H",
                             subcat = NULL,
                             custom_geneset = NULL) {
  # Debug: check rownames
  cat(sprintf("Processing %s (coef: %s)\n", fixed_eff, coef))
  cat(sprintf("Total genes: %d\n", nrow(glmmseq_obj@stats$coef)))
  cat(sprintf("Genes with NA rownames: %d\n", sum(is.na(rownames(glmmseq_obj@stats$coef)))))
  cat(sprintf("Genes with empty rownames: %d\n", sum(rownames(glmmseq_obj@stats$coef) == "")))
  
  out_dat <-
    gsea_input <- tibble(gene_symbol = rownames(glmmseq_obj@stats$coef),
                         score = glmmseq_obj@stats$coef[, coef]) %>%
    # Remove NA and empty gene symbols
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "", is.finite(score)) %>%
    dplyr::arrange(desc(score)) %>% # remove duplicate gene symobls
    distinct(gene_symbol, .keep_all = TRUE)
  
  cat(sprintf("Genes after filtering for coefficient ranking: %d\n", nrow(gsea_input)))

  gsea_pval <-
    tibble(
      gene_symbol = rownames(glmmseq_obj@stats$coef),
      score = glmmseq_obj@stats$coef[, coef] *
        -log10(glmmseq_obj@stats$pvals[, fixed_eff])
    ) %>%
    # Remove NA and empty gene symbols
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "", is.finite(score)) %>%
    dplyr::arrange(desc(score)) |>
    distinct(gene_symbol, .keep_all = TRUE)
  
  cat(sprintf("Genes after filtering for pval ranking: %d\n", nrow(gsea_pval)))
  gsea_input$score[gsea_input$score == Inf] <-
    max(gsea_input$score[gsea_input$score < Inf], na.rm = TRUE)
  gsea_pval$score[gsea_pval$score == Inf] <-
    max(gsea_pval$score[gsea_pval$score < Inf], na.rm = TRUE)
  
  # Replace negative infinities with the minimum non-infinite value
  
  
  gsea_input$score[gsea_input$score == -Inf] <-
    min(gsea_input$score[gsea_input$score > -Inf], na.rm = TRUE)
  gsea_pval$score[gsea_pval$score == -Inf] <-
    min(gsea_pval$score[gsea_pval$score > -Inf], na.rm = TRUE)
  gsea_input$score[gsea_input$score == Inf] <-
    max(gsea_input$score[gsea_input$score < Inf], na.rm = TRUE)
  gsea_pval$score[gsea_pval$score == Inf] <-
    max(gsea_pval$score[gsea_pval$score < Inf], na.rm = TRUE)
  
  #gsea_pval <- gsea_pval[is.finite(gsea_pval$score), ]
  #gsea_input <- gsea_input[is.finite(gsea_input$score), ]
  
  # Final check before GSEA
  cat(sprintf("Final gene counts - coef: %d, pval: %d\n", 
              nrow(gsea_input), nrow(gsea_pval)))
  
  # Check for duplicates
  if (any(duplicated(gsea_input$gene_symbol))) {
    cat("Warning: duplicated gene symbols found in gsea_input\n")
  }
  if (any(duplicated(gsea_pval$gene_symbol))) {
    cat("Warning: duplicated gene symbols found in gsea_pval\n")
  }
  
  tryCatch({
    gsea_res <-
      RNAscripts::run_msig_enricher(
        gset_list = list(gsea_input, gsea_pval),
        category = cate,
        subcategory = subcat,
        pvalueCutoff = 0.05,  
        custom_geneset = custom_geneset,
      )
    
    names(gsea_res) <- c("coef", "pvalue")
    gsea_res
  }, error = function(e) {
    cat(sprintf("Error in run_msig_enricher: %s\n", e$message))
    cat("First few rows of gsea_input:\n")
    print(head(gsea_input))
    cat("\nFirst few rows of gsea_pval:\n")
    print(head(gsea_pval))
    stop(e)
  })
}

run_gsea_across_genesets <- function(glmmseq_obj, gset_config) {
  gsea_res <- list()
  for (gset in names(gset_config)) {
    print( gset_config[[gset]])
      if(!is.null(gset_config[[gset]]$database)) {
        custom_set <- gset_config[[gset]]$database
        if (custom_set == "custom_senescence"){ 
          custom_set <- RNAscripts::senes
        } else if (gset_config[[gset]]$database == "MitoCarta") {
          custom_set <- RNAscripts::MitoPathways
        } else {
          custom_set <- NULL
        }
      } else {
        custom_set <- NULL
      }
    gsea_res[[gset]] <- run_gsea_glmmseq(
      fixed_eff = var,
      coef = coef,
      cate = gset_config[[gset]]$category,
      subcat = gset_config[[gset]]$subcategory,
      glmmseq_obj = glmmseq_obj,
      custom_geneset = custom_set
    )
  }
  gsea_res
}

if (exists("snakemake")) {
    glmmseq_file <- snakemake@input[["glmmseq_obj"]]
    enrichment_obj <- snakemake@output[["enrichment_obj"]]
    enrichments <- snakemake@params[["enrichments"]]
    coef <- snakemake@params[["coef_name"]]
    var <- snakemake@wildcards[["var"]]
    gset_config <- snakemake@config[["glmmseq"]][["enrichments"]]
} else {
  glmmseq_file <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/glmmseq/glmmseq/glmmseq_obj.rds.gz"
  coef <- "ageaged"
  var <- "age"
  gset_config <- list(
  Hallmark_gsea = list(
    category = "H",
    subcategory = NULL,
    use_gsea = TRUE,
    database = "MSigDB",
    id_class = "SYMBOL"
  ),
  Reactome_GSEA = list(
    category= NULL,
    subcategory= NULL,
    use_gsea= TRUE,
    database= "custom_senescence"
  )
)
}
write(c(coef,var), file = stderr())
glmmseq_list <- readRDS(glmmseq_file)
glmmseq_list$norm_counts <- glmmSeq::glmmQvals(glmmseq_list$norm_counts)

enrich_out <- run_gsea_across_genesets(glmmseq_obj = glmmseq_list$norm_counts, 
                                       gset_config = gset_config)

saveRDS(enrich_out, file = enrichment_obj)
