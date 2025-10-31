library(clusterProfiler)
library(magrittr)
library(furrr)
# library(tidyverse)
library(RNAscripts)

## Snakemake header
options <- furrr_options(seed = 6827)
if (exists("snakemake")) {
  diffexp_tb_path <- snakemake@input[["table"]]
  fpkm_path <- snakemake@input[["fpkm_path"]]
  contrast_groups <- snakemake@params[["contrast"]]
  gsea_use_stat <- snakemake@params[["gsea_use_stat"]]
  pvalue_threshold <- snakemake@config[["diffexp"]][["pval_threshold"]]
  LFC_threshold <- snakemake@config[["diffexp"]][["LFC_threshold"]]
  out_file <- snakemake@output[["gsea_result"]]
  organism <- snakemake@config[["organism"]]
  gsea_config <- snakemake@config[["gsea"]]
  plan(strategy = multicore, workers = snakemake@threads)
} else {
  conf <- yaml::read_yaml("./configs/tec_aging.yaml")
  BASE_ANALYSIS_DIR <- file.path(conf$dirs$BASE_ANALYSIS_DIR)

  gsea_config <- conf$gsea
  cond_id <- names(conf$diffexp$contrasts)[1]
  comp_id <- names(conf$diffexp$contrasts[[cond_id]])[1]

  contrast_groups <- conf$diffexp$contrasts[[cond_id]][[comp_id]]
  diffexp_tb_path <- file.path(BASE_ANALYSIS_DIR, glue::glue("results/diffexp/{cond_id}/{comp_id}.diffexp.tsv"))
  fpkm_path <- file.path(BASE_ANALYSIS_DIR, "fpkm/all.tsv")
  pvalue_threshold <- 0.05
  LFC_threshold <- 0.5
  organism <- "Mus musculus"
  gsea_use_stat <- TRUE
  plan(strategy = sequential)
}
org_db <- RNAscripts::get_org_db(organism)

diffxp_tb <- readr::read_tsv(diffexp_tb_path, col_names = c(
  "gene_id", "baseMean", "logFoldChange", "lfcSE", "stat", "pvalue",
  "padj"
), skip = 1)
fpkm <- readr::read_tsv(fpkm_path)

# Read Files
filer <- fpkm %>%
  dplyr::filter(gene %in% diffxp_tb$gene_id)

joined_df <- RNAscripts::join_tables(diffxp_tb, filer)

joined_df <- joined_df %>%
  dplyr::mutate(overexpressed_in = ifelse(logFoldChange > 0, contrast_groups[1], contrast_groups[2]))

if (gsea_use_stat) {
  joined_df <- joined_df %>% dplyr::mutate(gsea_stat = joined_df$stat)
} else {
  joined_df <- joined_df %>%
    dplyr::mutate(gsea_stat = -log10(pvalue) * logFoldChange)
}

joined_df <- joined_df %>% dplyr::arrange(desc(gsea_stat))
gene_list <- joined_df %>%
  dplyr::select(c(gname, gsea_stat))
ensemblgene_list <- joined_df %>%
  dplyr::select(c(gene, gsea_stat))
de_genes <- joined_df %>%
  dplyr::filter(padj < pvalue_threshold & abs(gsea_stat) > LFC_threshold) %>%
  dplyr::select(c(gene, stat, gsea_stat))

t_table <- RNAscripts::table_to_list(joined_df, 'gname', 'gene')

# t_table <- RNAscripts::get_entrezgenes_from_ensembl(filer %>% dplyr::pull(gene), input_type = 'ENSEMBL', org_db =
# org_db ) %>% RNAscripts::table_to_list(., 'ENTREZID', 'ENSEMBL')

enrich_data <- furrr::future_map(names(gsea_config), RNAscripts::run_gsea_query,
  gsea_genes = ensemblgene_list,
  de_genes = de_genes,
  gset_config = gsea_config,
  species = organism,
  org_db = org_db,
  t_table = t_table
)
names(enrich_data) <- names(gsea_config)

saveRDS(object = enrich_data, file = out_file)
