library(sva)
library(readr)
library(DESeq2)
library(biomaRt)
biomartCacheClear()
stable_get_bm <- function(species, host_url = "https://nov2020.archive.ensembl.org") {
    mart <- "www"
    rounds <- 0
    while (class(mart)[[1]] != "Mart") {
        mart <- tryCatch(
            {
                # done here, because error function does not
                # modify outer scope variables, I tried
                if (mart == "www") rounds <- rounds + 1
                # equivalent to useMart, but you can choose
                # the mirror instead of specifying a host
                biomaRt::useEnsembl(
                    biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset = glue::glue("{species}_gene_ensembl"),
                    mirror = mart,
                    host = host_url
                )
            },
            error = function(e) {
                # change or make configurable if you want more or
                # less rounds of tries of all the mirrors
                if (rounds >= 3) {
                    stop()
                }
                # hop to next mirror
                mart <- switch(mart,
                               useast = "uswest",
                               uswest = "asia",
                               asia = "www",
                               www = {
                                   # wait before starting another round through the mirrors,
                                   # hoping that intermittent problems disappear
                                   Sys.sleep(30)
                                   "useast"
                               },
                               host = host_url
                )
            }
        )
    }
    mart
}


#' ENSEMBL IDs to ENSEMBL Gene symbol
#'
#' @param ens_id_vector Vector of Ensembl IDs. Version numbers removed via regex
#' @param organism_name Name of organism as used in ensembl "mmusculus" 
#' @return
#' @export
#'
#' @examples
ensembl_to_symbol <- function(ens_id_vector, organism_name = "mmusculus") {
    ens_short <- stringr::str_extract(ens_id_vector,
                                      pattern = "^ENS[A-Z0-9]*")
    mart <- biomaRt::useEnsembl(  biomart = "ENSEMBL_MART_ENSEMBL",
                                  dataset = glue::glue("{organism_name}_gene_ensembl"),
                                  host = "https://nov2020.archive.ensembl.org")
    g2g <- biomaRt::getBM(
        attributes = c( "ensembl_gene_id",
                        "external_gene_name"),
        filters = "ensembl_gene_id",
        values = ens_short,
        mart = mart,  
    )
    symbol_vec <- setNames(ens_short, nm = ens_short)
    symbol_vec[g2g$ensembl_gene_id] <- g2g$external_gene_name
    
    symbol_vec
}

# snakemake Boilerplate
if (exists("snakemake", inherits = FALSE)) {
    counts <- snakemake@input[["counts"]]
    samples <- snakemake@params[["samples"]]
    batch_variable <- snakemake@params[["batch_variable"]]
    batch_correct_model <- snakemake@params[["batch_correct_model"]] |> as.formula()
    the_formula <- snakemake@config[["glmmseq"]][["formula"]] |> as.formula()
    batch_corrected_counts <- snakemake@output[["batch_corrected_counts"]]
    uncorrected_counts <- snakemake@output[["uncorrected_counts"]]
    reference_groups <- snakemake@config[["glmmseq"]][["reference_group"]]
    parallel <- TRUE
} else {
    counts <- "/desktop-home/heyer/projects/Vascular_Aging/Integrative_Analysis/join_RNA/count.tsv"
    samples <- "/desktop-home/heyer/projects/Vascular_Aging/Integrative_Analysis/join_RNA/full_metadata.tsv"
    batch_correct_model <- ~as.factor(age) + as.factor(Aplnr_KO) + as.factor(Apln_treatment) + as.factor(EC_status)
    batch_variable <- "experiment"  # Define a default batch variable if not using snakemake
    batch_corrected_counts <- "corrected_counts.rds"  # Define a default output path if not using snakemake
    uncorrected_counts <- "/desktop-home/heyer/projects/Vascular_Aging/Integrative_Analysis/join_RNA/dds_test.rds.gz"
    reference_groups <- list("age" = "young", "EC_status" = "healthy")
    the_formula <- as.formula("~ age * EC_status + Apln_treatment:age + 
                        Aplnr_KO + EC_status:Aplnr_KO + Apln_treatment +  (1 | condition) + experiment")
}

# Read counts matrix
count_mat <-read.table(counts,  header = TRUE,
                       row.names = 1,
                       check.names = FALSE
)

sample_mat <- readr::read_tsv(samples)

# Ensure all count_mat columns are present in sample_mat
stopifnot(all(colnames(count_mat) %in% sample_mat$sample))
# Reorder sample mat to match count mat
sample_mat <- sample_mat[sample_mat$sample %in% colnames(count_mat),]
keep <- rowSums(count_mat) >= nrow(count_mat) * 0.3
count_mat <- count_mat[keep,]
sample_mat[,all.vars(the_formula)] <- purrr::modify(sample_mat[,all.vars(the_formula)], as.factor )

for (var_name in names(reference_groups)) {
    sample_mat[[var_name]] <- relevel(sample_mat[[var_name]], ref = reference_groups[[var_name]])
}
# Check if all variables in model are present in metadata
stopifnot(all(all.vars(batch_correct_model) %in% colnames(sample_mat)))

m <- model.matrix(batch_correct_model, data = sample_mat)

uncorrected_dds <- DESeqDataSetFromMatrix(count_mat, colData = sample_mat, design = batch_correct_model)
rownames(uncorrected_dds) <-  ensembl_to_symbol(rownames(uncorrected_dds))
uncorrected_dds <- DESeq(uncorrected_dds, parallel = F)

saveRDS(uncorrected_dds, file = uncorrected_counts )
# Batch correction
corrected_counts <- sva::ComBat_seq(as.matrix(count_mat), batch = sample_mat |> dplyr::pull(!!batch_variable), covar_mod = m)

# Initialize DESeq2 object
corrected_dds <- DESeqDataSetFromMatrix(corrected_counts, colData = sample_mat, design = batch_correct_model)
#keep <- rowSums(counts(corrected_dds)) >= 30

corrected_dds <- DESeq(corrected_dds)
rownames(corrected_dds) <-  ensembl_to_symbol(rownames(corrected_dds))


# Write corrected counts matrix
saveRDS(corrected_dds, file = batch_corrected_counts)