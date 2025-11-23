library(sva)
library(readr)
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)  # for mouse
library(org.Hs.eg.db)  # for human

#' ENSEMBL IDs to Gene Symbol using AnnotationDbi
#'
#' @param ens_id_vector Vector of Ensembl IDs. Version numbers removed via regex
#' @param organism Name of organism: "mouse" or "human"
#' @return Named vector of gene symbols
#' @export
#'
#' @examples
ensembl_to_symbol <- function(ens_id_vector, organism = "mouse") {
    # Remove version numbers from ENSEMBL IDs
    ens_short <- stringr::str_extract(ens_id_vector, pattern = "^ENS[A-Z0-9]*")
    
    # Select appropriate database
    if (organism == "mouse") {
        db <- org.Mm.eg.db
    } else if (organism == "human") {
        db <- org.Hs.eg.db
    } else {
        stop("Organism must be 'mouse' or 'human'")
    }
    
    # Initialize output vector with ENSEMBL IDs as default
    symbol_vec <- setNames(ens_short, nm = ens_short)
    
    # Map ENSEMBL IDs to gene symbols
    mapped <- AnnotationDbi::mapIds(
        db,
        keys = ens_short,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
    )
    
    # Replace successfully mapped IDs, keep ENSEMBL ID for unmapped
    mapped[is.na(mapped)] <- names(mapped)[is.na(mapped)]
    symbol_vec[names(mapped)] <- mapped
    
    # Report mapping statistics
    n_mapped <- sum(!is.na(mapped) & mapped != names(mapped))
    message(sprintf("Mapped %d/%d genes to symbols (%.1f%%)", 
                    n_mapped, length(ens_short), 
                    100 * n_mapped / length(ens_short)))
    
    return(symbol_vec)
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