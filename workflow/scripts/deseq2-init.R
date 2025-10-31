# log <- file(snakemake@log[[1]], open="wt")
# sink(log)
# sink(log, type="message")
library(readr)
library("DESeq2")
parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  # setup parallelization
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}
if (!exists("snakemake")) {
  library(magrittr)
  the_yaml <- yaml::read_yaml("./configs/VascAge_config.yaml")
  BASE_DIR <- the_yaml$dirs$BASE_ANALYSIS_DIR
  cts <- file.path(BASE_DIR, "counts/all.tsv") %>%
    read.table(header = T, row.names = "gene", check.names = F)
  coldata <- readr::read_delim("./data/Vasc_age2020/Vascage_samples.tsv",
    header = TRUE,
    row.names = "gene",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  # coldata %>% dplyr::filter(cell_type == "tumor") ->coldata
  # cts <- cts[,rownames(small_coldata)]
  snakemake_conf <- yaml::read_yaml("./configs/VascAge_config.yaml")
  all_conditions <- names(snakemake_conf$diffexp$contrasts)
}
all_conditions <- names(snakemake@config$diffexp$contrasts)


# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]],
  header = TRUE, sep = "\t",
  row.names = 1,
  check.names = FALSE
)


#sample_mat <- readr::read_tsv(samples)

coldata <- read_tsv(snakemake@params[["samples"]],
  col_names = TRUE,
) 
stopifnot("sample" %in% colnames(coldata))
coldata <- coldata |> tibble::column_to_rownames(var = "sample")

## Reorder coldata rows to match cts col order (beacause deseq things)
if (!all(colnames(cts) == rownames(coldata))) {
  sample_ids <- colnames(cts)
  coldata <- coldata[sample_ids, ]
}
# Remove NAs from Data (not supported)
if (any(is.na(coldata[, c(all_conditions)]))) {
  na_index <- apply(coldata, 1, function(x) {
    any(is.na(x))
  })
  coldata <- coldata[!na_index, ]
  cts <- cts[, rownames(coldata)]
}
all_conditions <- names(snakemake@config$diffexp$contrasts)
for (x in all_conditions) {
  coldata[, x] <- as.factor(coldata[, x])
  coldata[, x] <- relevel(coldata[, x], ref = snakemake@config$diffexp$contrasts[[x]][[1]][[2]])
}

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  #  design = as.formula("~condition"))
  design = as.formula(snakemake@params[["model"]]),
)




# remove genes defined in config
excluded_genes <- snakemake@config[["excluded_genes"]]
# print(excluded_genes)
if (length(excluded_genes) > 0) {
  ex_index <- purrr::map_int(excluded_genes,
    grep, rownames(dds),
    value = F
  )
  if (length(ex_index) > 0) {
    dds <- dds[-ex_index, ]
  }
}
# remove uninformative rows
dds <- dds[rowSums(counts(dds)) > ncol(dds) / 2, ]

# normalization and preprocessing
dds <- DESeq(dds,
  parallel = parallel
)
if (!is.null(snakemake@config[["cpm_filter"]])) {
  cpm_threshold <- snakemake@config[["cpm_filter"]]
} else {
  cpm_threshold <- 10
}

cpm_filter <- apply(edgeR::cpm(counts(dds, normalized = T)), 1, function(x) {
  if(length(which(x > cpm_threshold)) > 0.2 * length(x)) {
    val <- 1
  } else {
    val <- 0 
  }
  as.logical(val)
})
dds <- dds[cpm_filter,]
saveRDS(dds, file = snakemake@output[[1]])
