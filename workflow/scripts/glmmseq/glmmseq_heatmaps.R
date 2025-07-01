
library(ComplexHeatmap)
library(circlize)
library(glmmSeq)
library(seriation)
library(DESeq2)
register_DendSer()
if (exists("snakemake")) {
    glmmseq_obj <- snakemake@input[["glmmseq_obj"]]
    corrected_counts <- snakemake@input[["batch_corrected_counts"]]
    png_file <- snakemake@output[["png_file"]]
    var <- snakemake@wildcards[["var"]]
    coef_id <- snakemake@params[["coef"]]
    col_data_to_plot <- snakemake@config[["pca"]][["labels"]] |> unlist()
    group_colors <- snakemake@config[["group_colors"]]
    print(group_colors)
} else {
    glmmseq_obj <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/glmmseq/glmmseq/glmmseq_obj.rds.gz"
    corrected_counts <- "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/glmmseq/counts/batch_corrected_counts.rds"
    coef_id <- "cre_statusTerf1KO"
    var <- "cre_status"
    col_data_to_plot <- c("EC_status", "cre_status")
    group_colors <- NULL
}

if (length(group_colors) == 0) {
  # Dynamically generate annotation_colors for all groups in col_data_to_plot
  annotation_colors <- list()
  color_palette <- function(n) {
    # Use RColorBrewer if available, else fallback to rainbow
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      cols <- RColorBrewer::brewer.pal(min(n, 8), "Set2")
      if (n > 8) cols <- c(cols, RColorBrewer::brewer.pal(n - 8, "Set3"))
      return(cols[1:n])
    } else {
      return(rainbow(n))
    }
  }
  for (group in col_data_to_plot) {
    group_levels <- NULL
    # Try to get levels from colData if available, else fallback to unique values
    if (exists("glmmseq_obj") && file.exists(glmmseq_obj)) {
      coldata <- tryCatch({
        readRDS(glmmseq_obj)$norm_counts@metadata
      }, error = function(e) NULL)
      if (!is.null(coldata) && group %in% colnames(coldata)) {
        group_levels <- unique(as.character(coldata[[group]]))
      }
    }
    if (is.null(group_levels)) {
      group_levels <- c("level1", "level2") # fallback if unknown
    }
    annotation_colors[[group]] <- setNames(color_palette(length(group_levels)), group_levels)
  }
} else {
  annotation_colors <- group_colors
}

plot_heatmap <- function(coef,vst_obj = vst_dds, glmm_obj =  glmmseq_norm_counts,
                        z_score = TRUE, use_vst = FALSE, 
                        ann_col = annotation_colors,
                        coldata_to_plot = c("age", "EC_status","Apln_treatment", "Aplnr_KO", "experiment"),
                        n_genes = 30, qval_filt = 0.01, coef_filt = 0.8, meanExp_cutoff = 7,
                        row_km = 3, col_km = NULL) {
    # Get qvals and coefs for one coeff from glmm_obj
    qvals <- glmm_obj@stats$qvals[,coef[1]]
    coef_values <- glmm_obj@stats$coef[,coef[2]]
    meanexp <- glmm_obj@stats$res[,"meanExp"]
    
    stopifnot(all.equal(names(coef_values) , names(qvals)))
    # get vst normalized expression data 
    vst_data <- assay(vst_obj)

    # retain only genes with qvals < 0.05 and ceef_value > abs(0.5)o
    selected_genes <- which(qvals < 0.05 )
    if (length(selected_genes) > 200) {
    selected_genes <- which(qvals < qval_filt & abs(coef_values) > coef_filt &
                              meanexp > meanExp_cutoff)
    }
    length(selected_genes)
    
    if (use_vst) {
      plot_data <- vst_data#[names(selected_genes),]
      #plot_data <- t(scale(t(plot_data)))
      coldata <- colData(vst_obj)
    } else {
      plot_data <- glmm_obj@countdata[selected_genes,]
      #plot_data <- t(scale(t(plot_data)))

      coldata <- glmm_obj@metadata[,coldata_to_plot]
    }
    plot_data <- plot_data[names(selected_genes),]
    if (z_score) {
      plot_data <- t(scale(t(plot_data)))
      col_vec <- circlize::colorRamp2(c(-2, 0, 2),  
                                      c("blue", "white", "red")
                                                                     
                                      )
    
    } else {
      col_vec <- circlize::colorRamp2(c(min(plot_data), median(plot_data), quq(plot_data)),
                                      c("blue", "white", "red"))
    }

    # Create Column Heatmapannotation from coldata(vst_dds
    cluster_data <- vst_data[selected_genes, ]
    #row_clust <- hclust(dist(cluster_data),method = "average")
    #col_clust <- hclust(dist(t(cluster_data)), method = "average")
    
    top_anno <- HeatmapAnnotation(df = coldata[, coldata_to_plot],
                                  col = ann_col)
    
    print(dim(plot_data))
    print(dim(cluster_data))
    if(nrow(plot_data) == 0 ) {
      print("hit")
      return()
    }
    o1 = seriate(dist(plot_data), method = "DendSer")
    o2 = seriate(dist(t(plot_data)), method = "DendSer")
    if (!is.null(col_km)) {
      col_dend <- dendextend::color_branches(as.dendrogram(o2[[1]]), col_km)
    } else {
      col_dend <- as.dendrogram(o2[[1]])
    }
    hmap <- ComplexHeatmap::Heatmap(
      plot_data,
      show_row_names = nrow(plot_data) <= 50,
      top_annotation = top_anno,
      show_column_names = FALSE,
   #   cluster_rows = as.dendrogram(o1[[1]]),
      cluster_columns = col_dend,
      name = "z-scaled expression",
      cluster_rows = dendextend::color_branches(as.dendrogram(o1[[1]]), row_km),
      row_split = row_km,
      column_split = col_km
      # cluster methods for rows and columns
      #     clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
      #      clustering_method_columns = 'ward.D2',
      #      clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
      #     clustering_method_rows = 'ward.D2', col = col_vec,
    )
    # Create row annotation for coefficient values
    anno_df <- data.frame(coef =  coef_values[selected_genes] )
    colnames(anno_df) <- coef[1]
    col_list <- list()
    col_list[[coef[1]]] <-colorRamp2(c(min(coef_values[selected_genes]), 0, 
                                       max(coef_values[selected_genes])), c("blue", "white", "red"))
    coef_values_anno <- rowAnnotation(
      df = anno_df,
      col = col_list,
      annotation_name_rot = 45
    )

    # Create a legend for the coefficient values
    coef_legend <- Legend(
    title = coef[1],
    col_fun = colorRamp2(c(min(coef_values[selected_genes]), 0, max(coef_values[selected_genes])), c("blue", "white", "red")),
     at = c(min(coef_values[selected_genes]), 0, max(coef_values[selected_genes])),
    labels = c(round(min(coef_values[selected_genes]), 2), 0, round(max(coef_values[selected_genes]), 2))
    )

    if (nrow(plot_data) > 1) {
      genes_to_plot <- names(head(sort(abs(coef_values[selected_genes]), decreasing = T),n_genes))
      names_ordered <- rownames(plot_data)[rownames(plot_data) %in% genes_to_plot]
      at_vec <- which(rownames(plot_data) %in% names_ordered)
        genelabels <- rowAnnotation(
        Genes = anno_mark(
        at = at_vec,
        labels = names_ordered,
        labels_gp = gpar(fontsize = 10, fontface = 'bold')),
        width = unit(2.0, 'cm') +
      max_text_width(
        rownames(plot_data)[seq(1, nrow(plot_data), 20)],
        gp = gpar(fontsize = 10,  fontface = 'bold')))

      png(png_file, width = 10, height = 10, units = 'in', res = 300) 
      draw(hmap + genelabels + coef_values_anno,
      heatmap_legend_side = 'left',
      annotation_legend_side = 'right',
      row_sub_title_side = 'left')
      #annotation_legend_list = list(coef_legend))
      dev.off()
    } else {
      png(png_file, width = 10, height = 10, units = 'in', res = 300)
      draw(hmap + coef_values_anno, annotation_legend_list = list(coef_legend))
      dev.off()
    }
}

# Read_glmmseq_obj
glmmseq_obj <- readRDS(glmmseq_obj)
glmmseq_obj$norm_counts <- glmmQvals(glmmseq_obj$norm_counts)
# Read DDs Object 
vst_dds <- readRDS(corrected_counts)
# run_vst
vst_dds <- vst(vst_dds[rownames(glmmseq_obj$norm_counts@countdata),])

plot_res <- plot_heatmap(
  coef = c(var, coef_id),
  vst_obj = vst_dds,
  glmm_obj = glmmseq_obj$norm_counts,
  z_score = TRUE,
  use_vst = TRUE,
  coldata_to_plot = col_data_to_plot,
  n_genes = 30,
  qval_filt = 0.01,
  coef_filt = 0.8,
  meanExp_cutoff = 7,
  row_km = 3,
  col_km = NULL
)
if(is.null(plot_res)) {
  plot_res <- plot_heatmap(
    coef = c(var, coef_id),
    vst_obj = vst_dds,
    glmm_obj = glmmseq_obj$norm_counts,
    z_score = TRUE,
    use_vst = TRUE,
    coldata_to_plot = col_data_to_plot,
    n_genes = 30,
    qval_filt = 0.01,
    coef_filt = 0,
    meanExp_cutoff = 0,
    row_km = 3,
    col_km = NULL
  )
}
if(is.null(plot_res)) {
  write.table("", file = png_file)
}
