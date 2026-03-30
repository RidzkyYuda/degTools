#' Run differential expression and TE analysis
#'
#' @param seurat_obj A Seurat object
#' @param cohort_name Character
#' @param cluster_name Character
#' @param metadata_df Data frame with gene annotations
#' @param outdir Output directory (default = current working dir)
#'
#' @return NULL
#' @export
run_DEG <- function(seurat_obj,
                    treatment,
                    control,
                    stat.test,
                    cohort_name,
                    cluster_name,
                    metadata_df,
                    outdir = ".") {

  options(future.globals.maxSize = 26000 * 1024^2)

  DefaultAssay(seurat_obj) <- "RNA"
  Idents(seurat_obj) <- "cell_type"

  # Perform DEG analysis using FindMarkers wilcoxon
  DEG <- Seurat::FindMarkers(
    seurat_obj,
    ident.1 = treatment, # define the treatment group
    ident.2 = control, # define the control group
    test.use = stat.test, # define the statistical test
    subset.ident = cluster_name, # set the cluster to be tested
    group.by = "condition", 
    assay = "RNA",
    verbose = FALSE
  )

  # Annotation added to the DEG results
  DEG_anno <- DEG %>%
  dplyr::mutate(gene = rownames(DEG)) %>%
  dplyr::left_join(metadata_df, by = "gene") %>%
  dplyr::mutate(
    feature_type = ifelse(is.na(feature_type), "Unknown", feature_type),
    is_signif = p_val_adj < 0.05 & !is.na(p_val_adj)
  )

  # Filter differentially expressed TE only
  DETE <- dplyr::filter(DEG_anno,
                        feature_type == "TE",
                        p_val_adj < 0.05)

  # Filter differentially expressed genes only
  DEGonly <- dplyr::filter(DEG_anno,
                        feature_type == "Gene",
                        p_val_adj < 0.05)

  prefix <- paste0(cohort_name, ".", cluster_name, ".clust")

  # save the results in excel table
  
  WriteXLS::WriteXLS(
    "DEG_anno",
    file.path(outdir, paste0(prefix, "_DEG.genes.and.TE.xlsx"))
  )

  WriteXLS::WriteXLS(
    "DEGonly",
    file.path(outdir, paste0(prefix, "_DEG.genes.only.xlsx"))
  )

  WriteXLS::WriteXLS(
    "DETE",
    file.path(outdir, paste0(prefix, "_DETE.xlsx"))
  )
  
  ################## create volcano plot for TE ########
  DETE_volcano <- dplyr::filter(DEG_anno, feature_type == "TE")

  volcano_plot_TE <- EnhancedVolcano::EnhancedVolcano(
    DETE_volcano,
    lab = DETE_volcano$gene,
    x = "avg_log2FC",
    y = "p_val_adj",
    FCcutoff = 0.1,
    pCutoff = 0.05,
    title = paste0(prefix, " TE volcano")
  )

  ggplot2::ggsave(
    file.path(outdir, paste0(prefix, "_TE.pdf")),
    volcano_plot_TE
  )

  ############### create volcano plot for genes #########
  DEG_volcano <- dplyr::filter(DEG_anno, feature_type == "Gene")
  volcano_plot_genes <- EnhancedVolcano::EnhancedVolcano(
    DEG_volcano,
    lab = DEG_volcano$gene,
    x = "avg_log2FC",
    y = "p_val_adj",
    FCcutoff = 0.58,
    pCutoff = 0.05,
    title = paste0(prefix, " Gene volcano")
  )

  ggplot2::ggsave(
    file.path(outdir, paste0(prefix, "_genes.pdf")),
    volcano_plot_genes
  )
}
