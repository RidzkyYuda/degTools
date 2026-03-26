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
                    cohort_name,
                    cluster_name,
                    metadata_df,
                    outdir = ".") {

  options(future.globals.maxSize = 26000 * 1024^2)

  DefaultAssay(seurat_obj) <- "RNA"
  Idents(seurat_obj) <- "cell_type"

  DEG <- Seurat::FindMarkers(
    seurat_obj,
    ident.1 = "Alcohol",
    ident.2 = "Control",
    test.use = "wilcox",
    subset.ident = cluster_name,
    group.by = "condition",
    assay = "RNA",
    verbose = FALSE
  )

  DEG$gene <- rownames(DEG)

  DEG_anno <- dplyr::left_join(DEG, metadata_df, by = "gene")

  DETE <- dplyr::filter(DEG_anno,
                        feature_type == "TE",
                        p_val_adj < 0.05)

  DEGonly <- dplyr::filter(DEG_anno,
                        feature_type == "Gene",
                        p_val_adj < 0.05)

  prefix <- paste0(cohort_name, ".", cluster_name, ".clust")

  WriteXLS::WriteXLS(
    DEG_anno,
    file.path(outdir, paste0("AllCohorts.", prefix, ".DEG.genes.and.TE.xlsx"))
  )

  WriteXLS::WriteXLS(
    DEGonly,
    file.path(outdir, paste0("AllCohorts.", prefix, ".DEG.genes.only.xlsx"))
  )

  WriteXLS::WriteXLS(
    DETE,
    file.path(outdir, paste0("AllCohorts.", prefix, ".DETE.xlsx"))
  )
  
  ################## create volcano plot for TE ########
  DETE_volcano <- dplyr::filter(DEG_anno, feature_type == "TE")

  volcano_plot_TE <- EnhancedVolcano::EnhancedVolcano(
    DETE_volcano,
    lab = DETE_volcano$gene,
    x = "avg_log2FC",
    y = "p_val_adj"
  )

  ggplot2::ggsave(
    file.path(outdir, paste0("AllCohorts.", prefix, ".TE.volcano.pdf")),
    volcano_plot_TE
  )

  ############### create volcano plot for genes #########
  DEG_volcano <- dplyr::filter(DEG_anno, feature_type == "Gene")
  volcano_plot_genes <- EnhancedVolcano::EnhancedVolcano(
    DEG_volcano,
    lab = DEG_volcano$gene,
    x = "avg_log2FC",
    y = "p_val_adj"
  )

  ggplot2::ggsave(
    file.path(outdir, paste0("AllCohorts.", prefix, ".genes.volcano.pdf")),
    volcano_plot_genes
  )
}
