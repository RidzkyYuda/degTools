#' Run differential expression and TE analysis
#'
#' @param seurat_obj A Seurat object
#' @param cohort_name Character
#' @param cluster_name Character (cell_type to subset)
#' @param metadata_df Data frame with gene/TE annotations (must contain "gene" and "feature_type")
#' @param outdir Output directory (default = current working dir)
#'
#' @return NULL
#' @export

run_DEG <- function(seurat_obj,
                    cohort_name,
                    cluster_name,
                    metadata_df,
                    outdir = ".") {

  ## Memory setting
  options(future.globals.maxSize = 26000 * 1024^2)

  ## -----------------------------
  ## 1. Subset cells safely
  ## -----------------------------
  obj_sub <- subset(seurat_obj, subset = cell_type == cluster_name)

  if (ncol(obj_sub) == 0) {
    stop("No cells found for cluster: ", cluster_name)
  }

  ## Check condition balance
  cond_table <- table(obj_sub$condition)
  message("Condition counts: ", paste(names(cond_table), cond_table, collapse = ", "))

  if (any(cond_table < 10)) {
    warning("Very small group size detected: ",
            paste(names(cond_table), cond_table, collapse = ", "))
  }

  ## -----------------------------
  ## 2. Set assay and identities
  ## -----------------------------
  DefaultAssay(obj_sub) <- "RNA"
  Idents(obj_sub) <- "condition"

  ## -----------------------------
  ## 3. Differential expression
  ## -----------------------------
  DEG <- Seurat::FindMarkers(
    obj_sub,
    ident.1 = "Alcohol",
    ident.2 = "Control",
    test.use = "wilcox",
    assay = "RNA",
    verbose = FALSE
  )

  if (nrow(DEG) == 0) {
    stop("FindMarkers returned no results.")
  }

  ## -----------------------------
  ## 4. Annotation + significance
  ## -----------------------------
  DEG_anno <- DEG %>%
    dplyr::mutate(gene = rownames(DEG)) %>%
    dplyr::left_join(metadata_df, by = "gene") %>%
    dplyr::mutate(
      feature_type = ifelse(is.na(feature_type), "Unknown", feature_type),
      is_signif = p_val_adj < 0.05 & !is.na(p_val_adj)
    )

  ## Sanity check
  message("Feature type distribution:")
  print(table(DEG_anno$feature_type))

  ## -----------------------------
  ## 5. Split outputs
  ## -----------------------------
  DEGonly <- DEG_anno %>%
    dplyr::filter(feature_type == "Gene", is_signif)

  DETE <- DEG_anno %>%
    dplyr::filter(feature_type == "TE", is_signif)

  ## -----------------------------
  ## 6. Export Excel
  ## -----------------------------
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

  ## -----------------------------
  ## 7. Volcano plots
  ## -----------------------------

  ## ---- TE Volcano ----
  DETE_volcano <- DEG_anno %>%
    dplyr::filter(feature_type == "TE") %>%
    dplyr::filter(!is.na(p_val_adj), !is.na(avg_log2FC))

  if (nrow(DETE_volcano) > 0) {

    # optional: label top TE only
    top_TE <- DETE_volcano %>%
      dplyr::arrange(p_val_adj) %>%
      dplyr::slice(1:20) %>%
      dplyr::pull(gene)

    volcano_plot_TE <- EnhancedVolcano::EnhancedVolcano(
      DETE_volcano,
      lab = DETE_volcano$gene,
      selectLab = top_TE,
      x = "avg_log2FC",
      y = "p_val_adj",
      pCutoff = 0.05,
      FCcutoff = 0.1,
      labSize = 3
    )

    ggplot2::ggsave(
      file.path(outdir, paste0("AllCohorts.", prefix, ".TE.volcano.pdf")),
      volcano_plot_TE,
      width = 7,
      height = 6
    )

  } else {
    message("No TE features available for volcano plot.")
  }

  ## ---- Gene Volcano ----
  DEG_volcano <- DEG_anno %>%
    dplyr::filter(feature_type == "Gene") %>%
    dplyr::filter(!is.na(p_val_adj), !is.na(avg_log2FC))

  if (nrow(DEG_volcano) > 0) {

    top_genes <- DEG_volcano %>%
      dplyr::arrange(p_val_adj) %>%
      dplyr::slice(1:20) %>%
      dplyr::pull(gene)

    volcano_plot_genes <- EnhancedVolcano::EnhancedVolcano(
      DEG_volcano,
      lab = DEG_volcano$gene,
      selectLab = top_genes,
      x = "avg_log2FC",
      y = "p_val_adj",
      pCutoff = 0.05,
      FCcutoff = 0.25,
      labSize = 3
    )

    ggplot2::ggsave(
      file.path(outdir, paste0("AllCohorts.", prefix, ".genes.volcano.pdf")),
      volcano_plot_genes,
      width = 7,
      height = 6
    )

  } else {
    message("No gene features available for volcano plot.")
  }

  message("Analysis completed for: ", cluster_name)
}
