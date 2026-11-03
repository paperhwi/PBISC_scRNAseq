#!/usr/bin/env Rscript
# CellDEEP metacell source DEG, with donor-level limma-voom as a robustness check.
# Direction of all logFC values: New LoRA minus Real, within cohort/condition/cell type.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(CellDEEP)
  library(edgeR)
  library(limma)
})

root <- "F:/PBISC/NewLora/analysis_20261102/celldeep_source_deg"
dir.create(file.path(root, "results"), recursive=TRUE, showWarnings=FALSE)

one_stratum <- function(cohort, stratum, index) {
  path <- file.path(root, cohort, stratum)
  meta <- read.csv(file.path(path, "metadata.csv"), stringsAsFactors=FALSE)
  genes <- scan(file.path(root, cohort, "genes.txt"), what="character", quiet=TRUE)
  counts <- as(readMM(file.path(path, "counts.mtx")), "dgCMatrix")
  rownames(counts) <- genes
  colnames(counts) <- meta$cell_id
  rownames(meta) <- meta$cell_id
  stopifnot(nrow(counts)==length(genes), ncol(counts)==nrow(meta))
  stopifnot(all(counts@x==round(counts@x)), all(meta$source %in% c("real", "new")))
  stopifnot(length(unique(meta$sample_id[meta$source=="real"])) ==
            length(unique(meta$sample_id[meta$source=="new"])))

  # CellDEEP pools ten random cells within source × celltype × biological sample.
  # The prepared input has exactly 30 cells per sample, so each sample yields
  # three metacells. CellDEEP/Seurat's Wilcoxon test treats those metacells as
  # independent; the sample-level analysis below checks this limitation.
  set.seed(20261103L + index)
  obj <- CreateSeuratObject(counts=counts, meta.data=meta,
                            min.cells=0, min.features=0)
  cd <- FindMarker.CellDEEP(
    object=obj, ident.1="new", ident.2="real", group.by="group_id",
    sample_id="sample_id", group_id="source", cluster_id="celltype",
    prepare=TRUE, test.use="wilcox", Pool=TRUE, readcounts="sum",
    n_cells=10, assay="RNA", min_cells_per_subgroup=25,
    cell_selection="random", name.only=FALSE, logfc.threshold=0,
    min.pct=0, full_list=TRUE)
  cd <- data.frame(gene=rownames(cd), cd, row.names=NULL, check.names=FALSE)
  cd$q_BH <- p.adjust(cd$p_val, method="BH")
  cd$deg_CellDEEP <- is.finite(cd$q_BH) & cd$q_BH<0.05 &
                     abs(cd$avg_log2FC)>=1

  # Pseudobulk uses one count vector per original biological sample.
  key <- paste(meta$source, meta$sample_id, sep="::")
  unit <- unique(key)
  incidence <- sparseMatrix(i=seq_along(key), j=match(key, unit), x=1,
                            dims=c(length(key), length(unit)))
  pb <- counts %*% incidence
  source <- factor(sub("::.*$", "", unit), levels=c("real", "new"))
  design <- model.matrix(~source)
  dge <- DGEList(counts=pb, group=source)
  keep <- filterByExpr(dge, design=design)
  pb_full <- data.frame(gene=genes, pb_logFC=NA_real_, pb_P=NA_real_,
                        pb_q_BH=NA_real_, pb_AveExpr=NA_real_,
                        pb_tested=keep, stringsAsFactors=FALSE)
  if (sum(keep)>=100) {
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=FALSE)
    fit <- eBayes(lmFit(v, design))
    tb <- topTable(fit, coef="sourcenew", number=Inf, sort.by="none")
    match_idx <- match(rownames(tb), pb_full$gene)
    pb_full$pb_logFC[match_idx] <- tb$logFC
    pb_full$pb_P[match_idx] <- tb$P.Value
    pb_full$pb_q_BH[match_idx] <- tb$adj.P.Val
    pb_full$pb_AveExpr[match_idx] <- tb$AveExpr
  }
  pb_full$deg_pseudobulk <- !is.na(pb_full$pb_q_BH) &
                            pb_full$pb_q_BH<0.05 & abs(pb_full$pb_logFC)>=1
  result <- merge(pb_full, cd, by="gene", all.x=TRUE, sort=FALSE)
  result$cohort <- cohort
  result$stratum <- stratum
  result$condition <- meta$group[1]
  result$celltype <- meta$celltype[1]
  result$source_effect <- "new_minus_real"
  out_file <- file.path(root, "results", paste0(stratum, "_deg.csv.gz"))
  con <- gzfile(out_file, open="wt")
  write.csv(result, con, row.names=FALSE, na="")
  close(con)

  valid <- is.finite(result$avg_log2FC) & is.finite(result$pb_logFC)
  overlap <- result$deg_CellDEEP & result$deg_pseudobulk &
             sign(result$avg_log2FC)==sign(result$pb_logFC)
  data.frame(cohort=cohort, stratum=stratum, condition=meta$group[1],
    celltype=meta$celltype[1], n_real_samples=sum(source=="real"),
    n_new_samples=sum(source=="new"),
    n_cells_real=sum(meta$source=="real"), n_cells_new=sum(meta$source=="new"),
    n_metacells_each=3*sum(source=="real"),
    n_genes_CellDEEP=sum(is.finite(result$p_val)),
    n_genes_pseudobulk=sum(result$pb_tested),
    n_DEG_CellDEEP=sum(result$deg_CellDEEP, na.rm=TRUE),
    n_DEG_pseudobulk=sum(result$deg_pseudobulk, na.rm=TRUE),
    n_overlap_same_direction=sum(overlap, na.rm=TRUE),
    logFC_pearson=if (sum(valid)>10) cor(result$avg_log2FC[valid],result$pb_logFC[valid]) else NA_real_)
}

summary_path <- file.path(root, "results", "summary.csv")
status <- list()
for (cohort in c("sepsis", "covid")) {
  plan <- read.csv(file.path(root, cohort, "strata.csv"), stringsAsFactors=FALSE)
  plan <- plan[plan$status=="ready", ]
  for (i in seq_len(nrow(plan))) {
    key <- plan$stratum[i]
    message("RUN ", key, " (", i, "/", nrow(plan), ")")
    started <- Sys.time()
    result <- tryCatch(one_stratum(cohort, key, i),
      error=function(e) {
        message("ERROR ", key, ": ", conditionMessage(e))
        data.frame(cohort=cohort, stratum=key, error=conditionMessage(e))
      })
    result$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units="secs"))
    status[[length(status)+1]] <- result
    write.csv(data.table::rbindlist(status, fill=TRUE), summary_path, row.names=FALSE)
    gc(verbose=FALSE)
  }
}
message("Saved ", summary_path)
