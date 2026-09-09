options(timeout = 3600)

suppressPackageStartupMessages({
  library(curatedTBData)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

study <- "GSE94438"
cat("Loading curated data for ", study, "...\n", sep = "")

# Availability record
dry <- tryCatch(
  curatedTBData(study, dry.run = TRUE, curated.only = TRUE),
  error = function(e) e
)
capture.output(dry, file = "GSE94438_dry_run_resources.txt")
if (inherits(dry, "error")) {
  stop("Dry-run failed. Review GSE94438_dry_run_resources.txt.")
}

downloaded <- curatedTBData(study, dry.run = FALSE, curated.only = TRUE)

# Handle either a named list return or a direct MultiAssayExperiment return.
if (methods::is(downloaded, "MultiAssayExperiment")) {
  gse <- downloaded
} else if (!is.null(downloaded[[study]])) {
  gse <- downloaded[[study]]
} else {
  stop("Could not locate a MultiAssayExperiment for GSE94438 in returned object.")
}

exps <- experiments(gse)
exp_names <- names(exps)

cat("Available experiments:\n")
print(exp_names)

if (!"assay_curated" %in% exp_names) {
  capture.output(
    {
      cat("Object class:\n"); print(class(gse))
      cat("\nAvailable experiments:\n"); print(exp_names)
      cat("\nObject summary:\n"); print(gse)
    },
    file = "GSE94438_object_audit.txt"
  )
  stop("assay_curated is unavailable. Review GSE94438_object_audit.txt.")
}

expr_obj <- exps[["assay_curated"]]
if (methods::is(expr_obj, "SummarizedExperiment")) {
  expr <- assay(expr_obj)
} else {
  expr <- as.matrix(expr_obj)
}
meta <- as.data.frame(colData(gse))

# Align expression samples with participant metadata.
if (all(colnames(expr) %in% rownames(meta))) {
  meta_aligned <- meta[colnames(expr), , drop = FALSE]
} else {
  smap <- as.data.frame(sampleMap(gse))
  write.table(smap, "GSE94438_sampleMap.tsv", sep = "\t",
              quote = FALSE, row.names = FALSE, na = "")
  map_rows <- smap[smap$assay == "assay_curated", , drop = FALSE]
  primary <- map_rows$primary[match(colnames(expr), map_rows$colname)]
  if (any(is.na(primary)) || !all(primary %in% rownames(meta))) {
    stop("Could not align assay_curated columns with colData using sampleMap.")
  }
  meta_aligned <- meta[primary, , drop = FALSE]
  rownames(meta_aligned) <- colnames(expr)
}

stopifnot(ncol(expr) == nrow(meta_aligned))

write.table(
  cbind(gene_symbol = rownames(expr), expr),
  gzfile("GSE94438_assay_curated.tsv.gz"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  cbind(sample_id = rownames(meta_aligned), meta_aligned),
  "GSE94438_colData.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

capture.output(
  {
    cat("GSE94438 curatedTBData export audit\n\n")
    cat("Object class:\n"); print(class(gse))
    cat("\nAvailable experiments:\n"); print(exp_names)
    cat("\nExpression dimensions:\n"); print(dim(expr))
    cat("\nMetadata dimensions:\n"); print(dim(meta_aligned))
    cat("\nMetadata columns:\n"); print(colnames(meta_aligned))
    cat("\nFirst six metadata rows:\n"); print(head(meta_aligned))
    cat("\nUnique values by metadata variable:\n")
    for (nm in colnames(meta_aligned)) {
      cat("\n--- ", nm, " ---\n", sep = "")
      print(head(unique(as.character(meta_aligned[[nm]])), 40))
    }
    cat("\nSession information:\n"); print(sessionInfo())
  },
  file = "GSE94438_export_audit.txt"
)

cat("\nExport complete:\n")
cat("  GSE94438_assay_curated.tsv.gz\n")
cat("  GSE94438_colData.tsv\n")
cat("  GSE94438_export_audit.txt\n")
