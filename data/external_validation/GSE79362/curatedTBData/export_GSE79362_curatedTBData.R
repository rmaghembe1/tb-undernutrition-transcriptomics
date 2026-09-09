options(timeout = 3600)

suppressPackageStartupMessages({
  library(curatedTBData)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

cat("Loading curated GSE79362 expression matrix and metadata...\n")

gse_list <- curatedTBData(
  "GSE79362",
  dry.run = FALSE,
  curated.only = TRUE
)

gse <- gse_list[["GSE79362"]]

if (is.null(gse)) {
  stop("GSE79362 object was not found in the downloaded result list.")
}

cat("\nObject class:\n")
print(class(gse))

exps <- experiments(gse)

cat("\nAvailable experiments:\n")
print(names(exps))

if (!"assay_curated" %in% names(exps)) {
  stop("assay_curated was not present. Available experiments: ",
       paste(names(exps), collapse = ", "))
}

expr <- exps[["assay_curated"]]
meta <- as.data.frame(colData(gse))

if (!is.matrix(expr)) {
  expr <- as.matrix(expr)
}

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nMetadata dimensions before alignment:\n")
print(dim(meta))

cat("\nMetadata columns:\n")
print(colnames(meta))

# Align clinical annotation to expression columns.
if (all(colnames(expr) %in% rownames(meta))) {
  meta <- meta[colnames(expr), , drop = FALSE]
} else {
  cat("\nExpression columns are not direct colData row names; inspecting sampleMap...\n")
  smap <- as.data.frame(sampleMap(gse))
  write.table(
    smap,
    file = "GSE79362_sampleMap.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )

  matched_primary <- smap$primary[match(colnames(expr), smap$colname)]

  if (any(is.na(matched_primary)) || !all(matched_primary %in% rownames(meta))) {
    stop("Could not align expression columns to metadata using sampleMap.")
  }

  meta <- meta[matched_primary, , drop = FALSE]
  rownames(meta) <- colnames(expr)
}

if (ncol(expr) != nrow(meta)) {
  stop("Expression columns and aligned metadata rows do not match.")
}

cat("\nMetadata dimensions after alignment:\n")
print(dim(meta))

cat("\nExpression column examples:\n")
print(head(colnames(expr)))

cat("\nAligned metadata row examples:\n")
print(head(rownames(meta)))

write.table(
  cbind(gene_symbol = rownames(expr), expr),
  file = gzfile("GSE79362_assay_curated.tsv.gz"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  cbind(sample_id = rownames(meta), meta),
  file = "GSE79362_colData.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

capture.output(
  {
    cat("GSE79362 curatedTBData export audit\n\n")

    cat("Object class:\n")
    print(class(gse))

    cat("\nAvailable experiments:\n")
    print(names(exps))

    cat("\nExpression dimensions:\n")
    print(dim(expr))

    cat("\nMetadata dimensions:\n")
    print(dim(meta))

    cat("\nMetadata columns:\n")
    print(colnames(meta))

    cat("\nFirst six aligned metadata rows:\n")
    print(head(meta))

    cat("\nUnique values by metadata variable:\n")
    for (nm in colnames(meta)) {
      cat("\n--- ", nm, " ---\n", sep = "")
      vals <- unique(as.character(meta[[nm]]))
      print(head(vals, 30))
    }

    cat("\nSession information:\n")
    print(sessionInfo())
  },
  file = "GSE79362_export_audit.txt"
)

cat("\nExport complete:\n")
cat("  GSE79362_assay_curated.tsv.gz\n")
cat("  GSE79362_colData.tsv\n")
cat("  GSE79362_export_audit.txt\n")
