#' @importFrom rlang inform
#' @importFrom grDevices rainbow
# Choose appropriate palette based on n
.updatePalette <- function(n, quiet = FALSE) {
  if (n > 40) {
    if (!quiet) inform(sprintf("Use `rainbow(%d)`.", n))
    color <- rainbow(n) # use rainbow() when more than 40 colurs
  } else if (n > 30) {
    if (!quiet) inform(sprintf("Use `c40(%d)`.", n))
    color <- c40()[1:n] # 40 colours
  } else {
    if (!quiet) inform(sprintf("Use `c30(%d)`.", n))
    color <- c30()[1:n]
  }
  color
}

#' @importFrom rlang abort
# Build export path
.build.path <- function(dir_path = NULL) {
  if (is.null(dir_path)) {
    dir_path <- getwd()
  } else {
    if (!dir.exists(dir_path)) abort("The provided path does not exist.")
  }
  dir_path
}

#' @importFrom rlang abort
# Calculate row means
.clac_rowMeans <- function(x) {
  if (requireNamespace("Matrix")) Matrix::rowMeans(x) else abort("Requires CRAN package 'Matrix'.")
}

#' @importFrom rlang abort
# Calculate row variance
.clac_rowVars <- function(x) {
  if (requireNamespace("DelayedMatrixStats")) DelayedMatrixStats::rowVars(x) else abort("Requires Bioconductor package 'DelayedMatrixStats'.")
}

#' @importFrom utils tail
# Retrieve the value from the p.value column of the last row of a data.frame
.get_lastp <- function(x) {
  ifelse(nrow(x) > 0, tail(x, 1)$p.value, NA)
}

# Return a vector containing unique labels and frequencies
.count_label <- function(x) {
  freq <- table(as.character(x))
  sprintf("%s (%d)", names(freq), as.numeric(freq))
}

#' @importFrom rlang abort
.check_is.null <- function(x) {
  var <- deparse(substitute(x))
  if (is.null(x)) abort(sprintf("Argument `%s` is required.", var))
}

#' @importFrom rlang abort
.is.sce <- function(x) {
  .check_is.null(x)
  var <- deparse(substitute(x))
  if (class(x) != "SingleCellExperiment") abort(sprintf("`%s` must be a SingleCellExperiment object.", var))
}

#' @importFrom rlang abort
.check_phase <- function(x) {
  var <- deparse(substitute(x))
  if (!is.list(x)) abort(sprintf("`%s` is not a list.", var))

  list_names <- c("phases", "scores", "normalized.scores")
  if (!identical(names(x), list_names)) abort("Please provide output from `scran::cyclone`.")
}

#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom rlang abort
.check_dimname <- function(sce, dimname) {
  .check_is.null(dimname)
  var <- deparse(substitute(sce))
  if (!dimname %in% reducedDimNames(sce)) abort(sprintf("Cannot find reducedDimNames `%s` in `%s`.", dimname, var))
}

#' @importFrom rlang abort
.check_wholenum <- function(x, abort = TRUE) {
  .check_is.null(x)
  var <- deparse(substitute(x))
  msg <- if (length(x) > 1L) "vector" else "scalar"
  if (!is.numeric(x) || !sum(x - floor(x)) == 0) {
    if (abort) abort(sprintf("Argument `%s` is invalid, integer %s required.", var, msg)) else FALSE
  } else {
    TRUE
  }
}

#' @importFrom rlang abort
.check_column_by <- function(df, x) {
  if (!is.null(x)) {
    if (length(x) > 1) {
      abort("`column_by` should contain name of a single column .")
    } else if (!x %in% colnames(df)) {
      abort(sprintf("Cannot find column named '%s' in result.", x))
    }
  }
}

#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @importFrom rlang abort
.check_assayname <- function(sce, exprs_by, return_value = FALSE) {
  var <- deparse(substitute(sce))
  if (!exprs_by %in% assayNames(sce)) abort(sprintf("Cannot find assay named '%s' in `%s`.", exprs_by, sce))
  if (return_value) {
    return(assay(sce, exprs_by))
  }
}

#' @importFrom SingleCellExperiment colData
#' @importFrom rlang abort
.check_feature <- function(sce, x, exprs_by) {
  .check_is.null(x)
  var <- deparse(substitute(sce))

  if (length(x) == 1) {
    if (x %in% colnames(colData(sce))) { # cell feature
      val <- colData(sce)[, x]
    } else if (x %in% rownames(sce)) { # gene expression
      val <- .check_assayname(sce, exprs_by, return_value = TRUE)
      val <- val[x, ]
    } else {
      abort(sprintf("Cannot find the specified `feature` in `%s`.", sce))
    }
  } else {
    if (length(x) != ncol(sce)) abort(sprintf("`feature` and `ncol(%s)` have different lengths.", sce))
    val <- x
    x <- "feature"
  }
  list(name = x, values = val)
}

#' @importFrom rlang abort
.check_res.type <- function(type, switch = NULL) {
  if (is.null(type)) abort("Unknown `object` type, supported source: findMarkers, edgeR and DESeq2.")

  if (switch == "is.list") { # check if a list
    if (!attr(type, "is.list")) abort("`object` is not an object of class `list` or `SimpleList`.")
  } else if (switch == "marker") { # check if from findMarkers
    if (type != "findMarkers") abort("`object` is not a findMarkers output.")
  } else if (switch == "deg") { # check if from edgeR/DESeq2
    if (!type %in% c("edgeR", "DESeq2")) abort("`object` is not a edgeR/DESeq2 output.")
  } else if (switch == "enrichr") {
    if (type != "EnrichR") abort("`object` is not a Enrichr output.")
  }
}

#' @importFrom rlang abort
.check_2nums <- function(x) {
  var <- deparse(substitute(x))
  if (length(x) != 2L || !is.numeric(x)) abort(sprintf("Wrong `%s` format.", var))
}

#' @importFrom rlang abort
.check_direction <- function(direction) {
  if (!direction %in% c("both", "up", "down")) abort("Wrong `direction` input.")
}

#' @importFrom rlang abort
# Determine the type of the input
.get_res.type <- function(x) {
  if (class(x) %in% c("list", "SimpleList")) {
    n <- length(Reduce(union, lapply(x, .get_type, is.list = TRUE)))
    if (n > 1) abort("The list contains results from more than one source.")
    .get_type(x[[1]], is.list = TRUE)
  } else {
    .get_type(x, is.list = FALSE)
  }
}

# Determine the type of the input, actual
.get_type <- function(x, is.list) {
  obj <- as.character(class(x))

  type <- NULL
  if (class(x) %in% c("TopTags", "DGEExact", "DGELRT")) {
    type <- structure("edgeR", is.list = is.list, obj = obj)
  } else if (class(x) == "DESeqResults") {
    type <- structure("DESeq2", is.list = is.list, obj = obj)
  } else if ("summary.logFC" %in% colnames(x) || "summary.AUC" %in% colnames(x)) {
    type <- structure("findMarkers",
      obj = obj, is.list = is.list,
      test.type = ifelse("summary.AUC" %in% colnames(x), "wilcox", "t or binom"),
      pval.type = ifelse("Top" %in% colnames(x), "any", "all or some")
    )
  } else if (length(intersect(c(
    "Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value",
    "Old.Adjusted.P.value", "Combined.Score", "Genes"
  ), colnames(x))) >= 8L) { # a list
    type <- structure("EnrichR", obj = obj, is.list = is.list)
  } else if (length(intersect(c(
    "Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value",
    "Old.Adjusted.P.value", "Combined.Score", "Genes"
  ), colnames(x[[1]]))) >= 8L) { # a list of lists
    type <- structure("EnrichR", obj = obj, is.list = is.list)
  }
  type
}

# Get test columns from findMarkers/DESeq2/edgeR results
.get_test.cols <- function(x, type) {
  if (type == "findMarkers") {
    summary <- if (attr(type, "test.type") == "wilcox") "summary.AUC" else "summary.logFC"
    cols <- if (attr(type, "pval.type") == "any") c("Top", summary, "p.value", "FDR") else c(summary, "p.value", "FDR")
  } else if (type == "edgeR") {
    # exactTest:	logFC	logCPM	PValue	FDR
    # glmTreat (prior.count==0):	logFC	logCPM	PValue	FDR
    # glmTreat (prior.count!=0):	logFC	unshrunk.logFC	logCPM	PValue	FDR
    # glmQLFTest:	logFC	logCPM	F	PValue	FDR
    # glmLRT:	logFC	logCPM	LR	PValue	FDR
    cols <- c("logFC", "logCPM", "PValue", "FDR")
  } else if (type == "DESeq2") {
    cols <- c("baseMean", "log2FoldChange", "pvalue", "padj")
  }
  .res2df(x, type)[, cols]
}

#' @importFrom rlang abort
# Convert findMarkers/DESeq2/edgeR results to DataFrame
.res2df <- function(x, type) {
  if (type == "edgeR") {
    if (requireNamespace("edgeR")) {
      if (attr(type, "obj") == "TopTags") {
        as.data.frame(x$table)
      } else {
        as.data.frame(edgeR::topTags(x, n = nrow(x))$table)
      }
    } else {
      abort("Requires Bioconductor package 'edgeR'.")
    }
  } else if (type == "DESeq2") {
    if (requireNamespace("DESeq2")) {
      df <- as.data.frame(x)
    } else {
      abort("Requires Bioconductor package 'DESeq2'.")
    }
  } else {
    as.data.frame(x)
  }
}

# Retrieve the index in findMarkers output DataFrame
.get_fmidx <- function(x, direction, top, fdr, test.type, pval.type, logfc, auc) {
  # Order smaller value first
  logfc <- sort(logfc)
  auc <- sort(auc)

  if (test.type == "wilcox") {
    if (pval.type == "any") {
      up <- which(x$Top <= top & x$summary.AUC > auc[2])
      dn <- which(x$Top <= top & x$summary.AUC < auc[1])
    } else {
      up <- which(x$FDR <= fdr & x$summary.AUC > auc[2])
      dn <- which(x$FDR <= fdr & x$summary.AUC < auc[1])
    }
  } else {
    if (pval.type == "any") {
      up <- which(x$Top <= top & x$summary.logFC > logfc[2])
      dn <- which(x$Top <= top & x$summary.logFC < logfc[1])
    } else {
      up <- which(x$FDR <= fdr & x$summary.logFC > logfc[2])
      dn <- which(x$FDR <= fdr & x$summary.logFC < logfc[1])
    }
  }

  idx <- sort(c(up, dn))
  if (direction == "up") {
    idx <- up
  } else if (direction == "down") {
    idx <- dn
  }
  idx
}

#' @importFrom enrichR setEnrichrSite
.resetEnrichrSite <- function(site) {
  if (is.null(getOption("enrichR.base.address"))) options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
  if (is.null(getOption("enrichR.live"))) options(enrichR.live = TRUE)
  if (is.null(getOption("modEnrichR.use"))) options(modEnrichR.use = TRUE)
  if (is.null(getOption("enrichR.sites.base.address"))) options(enrichR.sites.base.address = "https://maayanlab.cloud/")
  if (is.null(getOption("enrichR.sites"))) options(enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr", "OxEnrichr"))
  setEnrichrSite(site)
}

#' @importFrom rlang abort
# Same function in enrichR 3.0 but with bug fixed
.enrichment_prep_df <- function(df, showTerms, orderBy) {
  if (is.null(showTerms)) {
    showTerms <- nrow(df)
  } else {
    .check_wholenum(showTerms)
  }

  Annotated <- as.numeric(sub("^\\d+/", "", as.character(df$Overlap)))
  Significant <- as.numeric(sub("/\\d+$", "", as.character(df$Overlap)))

  # Build data frame
  df <- cbind(df, data.frame(
    Annotated = Annotated, Significant = Significant,
    stringsAsFactors = FALSE
  ))

  # Order data frame (P.value or Combined.Score)
  if (orderBy == "Combined.Score") {
    idx <- order(df$Combined.Score, decreasing = TRUE)
  } else if (orderBy == "Adjusted.P.value") {
    idx <- order(df$Adjusted.P.value, decreasing = FALSE)
  } else {
    idx <- order(df$P.value, decreasing = FALSE)
  }
  df <- df[idx, ]

  # Subset to selected number of terms
  if (showTerms <= nrow(df)) df <- df[1:showTerms, ]
  df
}
