#' Print marker gene stats from `findMarkers` output
#'
#' This function parses the `findMarkers()` output (from the `pkg{scran}
#' package) and print the number of marker genes that "passed" the `fdr`
#' or `top` threshold on to standard output in interactive sessions. It
#' provides a quick way to find out the number of candidate marker genes.
#'
#' @param object A named list of DataFrames returned by `findMarkers()`, each
#' of which contains a sorted marker gene list for the corresponding group of
#' comparison.
#' @param fdr A numeric scalar indicating the FDR threashold used to select
#' markers. This threshold has no effect when the input source is from
#' `findMarkers()` ran with `pval.type = "any"`. Default is 0.05.
#' @param top An integer scalar indicating the 'Top' threashold used to select
#' markers. This threshold is only applicable when the input source is from
#' `findMarkers()` ran with `pval.type = "any"`. Default is 200.
#' @param pval.type A string specifying how p-values are to be combined across
#' pairwise comparisons for a given group/cluster when using `findMarkers()`.
#' Default is `NULL`.
#' @param min.prop A numeric scalar specifying the minimum proportion of
#' significant comparisons per gene when using `findMarkers()`. Default is
#' `NULL`.
#' @param logfc A numeric vector of length 2 specifying the logFC threasholds
#' used to select up-/down-regulated markers/DEGs. This threshold has
#' no effect when the input source is from `findMarkers()` ran with
#' `test.type = "wilcox"`. Default is `c(0, 0)`.
#' @param auc A numeric vector of length 2 specifying the AUC threasholds
#' used to select up-/down-regulated markers. This threshold is only
#' applicable when the input source is from `findMarkers()` ran with
#' `test.type = "wilcox"`. Default is `c(0.3, 0.7)`.
#'
#' @return ‘simple’ diagnostic messages
#'
#' @details
#' The function can detect the `pval.type` and `test.type`, to a certain
#' degree, based on the resulting DataFrame structure. One can optionally
#' provide the `pval.type` and `min.prop` settings used to run `findMarkers()`
#' to allow `printMarkerStats()` to use them in the printed messages.
#'
#' @author I-Hsuan Lin
#'
#' @name printMarkerStats
#'
#' @seealso [scran::findMarkers()]
#'
#' @export
#' @importFrom rlang inform
#' @examples
#' library(SingleCellExperiment)
#'
#' # Load demo dataset
#' data(sce)
#'
#' # Use result saved after running `findMarkers()` with `pval.type = "any"`
#' printMarkerStats(metadata(sce)[["findMarkers1"]])
#'
#' # Use result saved after running `findMarkers()` with `pval.type = "all"`
#' printMarkerStats(metadata(sce)[["findMarkers2"]])
#'
#' # Use result saved after running `findMarkers()` with `pval.type = "any"`
#' # and `min.prop = 0.5`
#' printMarkerStats(metadata(sce)[["findMarkers3"]], min.prop = 0.5)
#'
#' # Use result saved after running `findMarkers()` with `test.type = "wilcox"`
#' printMarkerStats(metadata(sce)[["findMarkers4"]])
printMarkerStats <- function(object, fdr = 0.05, top = 200, pval.type = NULL,
                             min.prop = NULL, logfc = c(0, 0), auc = c(0.3, 0.7)) {
  type <- .get_res.type(object)
  .check_res.type(type, "marker")
  .check_2nums(logfc)
  .check_2nums(auc)

  test.type <- attr(type, "test.type")
  pval.type <- attr(type, "pval.type")

  if (pval.type == "any") {
    if (!is.null(min.prop)) {
      inform(sprintf(
        "Number of selected markers (Top %d genes of at least %.1f%% comparisons):",
        top, min.prop * 100
      ))
    } else {
      inform(sprintf("Number of selected markers (Top %d genes of any 1 comparison):", top))
    }
  } else {
    # min.prop is not applicable when pval.type = "all"
    if (!is.null(min.prop) & pval.type != "all") {
      inform(sprintf(
        "Number of selected markers (FDR at %.2f with 'min.prop' at %.2f):",
        fdr, min.prop
      ))
    } else {
      inform(sprintf("Number of selected markers (FDR at %.2f):", fdr))
    }
  }

  for (group in names(object)) {
    df <- .res2df(object[[group]], type)

    up <- .get_fmidx(df, direction = "up", top, fdr, test.type, pval.type, logfc, auc)
    up.n <- length(up)
    up.p <- .get_lastp(df[up, ])

    dn <- .get_fmidx(df, direction = "down", top, fdr, test.type, pval.type, logfc, auc)
    dn.n <- length(dn)
    dn.p <- .get_lastp(df[dn, ])

    n <- up.n + dn.n
    p <- if (is.na(up.p) & is.na(dn.p)) NA else max(c(up.p, dn.p), na.rm = TRUE)

    # If group is a number, add 'Cluster' before the number
    newg <- ifelse(grepl("^[0-9]+$", group), paste0("Cluster", group), group)

    inform(sprintf(
      "- %s: %d; Up = %d; Down = %d; Max. P-value = %.2g.", newg, n, up.n, dn.n, p
    ))
  }

  if (test.type != "wilcox") {
    inform(sprintf(
      "* Upregulated when logFC > %.1f and downregulated when logFC < %.1f.",
      logfc[2], logfc[1]
    ))
  } else {
    inform(sprintf(
      "* Upregulated when AUC > %.1f and downregulated when AUC < %.1f.",
      auc[2], auc[1]
    ))
  }
}

#' Retrieve marker genes from a `findMarkers` output DataFrame
#'
#' This function parses a `findMarkers()` output DataFrame and returns the
#' up- and down-regulated marker genes that passed the specified threshold.
#'
#' @param object A DataFrame that contain a sorted marker gene list for a
#' single group. It should be taken from a named list of DataFrames returned
#' by `findMarkers()`.
#' @param direction A string indicating the directionality of expression
#' changes of markers to retreive. Allowable character values are `"both"`,
#' `"up"` and `"down"`. Default is "both".
#' @param fdr A numeric scalar indicating the FDR threashold used to select
#' markers. This threshold has no effect when the input source is from
#' `findMarkers()` ran with `pval.type = "any"`. Default is 0.05.
#' @param top An integer scalar indicating the 'Top' threashold used to select
#' markers. This threshold is only applicable when the input source is from
#' `findMarkers()` ran with `pval.type = "any"`. Default is 200.
#' @param logfc A numeric vector of length 2 specifying the logFC threasholds
#' used to select up-/down-regulated markers/DEGs. This threshold has
#' no effect when the input source is from `findMarkers()` ran with
#' `test.type = "wilcox"`. Default is `c(0, 0)`.
#' @param auc A numeric vector of length 2 specifying the AUC threasholds
#' used to select up-/down-regulated markers. This threshold is only
#' applicable when the input source is from `findMarkers()` ran with
#' `test.type = "wilcox"`. Default is `c(0.3, 0.7)`.
#' @param max An integer scalar indicating the maximum number of markers to
#' return. Default is 1000.
#' @param column_by A string indicating the column name containing the gene
#' symbols or IDs in the input `object` and the selected information is
#' returned. When `column_by = NULL`, it retreives the row names from the
#' DataFrame. Default is `NULL`.
#'
#' @return A character vector
#'
#' @details
#' The function takes a DataFrame from the `findMarkers()` output (which is
#' a named list of DataFrames) and returns the marker genes that "passed" the
#' `fdr` or `top` threshold.
#'
#' By specifying `direction`, the function selects up-regulated markers when
#' `direction = "up"`, down-regulated markers when `direction = "down"` and
#' both up- and down-regulated markers when `direction = "both"` (default).
#' The definition of up/down-regulation depends on the `test.type` used when
#' running `findMarkers()`.
#'
#' @author I-Hsuan Lin
#'
#' @name getMarkers
#'
#' @seealso [scran::findMarkers()]
#'
#' @export
#' @importFrom utils head
#' @examples
#' library(SingleCellExperiment)
#'
#' # Load demo dataset
#' data(sce)
#'
#' # Retrieve first 1000 markers (default `max`) from the 1st DataFrame
#' # in 'findMarkers1', return rownames
#' getMarkers(metadata(sce)[["findMarkers1"]][[1]])
#'
#' # Retrieve first 1000 up-regulated markers from the 2nd DataFrame
#' # in 'findMarkers4', return rownames
#' getMarkers(metadata(sce)[["findMarkers4"]][[2]], direction = "up")
getMarkers <- function(object, direction = "both", fdr = 0.05, top = 200, logfc = c(0, 0),
                       auc = c(0.3, 0.7), max = 1000, column_by = NULL) {
  type <- .get_res.type(object)
  .check_res.type(type, "marker")
  .check_2nums(logfc)
  .check_2nums(auc)
  .check_direction(direction)

  test.type <- attr(type, "test.type")
  pval.type <- attr(type, "pval.type")

  df <- .res2df(object, type)
  .check_column_by(df, column_by)
  idx <- .get_fmidx(df, direction, top, fdr, test.type, pval.type, logfc, auc)

  # Return genes
  if (is.null(column_by)) head(rownames(df)[idx], max) else head(df[idx, column_by], max)
}

#' Retrieve differentially expressed genes from edgeR or DESeq2 analysis
#'
#' This function parses an object of class `TopTags`, `DGEExact` or `DGELRT`
#' from \pkg{edgeR} or `DESeqResults` class from \pkg{DESeq2}, and returns
#' the differentially expressed genes (DEGs) that passed the specified
#' thresholds.
#'
#' @param object An object of class `TopTags`, `DGEExact`, `DGELRT` or
#' `DESeqResults`.
#' @param direction A string indicating the directionality of expression
#' changes of DEGs to retreive. Allowable character values are `"both"`,
#' `"up"` and `"down"`. Default is "both".
#' @param fdr A numeric scalar indicating the FDR threashold used to select
#' DEGs. Default is 0.05.
#' @param logfc A numeric scalar indicating the logFC cutoffs used to select
#' DEGs. Default is `c(0, 0)`.
#' @param max An integer scalar indicating the maximum number of DEGs to
#' return. Default is 1000.
#' @param column_by A string indicating the column name containing the gene
#' symbols or IDs in the input `object` and the selected information is
#' returned. When `column_by = NULL`, it retreives the row names from the
#' DataFrame. Default is `NULL`.
#'
#' @return A character vector
#'
#' @details
#' The function parses the standard output objects from \pkg{edgeR} and
#' \pkg{DESeq2} to selects DEGs that "passed" the `fdr`, `logfc` and `max`
#' threshold.
#'
#' The accapted class types includes:
#'   - `DGEExact`, returned by `exactTest()` from \pkg{edgeR}.
#'   - `DGELRT`, returned by `glmLRT()`, `glmTreat()` or `glmQLFTest()`
#'     from \pkg{edgeR}.
#'   - `TopTags`, returned by `topTags()` from \pkg{edgeR}.
#'   - `DESeqResults` returned by `results()` from \pkg{DESeq2}.
#'
#' By specifying `direction`, the function selects DEGs with a positive and
#' negative logFC when `direction = "both"` (default), only positive logFC when
#' `direction = "up"`, and only negative logFC when `direction = "down"`.
#'
#' @author I-Hsuan Lin
#'
#' @name getDEGs
#'
#' @seealso [edgeR::exactTest()], [edgeR::glmLRT()], [edgeR::glmTreat()],
#' [edgeR::glmQLFTest()], [edgeR::topTags()], [DESeq2::results()]
#'
#' @export
#' @importFrom utils head
#' @examples
#' # Load demo dataset
#' data(res_deseq2)
#' data(res_edger)
#'
#' # Retrieve first 250 up-regulated DEGs from DESeq2's output,
#' # return rownames (Ensembl ID)
#' getDEGs(res_deseq2, direction = "up", max = 250)
#'
#' # Retrieve first 1000 DEGs at FDR = 10% from edgeR's output,
#' # return values from the `external_gene_name` column (Gene Symbol)
#' getDEGs(res_edger, fdr = 0.1, column_by = "external_gene_name")
getDEGs <- function(object, direction = "both", fdr = 0.05, logfc = c(0, 0), max = 1000,
                    column_by = NULL) {
  type <- .get_res.type(object)
  .check_res.type(type, "deg")
  .check_2nums(logfc)
  .check_direction(direction)

  # Order smaller value first
  logfc <- sort(logfc)

  df <- .res2df(object, type)
  .check_column_by(df, column_by)

  if (type == "edgeR") {
    up <- which(df$FDR < fdr & df$logFC > logfc[2])
    dn <- which(df$FDR < fdr & df$logFC < logfc[1])
  } else {
    up <- which(df$padj < fdr & df$log2FoldChange > logfc[2])
    dn <- which(df$padj < fdr & df$log2FoldChange < logfc[1])
  }

  idx <- sort(c(up, dn))
  if (direction == "up") {
    idx <- up
  } else if (direction == "down") {
    idx <- dn
  }

  # Return genes
  if (is.null(column_by)) head(rownames(df)[idx], max) else head(df[, column_by][idx], max)
}

#' Export results from findMarkers, edgeR or DESeq2 to tsv
#'
#' This function takes a list of objects containing results from
#' `findMarkers()`, \pkg{edgeR} or \pkg{DESeq2}, and save the content
#' to text file(s).
#'
#' @param object A named list of objects storing outputs returned by
#' `findMarkers()`, \pkg{edgeR} or \pkg{DESeq2}, each of which contains
#' results for the corresponding group of comparison.
#' @param concatenate Logical scalar indicating whether to concatenate
#' results stored in a list to a single file. Default is `FALSE`.
#' @param col_anno A character vector containing the names of the annotation
#' columns in the results stored in `object`, for example
#' `c("GeneID","Symbol")`. This is used to construct common annotation
#' column(s) when `concatenate = TRUE`. The row names of the result is used
#' to create a "Symbol" column in the concatenated output when
#' `col_anno = NULL`. Default is `NULL`.
#' @param direction A string indicating the directionality of expression
#' changes of markers/DEGs used to perform marker gene detection or
#' differential expression analysis. The string is added to the output file
#' name. Allowable character values are `"both"`, `"up"` and `"down"`.
#' Default  is `"both"`.
#' @param prefix A string indicating the prefix of output file. Default is
#' `"output"`.
#' @param dir_path The directory path for exported TSV. The file is exported
#' to the current working directory when `file_path = NULL`. Default is `NULL`.
#'
#' @return NULL
#'
#' @details When `concatenate = TRUE`, the function concatenate the test
#' statistics from results stored in the list and outputs a single file.
#'
#' @author I-Hsuan Lin
#'
#' @name exportResList
#'
#' @seealso NULL
#'
#' @export
#' @importFrom rlang inform
#' @importFrom rlang abort
#' @examples
#' library(SingleCellExperiment)
#'
#' # Load demo dataset
#' data(sce)
#'
#' # Use result saved after running findMarkers and 'pval.type = "any"'
#' # Save results into individual files in the per-session temporary directory
#' exportResList(metadata(sce)[["findMarkers1"]], dir_path = tempdir())
#'
#' # Save concatenate result in the per-session temporary directory
#' exportResList(metadata(sce)[["findMarkers4"]],
#'   concatenate = TRUE,
#'   prefix = "findMarkers4", dir_path = tempdir()
#' )
exportResList <- function(object, concatenate = FALSE, col_anno = NULL, direction = "both",
                          prefix = "output", dir_path = NULL) {
  type <- .get_res.type(object)
  .check_res.type(type, "is.list")
  inform(sprintf("Detecting %s input.", type))

  if (concatenate) {
    # Set up annotation columns
    first <- .res2df(object[[1]], type)

    if (!is.null(col_anno)) {
      # Check if all 'col_anno' columns can be found
      df <- if (all(col_anno %in% colnames(first))) first[, col_anno] else abort("Cannot find the specified columns in results.")
    } else {
      df <- data.frame(Symbol = rownames(first), row.names = rownames(first))
    }
    df <- df[order(rownames(df)), , drop = FALSE]
  }

  # Build the rest of the DataFrame
  for (group in names(object)) {
    x <- object[[group]]

    # If group is a number, add 'Cluster' before the number
    newg <- ifelse(grepl("^[0-9]+$", group), paste0("Cluster", group), group)

    if (concatenate) {
      test.cols <- .get_test.cols(x, type)
      test.cols <- test.cols[order(rownames(test.cols)), ]
      colnames(test.cols) <- paste0(newg, "_", colnames(test.cols)) # Append group/newg to columns
      df <- cbind(df, test.cols)
    } else {
      df <- .res2df(x, type)

      if (direction %in% c("up", "down")) {
        filename <- sprintf("%s_%s_%sregulated_%s.tsv", prefix, type, direction, as.character(newg))
      } else {
        filename <- sprintf("%s_%s_%s.tsv", prefix, type, as.character(newg))
      }

      inform(sprintf("Creating file: %s", filename))
      dir_path <- .build.path(dir_path)
      write.table(df,
        file = file.path(dir_path, filename), sep = "\t", quote = F,
        row.names = F, col.names = T
      )
    }
  }

  if (concatenate) {
    if (direction %in% c("up", "down")) {
      filename <- sprintf("%s_%s_%sregulated_concatenated.tsv", prefix, type, direction)
    } else {
      filename <- sprintf("%s_%s_concatenated.tsv", prefix, type)
    }

    inform(sprintf("Creating a concatenated file: %s", filename))
    dir_path <- .build.path(dir_path)
    write.table(df,
      file = file.path(dir_path, filename), sep = "\t", quote = F,
      row.names = F, col.names = T
    )
  }
}

#' Perform enrichR analysis on results from findMarkers, edgeR or DESeq2
#'
#' This function takes a list of objects containing `findMarkers()`,
#' \pkg{edgeR} or \pkg{DESeq2} results, and perform enrichment analysis on
#' the genes that passed the specified threshold.
#'
#' @param object A named list of objects storing outputs returned by
#' `findMarkers()`, \pkg{edgeR} or \pkg{DESeq2}, each of which contains
#' results for the corresponding group of comparison.
#' @param dbs A character vector containing the name of one or more gene-set
#' libraries used in the enrichment analysis.
#' @param site A string indicating the Enrichr Website to use. Available sites
#' are: `"Enrichr"`, `"FlyEnrichr"`, `"WormEnrichr"`, `"WormEnrichr"` and
#' `"FishEnrichr"`. Default is "Enrichr".
#' @param direction A string indicating the directionality of expression
#' changes of markers/DEGs to retreive. Allowable character values are
#' `"both"`, `"up"` and `"down"`. Default is "both".
#' @param fdr A numeric scalar indicating the FDR threashold used to select
#' markers/DEGs. This threshold has no effect when the input source is from
#' `findMarkers()` ran with `pval.type = "any"`. Default is 0.05.
#' @param top An integer scalar indicating the 'Top' threashold used to select
#' markers. This threshold is only applicable when the input source is from
#' `findMarkers()` ran with `pval.type = "any"`. Default is 200.
#' @param logfc A numeric vector of length 2 specifying the logFC threasholds
#' used to select up-/down-regulated markers/DEGs. This threshold has
#' no effect when the input source is from `findMarkers()` ran with
#' `test.type = "wilcox"`. Default is `c(0, 0)`.
#' @param auc A numeric vector of length 2 specifying the AUC threasholds
#' used to select up-/down-regulated markers. This threshold is only
#' applicable when the input source is from `findMarkers()` ran with
#' `test.type = "wilcox"`. Default is `c(0.3, 0.7)`.
#' @param max An integer scalar indicating the maximum number of markers/DEGs
#' to use in the analysis. Default is 1000.
#' @param min An integer scalar indicating the minimum number of markers/DEGs
#' required to perform the analysis. Default is 20.
#' @param column_by A string indicating the column name containing the gene
#' symbols or IDs in the input `object`'s DataFrames and the selected
#' information is returned. When `column_by = NULL`, it retreives the
#' row names from the DataFrame. Default is `NULL`.
#'
#' @return A named list of list of DataFrames
#'
#' @details This function uses `enrichr()` from the \pkg{enrichR} package
#' to submit gene symbols for enrichment analysis to the Enrichr website
#' (\url{https://maayanlab.cloud/Enrichr}). By default (with
#' `column_by = NULL`), this function extracts rownames from the input
#' DataFrames and submits the acquired strings for the analysis.
#'
#' The gene symbols are submitted to the main site (`site = "Enrichr"`)
#' with gene sets compiled from human and mouse genes. By change the
#' `site` setting, it can use other modEnrichr sites
#' (\url{https://maayanlab.cloud/modEnrichr/}) suitable for other model
#' organisms.
#'
#' @author I-Hsuan Lin
#'
#' @name runEnrichR
#'
#' @seealso [enrichR::enrichr()], [scran::findMarkers()], [edgeR::topTags()],
#' [DESeq2::results()]
#'
#' @export
#' @importFrom enrichR enrichr
#' @importFrom rlang inform
#' @importFrom utils flush.console
#' @examples
#' \dontrun{
#' # Load demo dataset
#' data(res_deseq2)
#'
#' # We construct a list of 3 DataFrames using the same DESeq2 output (`res_deseq2`)
#' # to mimic having a list containing results from multiple sets of comparisons.
#' res.de <- list(A_B = res_deseq2, A_C = res_deseq2, B_C = res_deseq2)
#'
#' # Select gene-set libraries
#' dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018")
#'
#' # Run enrichR using the D. melanogaster specific modEnrichr site, and
#' # specify the gene symbols are stored in the 'external_gene_name' column
#' res.ora <- runEnrichR(res.de,
#'   dbs = dbs, site = "FlyEnrichr",
#'   column_by = "external_gene_name"
#' )
#' }
runEnrichR <- function(object, dbs, site = "Enrichr", direction = "both", fdr = 0.05,
                       top = 200, logfc = c(0, 0), auc = c(0.7, 0.3), max = 1000,
                       min = 20, column_by = NULL) {
  type <- .get_res.type(object)
  .check_res.type(type, "is.list")
  .check_is.null(dbs)
  inform(sprintf("Detecting %s input.", type))

  dir.type <- if (direction == "both") "up- and down-" else paste0(direction, "-")
  dir.type <- paste0(dir.type, "regulated")

  resList <- list()
  genes <- NULL
  .resetEnrichrSite(site)
  for (group in names(object)) {
    x <- object[[group]]

    if (type == "findMarkers") {
      genes <- getMarkers(x,
        direction = direction, fdr = fdr, top = top, logfc = logfc,
        auc = auc, max = max, column_by = column_by
      )
    } else {
      genes <- getDEGs(x,
        direction = direction, fdr = fdr, logfc = logfc, max = max,
        column_by = column_by
      )
    }

    # If group is a number, add 'Cluster' before the number
    newg <- ifelse(grepl("^[0-9]+$", group), paste0("Cluster", group), group)

    # Run enrichr if there are at least 'min' genes satisfying the criteria; else skip
    if (length(genes) >= min) {
      inform(sprintf(
        "Running enrichR on '%s' with %d %s genes.", as.character(newg),
        length(genes), dir.type
      ))
      res <- enrichr(genes, databases = dbs)
      resList[[group]] <- res
    } else {
      inform(sprintf(
        "Skip '%s' with %d %s genes.", as.character(newg),
        length(genes), dir.type
      ))
    }
    if (interactive()) flush.console()
  }
  resList
}

#' Visualise Enrichr results as barplots
#'
#' This function takes a list of list of DataFrames containing Enrichr results
#' returned by `runEnrichR()` and create barplots from a selected gene-set
#' library.
#'
#' @param object A named list of list of DataFrames storing outputs returned
#' by `runEnrichR()`, each of which contains results for the corresponding
#' group of comparison.
#' @param db A string indicating the name of one gene-set library to plot.
#' @param showTerms An integer scalar indicating the number of terms to show.
#' Default is 20.
#' @param numChar An integer scalar indicating the number characters to keep
#' in the terms' descriptions. Default is 50.
#' @param y A string indicating the variable that should be mapped to the
#' y-axis. It can be `"Count"` and `"Ratio"`. Default is "Count".
#' @param order_by A string indicating how to order the Enrichr terms before
#' selecting the first `showTerms` terms to plot. It can be `"P.value"` or
#' `"Combined.Score"`. Default is "P.value".
#' @param theme_size A numeric scalar indicating the base font size.
#' Default is 16.
#' @param prefix A string indicating the prefix of output file. When
#' `prefix = NULL`, the plots are shown in the current graphics device using
#' [grid::grid.draw()]. The plots are saved to PDF files using `ggsave()`
#' when `prefix` is not `NULL`. Default is `NULL`.
#' @param dir_path The directory path for exported PDF. The files are exported
#' to the current working directory when `file_path = NULL`. Default is `NULL`.
#' @param ... Other arguments passed on to `ggsave()`.
#'
#' @return None
#'
#' @details The function prints barplots from each group of comparison using
#' [grid::grid.draw()] one after another, therefore the plotting feature is
#' suitable only in Jupyter Notebook or R Markdown. Otherwise, disable plotting
#' by specifying the `prefix` argument to save barplots to PDF files.
#'
#' @author I-Hsuan Lin
#'
#' @name plotEnrichR
#'
#' @seealso [runEnrichR()], [ggplot2::ggsave()]
#'
#' @export
#' @import ggplot2
#' @importFrom enrichR plotEnrich
#' @importFrom rlang inform
#' @importFrom rlang abort
#' @importFrom utils flush.console
#' @importFrom grid grid.draw
#' @examples
#' \dontrun{
#' # Load demo dataset
#' data(res_deseq2)
#'
#' # We construct a list of 3 DataFrames using the same DESeq2 output (`res_deseq2`)
#' # to mimic having a list containing results from multiple sets of comparisons.
#' res.de <- list(A_B = res_deseq2, A_C = res_deseq2, B_C = res_deseq2)
#'
#' # Select gene-set libraries
#' dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018")
#'
#' # Run enrichR using the D. melanogaster specific modEnrichr site, and
#' # specify the gene symbols are stored in the 'external_gene_name' column
#' res.ora <- runEnrichR(res.de,
#'   dbs = dbs, site = "FlyEnrichr",
#'   column_by = "external_gene_name"
#' )
#'
#' # Print plots on to standard output
#' plotEnrichR(res.ora, db = "GO_Biological_Process_2018", theme_size = 14)
#'
#' # Save plots to PDF files in the per-session temporary directory
#' plotEnrichR(res.ora,
#'   db = "GO_Biological_Process_2018", prefix = "Enrichr",
#'   dir_path = tempdir(), width = 12, height = 5
#' )
#' }
plotEnrichR <- function(object, db, showTerms = 20, numChar = 50, y = "Count",
                        order_by = "P.value", theme_size = 16, prefix = NULL,
                        dir_path = NULL, ...) {
  type <- .get_res.type(object)
  .check_res.type(type, "enrichr")
  .check_is.null(db)
  if (!is.null(prefix)) {
    dir_path <- .build.path(dir_path)
  }

  for (group in names(object)) {
    dbs <- names(object[[group]])
    if (!db %in% dbs) abort(sprintf("Gene-set library '%s' is not found in `object`.", db))

    # If group is a number, add 'Cluster' before the number
    newg <- ifelse(grepl("^[0-9]+$", group), paste0("Cluster", group), group)
    title <- sprintf("enrichR (%s) on '%s'", db, as.character(newg))

    p <- plotEnrich(object[[group]][[db]],
      showTerms = showTerms, numChar = numChar, y = y,
      orderBy = order_by, title = title
    ) + theme_bw(base_size = theme_size)

    if (is.null(prefix)) {
      # Return plot
      grid.draw(p)
      if (interactive()) flush.console()
    } else {
      filename <- sprintf("%s_%s_%s.pdf", prefix, as.character(newg), db)
      inform(sprintf("Creating file: %s", filename))
      # Save plot
      ggsave(file.path(dir_path, filename), plot = p, device = "pdf", ...)
    }
  }
}

#' Export results from Enrichr to tsv
#'
#' This function takes a list of list of DataFrames containing Enrichr results
#' returned by `runEnrichR()` and save all the results in the list object to
#' individual text files.
#'
#' @param object A named list of list of DataFrames storing outputs returned
#' by `runEnrichR()`, each of which contains results for the corresponding
#' group of comparison.
#' @param prefix A string indicating the prefix of output file.
#' Default is "enrichr".
#' @param showTerms An integer scalar indicating the number of terms to save
#' to file. All terms are saved when `showTerms = NULL`. Default is `NULL`.
#' @param columns An integer vector indicating the columns from each entry
#' of data to save to file. All columns are saved when `columns = c(1:9)`.
#' Default is `c(1:9)`.
#' 1-"Term", 2-"Overlap", 3-"P.value", 4-"Adjusted.P.value",
#' 5-"Old.P.value", 6-"Old.Adjusted.P.value", 7-"Odds.Ratio",
#' 8-"Combined.Score",  9-"Combined.Score".
#' @param dir_path The directory path for exported TSV. The files are exported
#' to the current working directory when `file_path = NULL`. Default is `NULL`.
#'
#' @return NULL
#'
#' @details We are not using `printEnrich()` in the \pkg{enrichR} package to
#' print results due to a bug not yet fixed in the current verion in CRAN at
#' the time of writing (v3.0).
#'
#' @author I-Hsuan Lin
#'
#' @name printEnrichR
#'
#' @seealso [runEnrichR()]
#'
#' @export
#' @importFrom rlang inform
#' @importFrom rlang abort
#' @importFrom utils write.table
#' @examples
#' \dontrun{
#' # Load demo dataset
#' data(res_deseq2)
#'
#' # We construct a list of 3 DataFrames using the same DESeq2 output (`res_deseq2`)
#' # to mimic having a list containing results from multiple sets of comparisons.
#' res.de <- list(A_B = res_deseq2, A_C = res_deseq2, B_C = res_deseq2)
#'
#' # Select gene-set libraries
#' dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018")
#'
#' # Run enrichR using the D. melanogaster specific modEnrichr site, and
#' # specify the gene symbols are stored in the 'external_gene_name' column
#' res.ora <- runEnrichR(res.de,
#'   dbs = dbs, site = "FlyEnrichr",
#'   column_by = "external_gene_name"
#' )
#'
#' # Save results to TSV files in the per-session temporary directory
#' printEnrichR(res.ora, prefix = "Enrichr", dir_path = tempdir())
#' }
printEnrichR <- function(object, prefix = "enrichr", showTerms = NULL,
                         columns = c(1:9), dir_path = NULL) {
  type <- .get_res.type(object)
  .check_res.type(type, "enrichr")
  .check_wholenum(columns)
  dir_path <- .build.path(dir_path)

  for (group in names(object)) {
    x <- object[[group]]

    # If group is a number, add 'Cluster' before the number
    newg <- ifelse(grepl("^[0-9]+$", group), paste0("Cluster", group), group)

    for (i in 1:length(x)) {
      db <- names(x)[i]
      df <- x[[i]]
      df <- .enrichment_prep_df(df, showTerms, orderBy = "P.value")
      df <- df[, !colnames(df) %in% c("Annotated", "Significant")]

      if (any(columns > ncol(df))) abort("Undefined columns selected.")

      filename <- sprintf("%s_%s_%s.tsv", prefix, as.character(newg), db)
      inform(sprintf("Creating file: %s", filename))
      write.table(df,
        file = file.path(dir_path, filename), sep = "\t", quote = F,
        row.names = F, col.names = T
      )
    }
  }
}
