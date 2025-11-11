#' Encode DNA sequences into nucleotide position features
#'
#' @description
#' Converts DNA 5-mer sequences (e.g., "ATCGA") into a structured data frame
#' with separate columns for each nucleotide position.
#'
#' @param seqs A character vector containing DNA 5-mer sequences (e.g., "ATCGA", "TTTGG").
#' @return A data frame with columns nt_pos1–nt_pos5.
#' @examples
#' dna_encoding(c("ATCGA", "TTTGG"))
#' @export
dna_encoding <- function(seqs) {
  if (!is.character(seqs)) stop("Input must be a character vector.")
  seqs <- toupper(seqs)
  if (any(nchar(seqs) != 5)) stop("All sequences must be 5 nucleotides long.")
  seq_split <- strsplit(seqs, "")
  df <- do.call(rbind, lapply(seq_split, function(x) {
    data.frame(nt_pos1 = x[1], nt_pos2 = x[2], nt_pos3 = x[3],
               nt_pos4 = x[4], nt_pos5 = x[5], stringsAsFactors = FALSE)
  }))
  df[] <- lapply(df, factor, levels = c("A", "T", "C", "G"))
  return(df)
}

#' Predict m6A modification for multiple candidate sites
#'
#' @description Predicts m6A modification probabilities and statuses
#' for multiple RNA sites using a trained random forest model.
#' @param ml_fit A fitted random forest model.
#' @param feature_df A data frame with required feature columns.
#' @param positive_threshold Threshold for classifying positive sites (default 0.5).
#' @return Data frame with predicted probabilities and statuses.
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' head(prediction_multiple(ml_fit, feature_df, positive_threshold = 0.6))
#' @import randomForest
#' @importFrom stats predict
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  required <- c("gc_content","RNA_type","RNA_region","exon_length",
                "distance_to_junction","evolutionary_conservation","DNA_5mer")
  stopifnot(all(required %in% names(feature_df)))

  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA","lincRNA","lncRNA","pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS","intron","3'UTR","5'UTR"))
  enc <- dna_encoding(feature_df$DNA_5mer)
  newdata <- cbind(feature_df, enc)
  needed <- c("gc_content","RNA_type","RNA_region","exon_length",
              "distance_to_junction","evolutionary_conservation",paste0("nt_pos",1:5))

  if (!is.null(ml_fit$terms)) {
    prob <- predict(ml_fit, newdata = newdata, type = "prob")
  } else {
    prob <- predict(ml_fit, newdata = newdata[, needed], type = "prob")
  }
  pos_prob <- prob[, "Positive"]
  feature_df$predicted_m6A_prob <- pos_prob
  feature_df$predicted_m6A_status <- ifelse(pos_prob > positive_threshold, "Positive", "Negative")
  return(feature_df)
}

#' Predict m6A modification for a single candidate site
#'
#' @description Predicts m6A probability and classification for one site.
#'
#' @param ml_fit A fitted random forest model.
#' @param gc_content Numeric GC content.
#' @param RNA_type Character; one of \code{"mRNA"}, \code{"lincRNA"}, \code{"lncRNA"}, \code{"pseudogene"}.
#' @param RNA_region Character; one of \code{"CDS"}, \code{"intron"}, \code{"3'UTR"}, \code{"5'UTR"}.
#' @param exon_length Numeric exon length.
#' @param distance_to_junction Numeric distance to nearest exon junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer Character 5-mer DNA sequence.
#' @param positive_threshold Numeric; probability threshold for labeling (default 0.5).
#'
#' @return Named vector with:
#' \itemize{
#'   \item \code{predicted_m6A_prob}: numeric probability for "Positive".
#'   \item \code{predicted_m6A_status}: "Positive" or "Negative".
#' }
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(ml_fit, 0.6, "mRNA", "CDS", 1500, 120, 0.8, "ATGCA", 0.5)
#'
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  df <- data.frame(gc_content, RNA_type, RNA_region, exon_length,
                   distance_to_junction, evolutionary_conservation, DNA_5mer)
  res <- prediction_multiple(ml_fit, df, positive_threshold)
  return(c(predicted_m6A_prob = res$predicted_m6A_prob[1],
           predicted_m6A_status = res$predicted_m6A_status[1]))
}
