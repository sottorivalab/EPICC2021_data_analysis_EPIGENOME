# This script gets the per-sample read counts for the recurrent peaks
library(GenomicFeatures)
library(GenomicRanges)
library(stringr)
source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"
regex_cutsites = ".*[1-9]_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

rec_data = readRDS("analysis/ATAC/MEGABULKS/datasets/edger_summary/final_reccurence_data.rds")
pat_ids = unique(gsub("_.*", "", list.files(file.path(dat_dir, "fits"))))
peaks_original = readRDS(file.path(dat_dir, "recurrent_peaks.rds"))
peak_matches_df = readRDS(file.path(dat_dir, "peak_matches_df.rds"))
peak_matches_df = peak_matches_df[,colnames(peak_matches_df) != "pan_E"]
peak_matches_df = peak_matches_df[rownames(peak_matches_df) %in% peaks_original$peaks_to_plot,]
annot_peak_matches = THmisc::annotation_from_barcode_epicc(paste0("EPICC_", colnames(peak_matches_df), "1_G1_C1"))
 
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

peaks = resize(peaks_original, width = 2000, fix = "center")
peaks_df = mcols(peaks)
colnames(peaks_df) = "peak"

sig_filter = rec_data$sig_matrix[peaks_original$peaks_to_plot,]
fc = rec_data$fc[peaks_original$peaks_to_plot,]
p = rec_data$p[peaks_original$peaks_to_plot,]
p_adj = rec_data$p_adj[peaks_original$peaks_to_plot,]

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


parse_file_set = function(x, regex_keep=NULL, suffix="") {
  
  if (!is.null(regex_keep)) 
    x = x[grepl(regex_keep, x)]
  
  # load data 
  res = list()
  for (i in seq_along(x)) {
    d = readRDS(x[i])
    res[[i]] = data.frame(peak=row.names(d), pvalue=d$pvalue, padj=d$padj)
  }
  
  # append ids to colnames
  names(res) = gsub("_.*", "", basename(x))
  for (i in seq_along(res)) {
    colnames(res[[i]])[-1] = paste0(names(res)[i], "_", suffix, colnames(res[[i]])[-1])
  }
  
  # reduce data
  Reduce(function(x, y) merge(x, y, by="peak", all=T), res)
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

regex_fits = c(
  "purity_only" = "_sample_type[+]purity_deseq_LRT[.]rds",
  "region_only" = "_sample_type[+]region_deseq_LRT[.]rds",
  "region" = "_purity[+]sample_type[+]region_deseq_LRT[.]rds",
  "region" = "_purity[+]sample_type[+]subtissue_deseq_LRT[.]rds"
)

# Get DESEQ results
regex_pats = paste(pat_ids, collapse="|")
fit_dir = file.path(dat_dir, "fits")

result_sets = list()
for (i in seq_along(regex_fits)) {
  fits = list.files(fit_dir, regex_fits[i], full.names = TRUE)
  d = parse_file_set(fits, regex_pats, paste0(names(regex_fits)[i], "_nf_"))
  
  if (!is.null(result_sets[[names(regex_fits)[i]]])) {
    result_sets[[names(regex_fits)[i]]] = cbind(result_sets[[names(regex_fits)[i]]], d[,colnames(d) != "peak"])
  } else {
    result_sets[[names(regex_fits)[i]]] = d
  }
}

# reduce and safe
merge_by_peak = function(x, y) merge(x, y, by="peak", all=T)
results = Reduce(merge_by_peak, c(list(peaks_df), result_sets))
saveRDS(results, file.path(dat_dir, "all_deseq_results.rds"))


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#MAKE SUBCLONAL MATRIX FOR HEATMAP
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

sc_res = results[,paste0(pat_ids,"_region_nf_padj")] %>%
  mapply(FUN=as.numeric) %>% 
  magrittr::set_rownames(results$peak) %>% 
  magrittr::set_colnames(pat_ids)

mod_cn = gsub("(-$)|[()]", "", gsub("(-normal_vs_(pure)?)| ", "-", colnames(fc)))
fcc = fc[rownames(sc_res), match(colnames(sc_res), mod_cn)] 
fcc[is.na(fcc)] = 0

annot_pm = gsub("[(]|[)]|(-cancer)", "", gsub(" |[.]", "-", annot_peak_matches$tumour_id))
peak_mt = t(apply(peak_matches_df, 1, tapply, annot_pm, any))[rownames(sc_res), colnames(sc_res)]

mt = match(rownames(sc_res), rec_data$summary$peak)
wh_gain = rec_data$summary$event_type[mt] == "gain"
wh_loss = rec_data$summary$event_type[mt] == "loss"

subclonal_matrix = matrix(FALSE, NROW(sc_res), NCOL(sc_res), dimnames = dimnames(sc_res))
subclonal_matrix[wh_loss,] = sc_res[wh_loss,] < 0.05 & fcc[wh_loss,] < 0
subclonal_matrix[wh_gain,] = sc_res[wh_gain,] < 0.05 & fcc[wh_gain,] > 0 & peak_mt[wh_gain,]
subclonal_matrix[is.na(subclonal_matrix)] = FALSE

saveRDS(subclonal_matrix, file.path(dat_dir, "subclonal_matrix.rds")) 
