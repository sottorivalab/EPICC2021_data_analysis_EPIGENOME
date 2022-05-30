# This script gets the per-sample read counts for the recurrent peaks
library(GenomicRanges)
source("setup_environment/0-source.R")
source("setup_environment/1-dirs_and_files.R")
source("analysis/ATAC/CLAIRE/Epigenetic/functions/horizon_plot_functions_commonE_CL6.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"
regex_peak_calls = "_nucleosome_free-region_[A-Z]-GRCh38-filtered_peaks[.]bed"
  
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

peaks = readRDS(file.path(dat_dir, "recurrent_peaks.rds"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

load_peaks = function(x, peaks) {
  d = pbapply::pblapply(x, read_filtered_peak_file) %>% 
    lapply(makeGRangesFromDataFrame, keep.extra.columns=TRUE) %>% 
    lapply(subsetByOverlaps, peaks)
  
  str_mt = stringr::str_match(basename(x) ,"(.*?)_.*region_(.*?)-GRCh38.*")
  names(d) = apply(str_mt[,c(2,3)], 1, paste, collapse = "_")
  
  return(d)
}

create_track = function(d) {
  tracks = list()
  for (i in seq_along(d)) {
    region = sub(".*_", "", names(d)[i])
    color = col[region]
    tracks[[names(d)[i]]] = Gviz::AnnotationTrack(d[[i]], name=NULL, fill=color)
  }
  return(tracks)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


bed_dir = file.path(.pipeline_ddir, "atac", "peaks")
beds = list.files(bed_dir, regex_peak_calls, full.names = TRUE,  recursive=TRUE)
beds_cancer = beds[grepl("_A|_B|_C|_D|_F|_G|_H|_I||pan", beds)]

# load peak beds
beds_gr = load_peaks(beds_cancer, peaks)
saveRDS(beds_gr, file.path(dat_dir, "beds.gr.rds"))

# subset to peaks overlapping peaks of region
peak_matches_df = do.call(what=cbind, lapply(beds_gr, overlapsAny, query = peaks))
dimnames(peak_matches_df) = list(peaks$peaks_to_plot, names(beds_gr))
saveRDS(peak_matches_df, file.path(dat_dir, "peak_matches_df.rds"))

# create peak tracks
peak_tracks = create_track(beds_gr)
saveRDS(peak_tracks, file.path(dat_dir, "peak_tracks.rds"))
