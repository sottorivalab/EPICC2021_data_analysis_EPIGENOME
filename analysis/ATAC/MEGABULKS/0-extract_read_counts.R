# This script gets the per-sample read counts for the recurrent peaks
library(GenomicRanges)
library(stringr)
source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dat_dir = "analysis/ATAC/MEGABULKS/datasets/counts"
idx = suppressWarnings(as.numeric(commandArgs(trailingOnly=TRUE)[1]))
peak_file = "created_datasets/peak_set.csv.gz"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

load_counts = function(x, peaks) {
  d = read_shifted_count_file(x)
  d_i = GRanges(seqnames(d), IRanges(start(d), start(d)))
  end(peaks) = end(peaks) - 1
  countOverlaps(peaks, d_i)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 
dir.create(dat_dir, FALSE, TRUE)

peaks = 
  readr::read_csv(peak_file)$peak %>% 
  strsplit(split = "[:-]") %>% 
  do.call(what =  rbind) %>% 
  data.frame() %>%  
  magrittr::set_colnames(c("chr", "start", "end")) %>% 
  as("GRanges")


# find input files
nf_cutsites = .atac_bed_files[!names(.atac_bed_files) %in% .excluded_samples]
if (!is.na(idx)) nf_cutsites = nf_cutsites[(seq_len(1000) %% 499) + 1 == idx]


# read counts form each file 
for (i in seq_along(nf_cutsites)) {
  ofile = file.path(dat_dir, paste0(names(nf_cutsites)[i], ".csv.gz"))
  if (file.exists(ofile)) next()
  peaks$n = load_counts(nf_cutsites[i], peaks)
  d = as.data.frame(peaks)[,c(1:3,6)]
  readr::write_csv(d, ofile)
}

# tabulate all counts 
if (is.na(idx)) {
  ifiles = list.files(dat_dir, "[.]csv[.]gz", full.names=TRUE)
  data = lapply(ifiles, function(x) readr::read_csv(x, show_col_types=FALSE, progress=FALSE)[,"n"])
  sid = sub("[.].*", "", basename(ifiles))
  
  d_mat = do.call(cbind, data)
  colnames(d_mat) = paste0(substr(sid, 7, 10), "_nucleosome_free-", sid)
  rownames(d_mat) = as.character(peaks)
  
  ofile = file.path(dat_dir, "data_matrix.rds")
  saveRDS(d_mat, ofile)
}
