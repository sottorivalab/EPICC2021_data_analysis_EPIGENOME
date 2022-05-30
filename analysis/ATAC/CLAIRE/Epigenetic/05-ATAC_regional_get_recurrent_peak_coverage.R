# This script gets the per-sample read counts for the recurrent peaks
library(GenomicFeatures)
library(GenomicRanges)
library(stringr)
source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

idat_dir = "analysis/ATAC/CLAIRE/Epigenetic/input_data"
dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"
regex_cutsites = ".*[1-9]_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"
idx = as.numeric(commandArgs(trailingOnly=TRUE)[1])

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

read_cuts_in_interval = function(f, gr) {
  rtracklayer::import(f, format = "BED", genome = "hg38", which = gr) %>%
    resize(width = 100, fix = "center")
}

read_concated_cuts_intervals = function(f, gr) {
  lapply(f, read_cuts_in_interval, gr = gr) %>% 
    GRangesList %>% unlist
}

calculate_normalised_coverage = function(d, nf, bins) {
  cov = GenomicRanges::coverage(d) / nf
  cov = cov[seqlevels(cov) %in% seqlevels(bins)]
  cov = GenomicRanges::binnedAverage(bins, cov, "score")
  return(cov)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

excluded_samples = na.omit(str_match(.excluded_samples, "\\s*(.*?)\\s*_C1$")[,2])

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

peaks = readRDS(file.path(dat_dir, "recurrent_peaks.rds"))
purity_data = readRDS("created_datasets//genotyping_estimates_per_sample.rds")
coldata = readRDS(file.path(idat_dir, "coldata.rds"))
coldata$purity = purity_data[paste0("EPICC_", rownames(coldata), "_C1"),"estimated_purity"]
coldata$purity_short = round(coldata$purity*100, digits=0)
coldata = coldata[which(coldata$purity > 0),]
normdata = readRDS("created_datasets//peak_coverage_data_nucleosome_free_filtered.rds")
coldata_region = readRDS(file.path(idat_dir, "coldata_region.rds"))
pat_ids = unlist(readr::read_csv("analysis/ATAC/CLAIRE/Epigenetic/input_data/patients.txt", col_names = FALSE))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

norm_files = list.files(file.path(dat_dir, "normfacs"), full.names = TRUE, recursive = TRUE)
norm_facs = lapply(norm_files, readRDS) %>% unlist()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# find input files
nf_cutsites = .atac_bed_files
names(nf_cutsites) = gsub("_C[0-9]+$", "", names(.atac_bed_files))
  
# drop excluded samples
exclude = names(nf_cutsites) %in% excluded_samples
nf_cutsites = nf_cutsites[!exclude]

# drop excluded samples
annot = THmisc::annotation_from_barcode(paste0(names(nf_cutsites), "_C1"))
exclude = !annot$patient %in% pat_ids | !names(nf_cutsites) %in% paste0("EPICC_", rownames(coldata))
nf_cutsites = nf_cutsites[!exclude]

# get recurrent peaks as GRanges
resizePeaks = resize(peaks, width = 10000, fix = "center")
bins = unlist(tile(resizePeaks, width = 1))
combos = rownames(coldata_region)

# group cut site files
gr = substr(names(nf_cutsites), 7, 12)
gr_cut_files = split(nf_cutsites, gr)
wh_norm = substr(gr, 6, 6) == "E"
gr_cut_files = c(gr_cut_files, list(panpatientE = nf_cutsites[wh_norm]))

if (!is.na(idx)) {
  if (idx > length(gr_cut_files)) quit("n")
  gr_cut_files = gr_cut_files[idx]
}

#
of1 = file.path(dat_dir, "cut_data", paste0(names(gr_cut_files), ".rds"))
of2 = file.path(dat_dir, "cpm_data", paste0(names(gr_cut_files), ".rds"))
fexist = file.exists(of1) & file.exists(of2)
gr_cut_files = gr_cut_files[!fexist]
if (length(gr_cut_files) == 0) quit("no")


# read in the cutsites
cat("Reading cut sites\n")
cut_data = pbapply::pblapply(gr_cut_files, read_concated_cuts_intervals, gr=resizePeaks)
norm_factor = lapply(gr_cut_files, function(x) mean(coldata[names(x), "libsize"] * norm_facs[names(x)] / 1e6))
cpm_data = mapply(calculate_normalised_coverage, cut_data, norm_factor, list(bins))          


# save data
for (i in seq_along(norm_factor)) {
  of = file.path(dat_dir, "cut_data", paste0(names(cut_data)[i], ".rds"))
  dir.create(dirname(of), FALSE, TRUE)
  saveRDS(cut_data[[i]], of)
  of = file.path(dat_dir, "cpm_data", paste0(names(cpm_data)[i], ".rds"))
  dir.create(dirname(of), FALSE, TRUE)
  saveRDS(cpm_data[[i]], of)
}
