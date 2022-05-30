# This script gets the per-sample read counts for the recurrent peaks
library(GenomicRanges)
library(stringr)
library(optparse)
source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"
regex_cutsites = ".*[1-9]_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Parse arguments ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

option_list = list(
  make_option(
    c("-p", "--patient"),
    type = "character",
    default = NULL,
    metavar = "character"
  ),
  make_option(
    c("-r", "--regions"),
    type = "character",
    default = NULL,
    metavar = "character"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

patient = opt$patient
peaks = readRDS(opt$regions)
peak_file = sub(".rds","", basename(opt$regions))
out_file = file.path(dat_dir, "readcounts", paste(patient,  peak_file, "nf_counts.rds", sep="_"))
if (file.exists(out_file)) quit(save = "no")
dir.create(dirname(out_file), FALSE, TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

load_counts = function(x, peaks) {
  gr = read_shifted_count_file(x) %>% resize(width = 100, fix = "center")
  countOverlaps(peaks, gr)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

excluded_samples = na.omit(str_match(.excluded_samples, "\\s*(.*?)\\s*_C1$")[,2])

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# find input files
nf_cutsites = .atac_bed_files
names(nf_cutsites) = gsub("_C[0-9]+$", "", names(.atac_bed_files))

# drop excluded samples
exclude = names(nf_cutsites) %in% excluded_samples | !grepl(patient, names(nf_cutsites))
nf_cutsites = nf_cutsites[!exclude]

# load per sample data
print(nf_cutsites)
peak_counts = pbapply::pblapply(nf_cutsites, load_counts, peaks=peaks)
values(peaks) = peak_counts
colnames(values(peaks)) = gsub("^EPICC_", "", names(nf_cutsites))

# save data
saveRDS(peaks, file = out_file)

