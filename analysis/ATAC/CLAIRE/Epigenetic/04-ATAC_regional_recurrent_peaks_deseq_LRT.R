# This script gets the normalisation factors for DESEQ (using counts in all peaks)
library(DESeq2)
source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir = "normfacs"
dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"
idat_dir = "analysis/ATAC/CLAIRE/Epigenetic/input_data"
peak_file = "recurrent_peaks"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Parse arguments ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Get arguments and coldata
args = c("C561", "~sample_type+region", "~sample_type", "effect_of_region_only_nfreads_LRT")
args = commandArgs(trailingOnly=TRUE)
print(args)
patient_ = sub("-.*", "", args[1])
subtissue_ = sub("^[a-zA-Z0-9]*(-)?", "", args[1])
out_prefix = paste0(patient_, ifelse(subtissue_ == "", "", paste0("-", gsub("_", "-", subtissue_))))
#if (subtissue_ == "") subtissue_ = "cancer"
model = as.formula(args[2])
reduced = as.formula(args[3])
dir = args[4] 

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

count_file = file.path(dat_dir, "readcounts", paste(patient_, peak_file, "nf_counts.rds", sep="_"))
mname = gsub(" ", "", as.character(model)[2])
fprefix = gsub("^C[0-9]+_", "", sub("[.]rds$", "", basename(count_file)))
out_file = file.path(dat_dir, "fits", paste0(out_prefix, "_", fprefix, "_", mname, "_deseq_LRT.rds"))
dir.create(dirname(out_file), FALSE, TRUE)
if (file.exists(out_file)) quit(save="n")

purity_data = readRDS("created_datasets/genotyping_estimates_per_sample.rds")
counts = readRDS(count_file)
coldata = readRDS(file.path(idat_dir, "coldata.rds"))
coldata$purity = purity_data[paste0("EPICC_", rownames(coldata), "_C1"),"estimated_purity"]
coldata$purity_short = round(coldata$purity*100, digits=0)
peaks_filt = peaks = readRDS(file.path(dat_dir, "recurrent_peaks.rds"))
nfpeakSizeFactors = readRDS(file.path(dat_dir, "normfacs", paste0(patient_, "_nf_counts_patientSpecPeakSizeFactors.rds")))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Remove samples with no purity info/ adenomas from coldata
wh_keep = with(coldata, !is.na(purity) & purity > 0 & patient == patient_ & (subtissue %in% c(subtissue_, "normal", "cancer")))
coldata = coldata[wh_keep, ]

S4Vectors::mcols(counts) = S4Vectors::mcols(counts)[,colnames(S4Vectors::mcols(counts)) %in% rownames(coldata)]
coldata = coldata[colnames(S4Vectors::mcols(counts)),]
coldata$subtissue[coldata$subtissue == "normal"] = "cancer"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Make the DESeq2 objects
count_matrix = as.matrix(mcols(counts))
rownames(count_matrix) = paste0(seqnames(counts), ":", start(counts), "-",end(counts))
regions = as.character(unique(coldata$region))
count_matrix = unique(count_matrix)

# Select only recurrent peaks
dds_fil = 
  DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix[rownames(count_matrix) %in% peaks$peaks_to_plot,],
    colData = coldata,
    design = model
  )


# Do DESEQ
dds_fil$sizeFactor = as.numeric(nfpeakSizeFactors[colnames(dds_fil)])
dds_fil = DESeq2::estimateDispersions(dds_fil)
dds_fil = DESeq2::DESeq(dds_fil, test="LRT", reduced=reduced)
res = DESeq2::results(dds_fil)

saveRDS(res, out_file)
