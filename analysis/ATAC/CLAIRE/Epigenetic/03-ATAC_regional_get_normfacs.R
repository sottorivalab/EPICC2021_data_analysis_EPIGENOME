# This script gets the normalisation factors for DESEQ (using counts in all peaks)
source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir = "normfacs"
idat_dir = "analysis/ATAC/CLAIRE/Epigenetic/input_data"
dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"
peak_file = "all_peaks"
  
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Parse arguments ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Get arguments and coldata
args = commandArgs(trailingOnly=TRUE)
patient_ = args[1]
out_file = file.path(dat_dir, dir, paste0(patient_, "_nf_counts_patientSpecPeakSizeFactors.rds"))
if (file.exists(out_file)) quit(save="n")
dir.create(dirname(out_file), FALSE, TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

purity_data = readRDS("created_datasets/genotyping_estimates_per_sample.rds")
counts = readRDS(file.path(dat_dir, "readcounts", paste(patient_, peak_file, "nf_counts.rds", sep="_")))
coldata = readRDS(file.path(idat_dir, "coldata.rds"))
coldata$purity = purity_data[paste0("EPICC_", rownames(coldata), "_C1"),"estimated_purity"]
coldata$purity_short <- round(coldata$purity*100, digits=0)
normdata = readRDS("created_datasets/peak_coverage_data_nucleosome_free_filtered.rds")
coldata$subtissue[coldata$subtissue == "normal"] = "cancer"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Remove samples with no purity info/ adenomas from coldata
wh_keep = with(coldata, !is.na(purity) & purity > 0 & patient == patient_)
coldata = coldata[wh_keep, ]

S4Vectors::mcols(counts) = S4Vectors::mcols(counts)[,colnames(S4Vectors::mcols(counts)) %in% rownames(coldata)]
coldata = coldata[colnames(S4Vectors::mcols(counts)),]

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# get counts, just for sample info
all_samp = names(S4Vectors::mcols(counts))
model = "~purity+sample_type+subtissue"
if (length(unique(coldata$subtissue)) < 2) model = gsub("[+]subtissue", "", model)
model = as.formula(model)
              
# subset to counts detected in case
regex = paste0(patient_, ".region_[ABCDFGHI]")
count_matrix = normdata$data[grepl(regex, normdata$peak_calls$name),]
annot = THmisc::annotation_from_barcode(colnames(count_matrix), TRUE)
colnames(count_matrix) = gsub("EPICC_", "", annot$tissue_barcode)
count_matrix = count_matrix[,colnames(count_matrix) %in% rownames(coldata)]

# Get normfacs from  peak set
coldata = coldata[colnames(count_matrix),]
dds = DESeq2::DESeqDataSetFromMatrix(count_matrix, colData=coldata, design=model)
peakSizeFactors = DESeq2::estimateSizeFactors(dds)$sizeFactor
saveRDS(peakSizeFactors, out_file)

