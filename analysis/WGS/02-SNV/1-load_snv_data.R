source("setup_environment/0-source.R")
library(VariantAnnotation)
library(THmisc)
library(BSgenome.Hsapiens.UCSC.hg38)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options  ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

vcf_dir = "analysis/WGS/02-SNV/datasets/vcf_datasets"
idx = suppressWarnings(as.numeric(commandArgs(trailingOnly=TRUE)[1]))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(vcf_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Loading of VCF files as R objects ####

all_vcfs = c(.vcfs_uf, .vcfs)
if (!is.na(idx)) all_vcfs = all_vcfs[idx]

for (vcf_file in all_vcfs) {
  out_file = file.path(vcf_dir, gsub("[.]vcf.*", ".rds", basename(vcf_file)))
  mtime = file.mtime(vcf_file)
  print(vcf_file)
  
  if (!file.exists(out_file) | isTRUE(mtime > file.mtime(out_file))) {
    data = epicc_vcf_loading_wrapper(vcf_file)
    saveRDS(data, out_file, version = 2)
    Sys.setFileTime(out_file, mtime)
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# md5sum all files in the vcf dir.  ####

if (is.na(idx)) {
  vcf_dir %>% 
    list.files("[.]rds$", full.names = TRUE) %>% 
    data.frame(file=.) %>% 
    dplyr::mutate(md5sum=tools::md5sum(file)) %>% 
    dplyr::mutate(file=basename(file)) %>%
    readr::write_tsv(file.path(vcf_dir, "md5sums.tsv"), col_names=FALSE)
}
