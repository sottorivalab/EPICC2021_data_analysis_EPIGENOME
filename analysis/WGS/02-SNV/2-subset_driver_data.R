source("setup_environment/0-source.R")
library(cowplot)
library(ggplot2)
library(VariantAnnotation)
library(dplyr)
library(reshape2)
library(seqinr)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options  ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

vcf_dir = "analysis/WGS/02-SNV/datasets/vcf_datasets"
dat_dir = "analysis/WGS/02-SNV/datasets"

gene_lists =
  c(
    "IntOGen-DriverGenes_COREAD",
    "coad_driver_genes_tcga_2012",
    "CRC_drivers_Cross2018_Tier1",
    "driver_genes_martincorena_2017",
    "driver_genes_tarabichi_2018"
  )

genes_from_lists = THmisc::gene_lists[gene_lists]

genes_from_lists[["chromatin_modifers_selected"]] = 
  grep(
    c("^(KMT.*)|(KDM.*)|(TET.*)|(DOT.*)|(KAT.*)|(DNMT.*)"), 
    THmisc::gene_lists$all_genes, 
    value = TRUE
  )

genes_from_lists$chromatin_modifer_genes = 
  "external_datasets/chromatin_modifer_genes.txt" %>% 
  readr::read_csv(col_names = FALSE) %>% 
  unlist() %>% magrittr::set_names(NULL)

drivers_to_export = c(
  THmisc::gene_lists$`IntOGen-DriverGenes_COREAD`
)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)
  
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  
.is_only_driver = function(x) {
  all(unlist(strsplit(x, ";")) %in% unlist(genes_from_lists))
}

.is_driver = function(x)  {
  any(unlist(strsplit(x, "[:,&; ]")) %in% drivers_to_export)
}

.load_vcf = function(x, ...) {
  print(x)
  readRDS(x) %>% # only keep genes in gene_lists to reduce size
    THmisc::vcf_to_data_frame(..., include_annot=TRUE) %>%
    dplyr::filter(!(variant == "" | is.na(variant))) %>% # drop irrelevant annotations (e.g., low or no impact)
    cbind(THmisc::annotation_from_barcode(.$sample), .) # append sample annotations 
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

vcf_driver_data_frames = # Find and load vcfs, might take a while...
  list.files(vcf_dir, "[.]rds$", full.names = TRUE) %>% 
  magrittr::set_names(gsub("[.]rds", "", basename(.))) %>%
  pbapply::pblapply(.load_vcf, genes=unlist(genes_from_lists))

# assert that only correct genes included
d_all = do.call(what=rbind, vcf_driver_data_frames)
print(d_all[!sapply(d_all$gene, .is_only_driver),])
stopifnot(all(sapply(d_all$gene, .is_only_driver)))

saveRDS(vcf_driver_data_frames, file.path(dat_dir, "driver_gene_data.rds"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

vcf_driver_data_frames %>% # used as a look-up table for labels...
  do.call(what=rbind) %>%
  dplyr::filter(sapply(gene, .is_driver)) %>% 
  dplyr::select(variant=variant, id=mutation) %>%
  unique() %>%
  magrittr::set_rownames(NULL) %>%
  saveRDS(file=file.path(dat_dir, "all_drivers_list.rds"))

