source("setup_environment/0-source.R")
library(VariantAnnotation)
library(dndscv)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

idx = suppressWarnings(as.numeric(commandArgs(trailingOnly=TRUE)[1]))
vcf_dir = "analysis/WGS/02-SNV/datasets/vcf_datasets"
dat_dir = "analysis/WGS/02-SNV/datasets/dnds"
fig_dir = "analysis/WGS/02-SNV/plots/dnds"
ccf_cutoff = 0.25

gene_lists = c(
  "all_genes",
  "IntOGen-DriverGenes_COREAD",
  "IntOGen-DriverGenes_COREAD_sc_drivers_rm",
  "CRC_drivers_Cross2018_Tier1",
  "coad_driver_genes_tcga_2012",
  "driver_genes_martincorena_2017",
  "driver_genes_tarabichi_2018",
  "chromatin_modifers_selected",
  "chromatin_modifer_genes"
)

inference_neutral = 
  c("C528", "C530", "C532", "C536",
    "C537", "C543", "C544", "C548", 
    "C550", "C552", "C554", "C555", 
    "C560", "C561")

inference_selection = 
  c("C516", "C524", "C525", "C531",
    "C538", "C539", "C542", "C549", 
    "C551", "C559", "C562","C518")


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Modify gene lists ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

gene_lists = 
  THmisc::gene_lists[gene_lists] %>% 
  magrittr::set_names(gene_lists)

gene_lists[["chromatin_modifers_selected"]] = 
  sapply(c("^KMT.*", "^KDM.*", "^TET.*", "^DOT.*", "^KAT.*","DNMT.*"),
         grep, gene_lists$all_genes, value=TRUE) %>% unlist() %>% as.character()

gene_lists$`IntOGen-DriverGenes_COREAD` = 
  gene_lists$`IntOGen-DriverGenes_COREAD` %>% 
  (function(x) x[!x %in% c("PARP4", "LRP1B", "KMT2C")])

gene_lists$`IntOGen-DriverGenes_COREAD_sc_drivers_rm` = 
  gene_lists$`IntOGen-DriverGenes_COREAD` %>% 
  (function(x) x[!x %in% c("KRAS","PIK3CA","SMAD4","ARID1A","APC","BIRC6")])

gene_lists$chromatin_modifer_genes = 
  "external_datasets/chromatin_modifer_genes.txt" %>% 
  readr::read_csv(col_names = FALSE) %>% 
  unlist()

gene_lists$chromatin_modifer_genes_reactome = 
  "external_datasets/reactome_cmg_genes.tsv" %>% 
  readr::read_tsv(col_names = FALSE) %>% 
  dplyr::select(X3) %>% 
  unlist() %>% 
  gsub(pattern = ".* ", replacement = "")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# prepare output dir  ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load and lift data to hg19 ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Lifting data.\n")
vcf_files = list.files(vcf_dir, "[.]rds$", full.names=TRUE, recursive = TRUE)
if (!is.na(idx)) vcf_files = vcf_files[idx]
stopifnot(!any(duplicated(basename(vcf_files))))
for (i in seq_along(vcf_files)) {
  # print input file id   
  file_id = gsub("[.]rds", "", basename(vcf_files[i]))
  print(file_id)
  
  # create/check output dir
  out_dir = file.path(dat_dir, "lifted_vcfs", file_id)
  if (dir.exists(out_dir)) next()
  dir.create(out_dir, FALSE, TRUE)
  
  # lift and save data
  d_lifted = parse_rds_file_for_dndscv(vcf_files[i])
  for (tid in unique(d_lifted$tissue)) {
    for (mstat in unique(d_lifted$status)) {
      d_lifted %>% 
        dplyr::filter(tissue == tid & status == mstat) %>% 
        saveRDS(file.path(out_dir, paste0(tid, ".", mstat, ".rds")))
    }
  }
}

if (!is.na(idx)) 
  quit("no")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Grouping of mutation data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

get_annot = function(x) {
  strsplit(basename(x), "[.]") %>%
    do.call(what = rbind) %>%
    data.frame() %>%
    magrittr::set_colnames(c("case", "tissue", "status", "suffix")) %>%
    dplyr::mutate(is_filt = grepl("_filt", x)) %>% 
    # always pool adenoma variants
    dplyr::mutate(tissue_short = gsub(" .*", "", tissue)) %>% 
    dplyr::mutate(status = ifelse(tissue_short == "adenoma", "all", status)) %>% 
    dplyr::mutate(
      is_msi = 
        (case %in% .msi_positiv & tissue_short == "cancer") |
        (case %in% .msi_positiv_adenoma & tissue_short == "adenoma")
    ) %>% 
    dplyr::mutate(msi_status = ifelse(is_msi, "MSI", "MSS")) %>% 
    dplyr::mutate(filter_status = ifelse(is_filt, "filtered", "unfiltered"))
}

input_files = 
  file.path(dat_dir, "lifted_vcfs") %>% 
  list.files(full.names=TRUE, recursive = TRUE) %>%
  data.frame(file=.) %>% 
  cbind(., get_annot(.$file))


# split by all properties
group_id = with(input_files,paste(filter_status,tissue_short,status,msi_status,sep="."))
file_groups = split(input_files, group_id)

# split by all properties, but MSI status
group_id = with(input_files,paste(filter_status,tissue_short,status,"all",sep="."))
file_groups = c(file_groups, split(input_files, group_id))

# split by all properties, but MSI status
group_id = with(input_files,paste(filter_status,tissue_short,"all","all",sep="."))
file_groups = c(file_groups, split(input_files, group_id))


# split by all properties, but MSI status in neutral and selected cases
get_case_inference_status = function(x) {
  case_when(
    x %in% inference_neutral ~ "neutral",
    x %in% inference_selection ~ "selection",
    TRUE ~ as.character(NA)
  )
}

input_files_mod  = input_files %>% 
  dplyr::filter(status == "subclonal" & tissue_short == "cancer") %>% 
  dplyr::mutate(inf_status = get_case_inference_status(case)) %>%
  dplyr::filter(!is.na(inf_status)) %>%
  dplyr::mutate(tissue_short = paste0(tissue_short, "-", inf_status)) %>% 
  dplyr::mutate(group_id = paste(filter_status,tissue_short,status,"all",sep="."))

group_id = with(input_files_mod, paste(filter_status,tissue_short,status,"all",sep="."))
file_groups = c(file_groups, split(input_files_mod, group_id))

group_id = with(input_files_mod, paste(filter_status,tissue_short,status,msi_status,sep="."))
file_groups = c(file_groups, split(input_files_mod, group_id))


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# do dnds_analysis for each set and list of driver genes ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

remove_env_references_dnds_fit = function(x) { 
  attr(x$nbreg$terms, ".Environment") <- NULL
  attr(x$nbregind$terms, ".Environment") <- NULL
  attr(x$poissmodel$terms, ".Environment") <- NULL
  return(x)
}

for (i in seq_along(file_groups)) {
  try({
    
    file_group = names(file_groups)[i]
    print(file_group)
    
    out_file = file.path(dat_dir, "dndscv", paste0(file_group, ".rds"))
    if (file.exists(out_file)) next()
    
    d = file_groups[[i]]$file %>% 
      lapply(readRDS) %>% 
      do.call(what = rbind) %>%
      curate_data_for_dndscv()
    
    dnds_res = 
      dndscv_wrapper(d, outmats=T, max_muts=Inf, max_coding=Inf) %>% 
      remove_env_references_dnds_fit()
    
    dir.create(dirname(out_file), FALSE, TRUE)
    saveRDS(dnds_res, out_file)
  })
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Calc geneset dnds ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = file.path(dat_dir, "dndscv") %>% 
  list.files(full.names = TRUE, recursive = TRUE)


for (ifile in ifiles) {
  
  print(ifile)
  out_filenames = paste0(names(gene_lists), ".", basename(ifile))
  out_files = file.path(dat_dir, "genesetdnds", out_filenames)
  if (all(file.exists(out_files))) next()
  
  dnds_res = readRDS(ifile)
  
  for (i in seq_along(gene_lists)) {
    try({
      if (file.exists(out_files[i])) next()
      geneset_dnds = genesetdnds(dnds_res, gene_lists[[i]])
      dir.create(dirname(out_files[i]), FALSE, TRUE)
      saveRDS(geneset_dnds$globaldnds_geneset, out_files[i])
    })
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Variant number ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = file.path(dat_dir, "dndscv") %>% list.files(full=TRUE, recursive=TRUE)
ifiles = head(ifiles)
data_n_variants = load_n_variant_data_dnds(ifiles, gene_lists)
out_file = file.path(dat_dir, "number_of_variants_dnds.rds")
saveRDS(data_n_variants, out_file)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Save dnds annotations ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = 
  file.path(dat_dir, "dndscv") %>% 
  list.files("unfiltered[.](adenoma|cancer)[.].*[.]MS[SI].rds", full=TRUE, recursive=TRUE)

for (ifile in ifiles) {
  out_file = file.path(dat_dir, "dndscv_annot", basename(ifile))
  if (file.exists(out_file)) next()
  dir.create(dirname(out_file), FALSE, TRUE)
  
  annot_parts = strsplit(basename(ifile), "[.]")[[1]]
  oring_mut_dir = file.path(dat_dir, "lifted_vcfs")
  oring_mut_files = list.files(oring_mut_dir, annot_parts[[2]], full=TRUE, rec=TRUE) 
  oring_mut_files = grep(paste0("_", annot_parts[[1]], "_"), oring_mut_files, value = TRUE)

  
  annot_dndscv = readRDS(ifile)$annotmuts
  
  orig_muts_data = 
    lapply(oring_mut_files, readRDS) %>% 
    do.call(what=rbind) %>% 
    dplyr::mutate(mut=alt) %>% 
    dplyr::mutate(chr=gsub("chr", "", chr))
  

  merge_by = c("sampleID","chr","pos","ref","mut")
  d_merged = merge(annot_dndscv, orig_muts_data, by=merge_by, all.x=TRUE)
  
  saveRDS(d_merged, out_file)
}

