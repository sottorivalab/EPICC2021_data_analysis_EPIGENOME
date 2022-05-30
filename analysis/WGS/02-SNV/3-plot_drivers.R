source("setup_environment/0-source.R")
library(cowplot)
library(ggplot2)
library(VariantAnnotation)
library(dplyr)
library(reshape2)
library(seqinr)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

driver_data_file = "analysis/WGS/02-SNV/datasets/driver_gene_data.rds"
fig_dir = "analysis/WGS/02-SNV/plots"
out_dir = "analysis/WGS/02-SNV/datasets/driver_tables"
ccf_cutoff = 0.2

gene_lists =
  c(
    "IntOGen-DriverGenes_COREAD",
    "coad_driver_genes_tcga_2012",
    "CRC_drivers_Cross2018_Tier1",
    "driver_genes_martincorena_2017",
    "driver_genes_tarabichi_2018"
  )

mmr_genes = c("PMS2", "MLH1", "MSH6", "MSH2", "MSH3", "MLH3", "POLE")
genes_from_lists = THmisc::gene_lists[gene_lists]

genes_from_lists$`IntOGen-DriverGenes_COREAD` = 
  genes_from_lists$`IntOGen-DriverGenes_COREAD` %>% 
  (function(x) x[!x %in% c("LRP1B","KMT2C","PARP4")])

genes_from_lists$`IntOGen-DriverGenes_COREAD_plusMMR` = 
  c(genes_from_lists$`IntOGen-DriverGenes_COREAD`, mmr_genes)

genes_from_lists$`CRC_drivers_Cross2018_Tier1_plusMMR` = 
  c(genes_from_lists$`CRC_drivers_Cross2018_Tier1`, mmr_genes)

genes_from_lists[["chromatin_modifers_selected"]] = 
  grep(
    c("^(ARID.*)|(KMT.*)|(KDM.*)|(TET.*)|(DOT.*)|(KAT.*)|(DNMT.*)"), 
    THmisc::gene_lists$all_genes, 
    value = TRUE
  )

genes_from_lists$chromatin_modifer_genes = 
  "external_datasets/chromatin_modifer_genes.txt" %>% 
  readr::read_csv(col_names = FALSE) %>% 
  unlist() %>% magrittr::set_names(NULL)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

update_annots = function(x) {
  annot = THmisc::annotation_from_barcode(x$sample)
  cbind(x[,!colnames(x) %in% colnames(annot),drop=FALSE], annot)
}

vcf_driver_data_frames = readRDS(driver_data_file) %>% lapply(update_annots)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Create per case heatmaps  ####

for (i in seq_along(vcf_driver_data_frames)) {
  sid = names(vcf_driver_data_frames)[i]
  
  for (j in seq_along(genes_from_lists)) {
    
    data = 
      vcf_driver_data_frames[[i]] %>% 
      dplyr::mutate(group=region) %>% 
      dplyr::filter(!(tissue_type == "normal" & sample_type == "G")) %>% 
      dplyr::filter(mutation %in% mutation[NV > 0])
    
    out_file = 
      file.path(
        fig_dir, 
        "driver_vaf_heatmaps",
        names(genes_from_lists)[j], 
        ifelse(grepl("_unfiltered", sid), "Unfiltered", "Filtered"),
        paste0(sid, ".pdf")
      )

    plot = 
      THmisc::save_driver_vaf_heatmap(
        d = data, 
        out_file = out_file, 
        genes = genes_from_lists[[j]]
      )
    
    saveRDS(plot, gsub("[.]pdf$", ".rds", out_file))
    
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Mutation heatmaps across all cases ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

wh_uf = grepl("_unfiltered", names(vcf_driver_data_frames))

vcf_driver_data_frames_mod = 
  vcf_driver_data_frames[wh_uf] %>%
  lapply(function(d){
    wh = is.na(d$CCF)
    d$CCF[wh] = d$VAF[wh] * 2 / get_purity(d$sample_barcode[wh])
    d[d$tissue_type != "normal" & !d$sample %in% .excluded_samples,]
  })

make_patient_gr = function(msi, pat) {
  level_order = c("MSS\nAdenoma","MSS\nCancer","MSI\nAdenoma","MSI\nCancer")
  levels = paste0(msi, "\n", Hmisc::capitalize(as.character(pat)))
  factor(levels, level_order, ordered = TRUE)
}

make_final_type = function(type, variant, split_indel_trunc=FALSE) {
  type = as.character(type)
  type[type == "MNV"] = "SNV"
  if (split_indel_trunc) {
    type[grepl("[*]", variant)] = "Truncating"
    type = factor(type, c("InDel","SNV","Truncating"), ordered = TRUE)
  } else {
    type[type == "InDel" | grepl("[*]", variant)] = "InDel/Truncating"
    type = factor(type, c("InDel/Truncating","SNV"), ordered = TRUE)
  }
  return(type)
}

gene_to_gr = function(x) { # labeller function for gene groups
  msi_gene = x %in% mmr_genes
  ifelse(msi_gene, "MMR", "Drivers")
}


d_for_heatmap = # modify patient groups
  do.call(rbind, vcf_driver_data_frames_mod) %>% 
  dplyr::mutate(patient_gr = make_patient_gr(msi_status, tissue_type)) %>%
  dplyr::mutate(type = make_final_type(type, variant)) %>% 
  dplyr::mutate(patient = gsub("[.][a-z]+", " ", tumour_id)) # replaces tissue id from barcode

# modify curated mutations:
ids = with(d_for_heatmap, paste0(tumour_id, "-", mutation))
d_for_heatmap[ids %in% .variants_clonal_hg38, c("CCF","VAF")] = 1
d_for_heatmap[ids %in% .variants_drop_hg38, c("CCF","VAF")] = 0

for (j in seq_along(genes_from_lists)) {
  
  if (NROW(d_for_heatmap) == 0)
    next()
  
  driver_plot = 
    THmisc:::plot_gene_mut_heatmap(
      d_for_heatmap,
      ccf_cutoff,
      genes = genes_from_lists[[j]],
      gene_to_gr = gene_to_gr, 
      xlab = "Patient"
    )
  
  out_file = 
    file.path(
      fig_dir, 
      "driver_panel",
      paste0( names(genes_from_lists)[j], ".pdf")
    )
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  
  n_pats = length(unique(paste0(d_for_heatmap$patient, d_for_heatmap$tissue_type)))
  n_genes = length(unique(genes_from_lists[[j]]))
  n_types = length(unique(d_for_heatmap$patient_gr))
  n_gensets = length(unique(gene_to_gr(genes_from_lists[[j]])))
  
  width = n_pats * 0.15 + 1.612679 + 0.5
  height =  n_genes * 0.122 + 4 + (n_gensets-1) * 0.4
  
  if (basename(out_file) == "IntOGen-DriverGenes_COREAD_plusMMR.pdf") {
    height = 14.5
  } else if (basename(out_file) == "CRC_drivers_Cross2018_Tier1.pdf") {
    height = 5.45
  } else if (basename(out_file) == "CRC_drivers_Cross2018_Tier1_plusMMR.pdf") {
    height = 6.72
    stop()
  }
  
  ggplot2::ggsave(
    out_file, 
    driver_plot, 
    width = width,
    height =  height, 
    limitsize = FALSE
  )
  
  saveRDS(driver_plot, gsub("[.]pdf$", ".rds", out_file))
}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Alternative mutation heatmaps across all cases ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

d_for_heatmap_alt = # modify patient groups
  d_for_heatmap %>% 
  dplyr::mutate(type = THmisc::get_mutation_type(mutation)) %>% 
  dplyr::mutate(type = make_final_type(type, variant, TRUE)) %>% 
  dplyr::filter(msi_status == "MSS") %>% 
  dplyr::mutate(type = factor(type, rev(c("SNV","Truncating")))) %>% 
  dplyr::filter(type != "InDel")

driver_plot = 
  THmisc:::plot_gene_mut_heatmap(
    d_for_heatmap_alt,
    ccf_cutoff,
    genes = genes,
    gene_to_gr = gene_to_gr
  )

out_file = file.path(fig_dir, "driver_panel", "final_all_cmg.pdf")
ggplot2::ggsave(out_file, driver_plot, width=6, height=11.05)
saveRDS(driver_plot, gsub("[.]pdf$", ".rds", out_file))


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Dump datasets of the drivers ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

driver_table = 
  driver_data = 
  d_for_heatmap %>% 
  dplyr::mutate(gene=sapply(strsplit(gene,";"), unique)) %>% 
  dplyr::filter(gene %in% genes_from_lists$`IntOGen-DriverGenes_COREAD`) %>%
  group_by(tumour_id, msi_status, mutation, variant, gene) %>% 
  summarise(mutated = paste0(sum(CCF > ccf_cutoff), "/", length(CCF)),
            n_mut = sum(CCF > ccf_cutoff),
            clonal = all(CCF > ccf_cutoff)) %>%
  dplyr::filter(n_mut > 0) %>%
  data.frame() %>% 
  dplyr::select(
    Patient = tumour_id,
    MSI_status = msi_status,
    Mutation = mutation,
    Consequence = variant,
    Gene = gene, 
    Mutated_samples = mutated,
    Clonal = clonal
  )

WriteXLS::WriteXLS(driver_table, file.path(out_dir, "all_drivers.xlsx"))
readr::write_tsv(driver_table, file.path(out_dir, "all_drivers.tsv"))

driver_sc = driver_table %>% dplyr::filter(!Clonal & grepl("cancer", Patient))
WriteXLS::WriteXLS(driver_sc, file.path(out_dir, "subclonal_driver.xlsx"))
readr::write_tsv(driver_sc, file.path(out_dir, "subclonal_driver.tsv"))

