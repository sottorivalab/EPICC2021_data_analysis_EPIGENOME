source("setup_environment/0-source.R")
library(VariantAnnotation)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

vcf_dir = "analysis/WGS/02-SNV/datasets/vcf_datasets"
dat_dir = "analysis/WGS/02-SNV/datasets/dnds"
fig_dir = "analysis/WGS/02-SNV/plots/dnds"
ccf_cutoff = 0.2
n_cores = 8

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# prepare output dir  ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load dNdS results ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

geneset_results = file.path(dat_dir, "genesetdnds") %>% 
  list.files(full.names = TRUE, recursive = TRUE) %>% 
  lapply(load_genset_result_file) %>% 
  do.call(what=rbind) %>% 
  dplyr::mutate(msi_status = factor(msi_status, c("MSS","MSI","all"))) %>%
  dplyr::filter(cilow>0 | is.finite(cihigh)) %>% 
  dplyr::mutate(group = paste0(tissue, ".", clonal_status)) %>% 
  dplyr::mutate(variant_group=factor(group, names(dnds_group_labels), dnds_group_labels)) %>% 
  dplyr::mutate(gene_set_label=factor(factor(geneset, names(gene_set_labels), gene_set_labels)))


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plot dNdS results ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

for (fstat in unique(geneset_results$filter_stat)) {
  
  tile = paste0("dndncv - EPICC somatic variants (", fstat, ")")
  d_cur = geneset_results %>% dplyr::filter(filter_stat == fstat)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # plot of all gene sets ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  plot_all_dnds_data = plot_dnds_results(d_cur, tile)
  out_file = file.path(fig_dir, paste0("tumour_somatic_dndscv_", fstat, ".pdf"))
  ggsave(out_file, plot_all_dnds_data, width=18, height=7)
  
  out_file = file.path(fig_dir, paste0("tumour_somatic_dndscv_", fstat, ".rds"))
  saveRDS(plot_all_dnds_data, out_file)
  
    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # subset of data for chromatin paper ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  genesets = c("All", "Chromatin modifiers (all)")
  labels = c("All", "Chromatin modifiers")
  cols = c("#e63946","#457b9d")
  sample_group = c("Adenoma", "Cancer\nclonal", "Cancer\nsubclonal")
  
  d_cur_m = d_cur %>% 
    dplyr::filter(group != "cancer.all") %>%
    dplyr::filter(msi_status != "all") %>% 
    dplyr::filter(gene_set_label %in% genesets) %>% 
    dplyr::filter(name %in% c("wmis","wtru")) %>%
    dplyr::filter(variant_group %in% sample_group)
    
  plt = 
    plot_dnds_results(d_cur_m, tile) + 
    scale_y_log10(limits=c(NA, NA)) + 
    scale_color_manual(values = cols,  breaks = genesets, labels=labels) + 
    theme(axis.text.x = element_text(hjust=0.5))
  
  out_file = file.path(fig_dir, paste0("cmg_dnds_", fstat, ".pdf"))
  ggsave(out_file, plt, width=7, height=4.0)
  
  out_file = file.path(fig_dir, paste0("cmg_dnds_", fstat, "_small.pdf"))
  ggsave(out_file, plt, width=5.76, height=3.15)
  
  out_file = file.path(fig_dir, paste0("cmg_dnds_", fstat, ".rds"))
  saveRDS(plt, out_file)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # subset of data for inference paper ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  genesets = c("All","IntOGen","Martincorena et al. (2017)")
  cols = c("#e63946","#457b9d","#83c5be")
  sample_group = c("Adenoma", "Cancer\nclonal", "Cancer\nsubclonal")
  
  d_cur_m = d_cur %>% 
    dplyr::filter(group != "cancer.all") %>%
    dplyr::filter(msi_status != "all") %>% 
    dplyr::filter(gene_set_label %in% genesets) %>% 
    dplyr::filter(name %in% c("wmis","wtru")) %>%
    dplyr::filter(variant_group %in% sample_group)
  
  plt = 
    plot_dnds_results(d_cur_m, tile) + 
    scale_y_log10(limits=c(NA, NA)) + 
    scale_color_manual(values = cols,  breaks = genesets)
  
  out_file = file.path(fig_dir, paste0("inference_paper_dnds_", fstat, ".pdf"))
  ggsave(out_file, plt, width=8, height=4.5)
  
  out_file = file.path(fig_dir, paste0("inference_paper_dnds_", fstat, "_small.pdf"))
  ggsave(out_file, plt, width=5.76, height=3.5)
  
  out_file = file.path(fig_dir, paste0("inference_paper_dnds_", fstat, ".rds"))
  saveRDS(plt, out_file)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # subset of all driver gens ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  genesets = c(
    "All",
    "TCGA (2012)",
    "IntOGen",
    "Cross et al. (2018)",
    "Martincorena et al. (2017)",
    "Tarabichi et al. (2018)"
  )
  
  cols = c("#e63946","#1d3557","#457b9d","#83c5be","#fca311")
  
  d_cur_m = d_cur %>% 
    dplyr::filter(group != "cancer.all") %>%
    dplyr::filter(msi_status != "all") %>% 
    dplyr::filter(gene_set_label %in% genesets) %>% 
    dplyr::filter(name %in% c("wmis","wtru"))
  
  plt = 
    plot_dnds_results(d_cur_m, tile) + 
    scale_y_log10(limits=c(NA, NA))
  
  out_file = file.path(fig_dir, paste0("driver_gene_dnds_", fstat, ".pdf"))
  ggsave(out_file, plt, width=10, height=4.5)
  
  out_file = file.path(fig_dir, paste0("driver_gene_dnds_", fstat, ".rds"))
  saveRDS(plt, out_file)
  
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plot dNdS results (split on inference) ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

for (fstat in unique(geneset_results$filter_stat)) {
  
  new_group_labels = c(
    "cancer.subclonal" = "Cancer subclonal\n(all)",
    "cancer-selection.subclonal" = "Cancer subclonal\n(selected)",
    "cancer-neutral.subclonal" = "Cancer subclonal\n(neutral)"
  )
  
  selected_genesets = c(
    "All", 
    "IntOGen",
    "Cross et al. (2018)",
    "Martincorena et al. (2017)",
    "Tarabichi et al. (2018)"
  )
  
  tile = paste0("dndncv - ABC-SMC Classification (", fstat, ")")
  
  d_cur = geneset_results %>% 
    dplyr::filter(filter_stat == fstat) %>% 
    dplyr::mutate(variant_group=factor(group, names(new_group_labels), new_group_labels)) %>% 
    dplyr::filter(gene_set_label %in% selected_genesets)
    
  plt = plot_dnds_results(d_cur, tile)
  out_file = file.path(fig_dir, paste0("abcsmc_dnds_full_", fstat, ".pdf"))
  ggsave(out_file, plt, width=12, height=8)
  
  out_file = file.path(fig_dir, paste0("abcsmc_dnds_full_", fstat, ".rds"))
  saveRDS(plt, out_file)
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Selected subset of dnds data ####
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  genesets = c("All","IntOGen")
  cols = c("#e63946","#457b9d")

  d_cur = geneset_results %>% 
    dplyr::filter(filter_stat == fstat) %>% 
    dplyr::mutate(variant_group=factor(group, names(new_group_labels), new_group_labels)) %>% 
    dplyr::filter(gene_set_label %in% genesets) %>% 
    dplyr::filter(msi_status != "all") %>% 
    dplyr::filter(name %in% c("wmis","wtru"))

  plt = plot_dnds_results(d_cur, tile) +  scale_y_log10(limits=c(0.03, 100)) + theme(axis.text.x = element_text(hjust=0.5))
  out_file = file.path(fig_dir, paste0("abcsmc_dnds_paper_", fstat, ".pdf"))
  ggsave(out_file, plt, width=7, height=5)
  
  out_file = file.path(fig_dir, paste0("abcsmc_dnds_full_", fstat, ".rds"))
  saveRDS(plt, out_file)
  
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Variant heatmap of CMG genes ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

data = rbind(
  readRDS("analysis/WGS/02-SNV/datasets/dnds/dndscv_annot/unfiltered.adenoma.all.MSS.rds"),
  readRDS("analysis/WGS/02-SNV/datasets/dnds/dndscv_annot/unfiltered.cancer.clonal.MSS.rds"),
  readRDS("analysis/WGS/02-SNV/datasets/dnds/dndscv_annot/unfiltered.cancer.subclonal.MSS.rds")
) %>% curate_data_for_dndscv()

cmg_genes = 
  "external_datasets/chromatin_modifer_genes.txt" %>% 
  readr::read_csv(col_names = FALSE) %>% 
  unlist() 

plot_heatmap = 
  plot_variant_heatmap_dndscv(
    data, 
    cmg_genes[cmg_genes %in% data$gene[data$impact %in% c("Missense","Nonsense")]]
  )

out_file = file.path(fig_dir, "chromatin_modifiers_fill_by_fraction_mutated.pdf")
ggsave(out_file, plot_heatmap, width=6, height=9.865)

out_file = file.path(fig_dir, "chromatin_modifiers_fill_by_fraction_mutated.rds")
saveRDS(plot_heatmap, out_file)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plot number of variants per type ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

data_n_variants = readRDS(file.path(dat_dir, "number_of_variants_dnds.rds"))
plot_n_variants = plot_n_variants_per_type_group(data_n_variants)
out_file = file.path(fig_dir, "number_of_variants.pdf")
ggsave(out_file, plot_n_variants, width=9.2, height=9)
