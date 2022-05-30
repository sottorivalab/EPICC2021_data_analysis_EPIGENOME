source("setup_environment/0-source.R")
library(VariantAnnotation)
library(THmisc)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

fig_dir = "analysis/WGS/02-SNV/plots" # output dir for figures
vcf_dir = "analysis/WGS/02-SNV/datasets/vcf_datasets"
mut_mapping_file = "analysis/WGS/02-SNV/datasets/all_drivers_list.rds"
geno = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
data_mut_map = readRDS(mut_mapping_file)
mutation_map = with(data_mut_map, magrittr::set_names(gsub(";.*", "", variant), id))

construct_out_path = function(x, plot) {
  filt_status = ifelse(grepl("_filtered", x), "filtered", "unfiltered")
  file_name = gsub("[.]rds$", ".pdf", basename(x))
  out_file = file.path(fig_dir, plot, filt_status, file_name)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Heatmaps of data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = list.files(vcf_dir, "[.]rds$", full.names = TRUE)
for (vcf_file in ifiles) {
  
  out_file = gsub("[.]pdf", ".png", construct_out_path(vcf_file,  "heatmap"))
  if (file.exists(out_file)) next()
  title = paste0("Somatic mutations - ", gsub("_.*", "", basename(out_file)))
  
  # load data
  data = readRDS(vcf_file)
  
  # if filtered file add missing drivers from unfiltered file
  if (grepl("_filtered", basename(vcf_file))) {
    uf_file = gsub("_filtered", "_unfiltered", vcf_file)
    d_uf = readRDS(uf_file)
    wh_add = rownames(d_uf) %in% names(mutation_map) & !rownames(d_uf) %in% names(data)
    data = rbind(data, d_uf[wh_add,])
  }
  
  pl = # plot cna state data
    THmisc::vcf_to_data_frame(data) %>% 
    cbind(., THmisc::annotation_from_barcode_epicc(as.character(.$sample))) %>% 
    dplyr::filter(tissue_type != "normal") %>% 
    add_cluster_infos(small_frac=0.01, drop_all_sc=TRUE) %>%
    dplyr::mutate(group=region) %>% 
    plot_mutation_heatmap(annot = mutation_map, title=title, value="CCF") + 
    theme(axis.text.y=element_text(size=8))
  
  # save figure
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_file, pl, width=7, height=14)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plots of VAF distribution across CN states ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = list.files(vcf_dir, "[.]rds$", full.names = TRUE)
for (vcf_file in ifiles) {
  
  out_file = construct_out_path(vcf_file,  paste0("histograms_per_cn/", c("VAF","CCF")))
  if (all(file.exists(out_file))) next()
  
  data = THmisc::vcf_to_data_frame(vcf_file)
  
  # get sample information
  samples = unique(data$sample)
  purity = magrittr::set_names(get_purity(samples), samples)
  ploidy = magrittr::set_names(get_ploidy(samples), samples)
  
  for (what in c("CCF","VAF")) {
    
    out_file = construct_out_path(vcf_file,  paste0("histograms_per_cn/", what))
    if (file.exists(out_file)) next()
    
    plot_list = 
      plot_vaf_distribution(
        data,
        what = what,
        purity = purity,
        ploidy = ploidy, 
        max_cn = 8
      )
    
    # save all plots to one file:
    grobs = do.call(grid::gList, plot_list)
    p = gridExtra::marrangeGrob(grobs, nrow=1, ncol=1)
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    ggsave(out_file, p, width=9, height=6.5)
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plots of CN states across samples ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = list.files(vcf_dir, "[.]rds$", full.names = TRUE)
for (vcf_file in ifiles) {
  
  out_file = construct_out_path(vcf_file,  "cn_states")
  if (file.exists(out_file)) next()
  
  # plot cna state data
  data = THmisc::vcf_to_data_frame(vcf_file)
  pl = plot_cn_states_across_samples(data)
  
  # save figure
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_file, pl, width=20, height=6)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Histograms of data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = list.files(vcf_dir, "[.]rds$", full.names = TRUE)
for (vcf_file in ifiles) {
  
  out_file = construct_out_path(vcf_file,  paste0("histograms/", c("VAF","CCF")))
  if (all(file.exists(out_file))) next()
  
  filt_status = ifelse(grepl( "_filtered", vcf_file), "filtered", "unfiltered")
  title = paste0(gsub("_.*", "", basename(vcf_file)), " (", filt_status, ")")
  
  data_plot = 
    THmisc::vcf_to_data_frame(vcf_file) %>%
    THmisc::add_cluster_infos() %>% 
    cbind(., THmisc::annotation_from_barcode_epicc(as.character(.$sample))) %>% 
    dplyr::filter(tissue_type != "normal")
    
  for (value in c("VAF","CCF")) {
    out_file = construct_out_path(vcf_file,  paste0("histograms/", value))
    if (file.exists(out_file)) next()
    pl = THmisc::plot_mutation_histogram(data_plot, value, max_value = 3, title)
    n_rows = ceiling(length(unique(pl$data$sample))/4)
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    ggsave(out_file, pl, width=12, height=1+4.5/2*n_rows)
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Waterfall plots of mutations ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

ifiles = list.files(vcf_dir, "[.]rds$", full.names = TRUE)
for (vcf_file in ifiles) {
  out_file = construct_out_path(vcf_file, "mm_distances")
  if (all(file.exists(out_file))) next()
  filt_status = ifelse(grepl( "_filtered", vcf_file), "filtered", "unfiltered")
  title = paste0(gsub("_.*", "", basename(vcf_file)), " (", filt_status, ")")
  pl = readRDS(vcf_file) %>% THmisc::plot_inter_mutation_distance(title, geno)
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_file, pl, width=19, height=5)
}

