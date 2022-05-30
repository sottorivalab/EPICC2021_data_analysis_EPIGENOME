source("setup_environment/0-source.R")


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

vcf_subdir = "analysis/WGS/02-SNV/datasets/vcf_datasets"
fig_dir = "analysis/WGS/02-SNV/plots/msi_status"
dat_dir = "analysis/WGS/02-SNV/datasets/msi_status"
.msi_positiv = c("C536","C548","C516","C518","C552","C562")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Functions ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

parse_msisensor_file =
  function(f) {
    tryCatch(
      cbind(file = f, read.delim(f, stringsAsFactors = FALSE)) %>% 
        magrittr::set_colnames(c("file","n","n_somatic","frac_somatic")) %>% 
        dplyr::mutate(file=as.character(file)) %>%
        dplyr::mutate(sample=gsub("^C[0-9]*_", "", gsub("[_-]GRCh38$", "", basename(file)))) %>%
        cbind(THmisc::annotation_from_barcode(.$sample)),
      error = function(e) {
        print(paste0("Failed reading file ", f))
        NULL
      }
    )
  }

get_mutation_burden = function(f, clonal=FALSE, min_freq=0.05) {
  if (is.character(f)) d = readRDS(f) else d = f
  
  d = d[,!colnames(d) %in% .excluded_samples]
  annot = THmisc::annotation_from_barcode(colnames(d))
  wh_samples = annot$tissue_type == "cancer"
  d_mut = geno(d)$VAF[,wh_samples] > min_freq
  
  get_freq_types = function(x) {
    ids = c("SNV","MNV","InDel")
    table(THmisc::get_mutation_type(x))[ids] %>% 
      matrix(nrow=1) %>% 
      data.frame() %>% 
      magrittr::set_colnames(ids)
  }
  
  list(subclonal=any, clonal=all) %>% 
    lapply(function(x) get_freq_types(names(which(apply(d_mut, 1, x))))) %>% 
    reshape2::melt(measure.vars=c()) %>% 
    dplyr::mutate(status=L1, .keep="unused")
}

plot_structure = # 
  list(
    geom_jitter(width=0.25, height = ),
    xlab(""),
    scale_color_brewer(palette="Set1", direction = -1),
    background_grid(size.major=0.5),
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)),
    facet_grid(~msi_status, scales="free", space = "free")
  )

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# prepare output dir:
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load MSI-sensor data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Getting MSIsensor results ...\n")

# load all msisensor data
msisensor_data = 
  .msisensor_files %>%
  pbapply::pblapply(parse_msisensor_file) %>%
  do.call(what=rbind) %>%
  dplyr::filter(!sample_barcode %in% .excluded_samples & tissue_type != "normal") %>% 
  dplyr::select(-file) %>% 
  dplyr::mutate(msi_status=ifelse(patient %in% .msi_positiv, "MSI", "MSS")) %>% 
  dplyr::mutate(msi_status=factor(msi_status, c("MSS","MSI"), ordered = TRUE))

saveRDS(msisensor_data, file.path(dat_dir, "msisensor_data.rds"))

msisensor_avg_mut_per_pat = # average number of mutated MS per patient
  msisensor_data %>% 
  dplyr::filter(tissue_type == "cancer") %>%
  (function(x) split(x$frac_somatic, x$patient)) %>%
  sapply(mean) %>%
  sort()

saveRDS(msisensor_avg_mut_per_pat, file.path(dat_dir, "msisensor_avg_mut_per_pat.rds"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load VCF burden data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cat("Getting variant counts...\n")

somatic_variants_per_patient =
  vcf_subdir %>% 
  list.files("_filtered_.*[.]rds$", full.names=TRUE, recursive=TRUE) %>% 
  magrittr::set_names(basename(.)) %>% 
  pbapply::pblapply(get_mutation_burden) %>% 
  reshape2::melt(measure.vars=c()) %>% 
  dplyr::mutate(file_id=L1) %>% 
  dplyr::mutate(case=gsub("_.*", "", file_id)) %>% 
  dplyr::mutate(filter_status=grepl("_filt", file_id)) %>%
  dplyr::mutate(indel_snv_ratio = InDel/SNV) %>%
  dplyr::mutate(msi_status=ifelse(case %in% .msi_positiv, "MSI","MSS")) %>% 
  dplyr::mutate(msi_status=factor(msi_status, c("MSS","MSI")))

saveRDS(somatic_variants_per_patient, file.path(dat_dir, "somatic_variants_per_patient.rds"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Create plots ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

plt_mis_sensor_data = # plot in numeric order
  msisensor_data %>% 
  dplyr::mutate(patient=factor(patient, names(msisensor_avg_mut_per_pat))) %>% 
  dplyr::mutate(tissue_type = factor(tissue_type, levels(tissue_type), Hmisc::capitalize(levels(tissue_type)))) %>% 
  ggplot(aes(x=patient, y=frac_somatic/100, color=tissue_type, shape=tissue_type)) +
  plot_structure +
  geom_hline(yintercept=0.25, linetype=2, alpha=0.5) + 
  scale_color_grey(start = 0.6, end=0.2) + 
  ylab("Mutated MS") + 
  labs(shape="Tissue", color="Tissue")

out_file = file.path(fig_dir, "msisensor_frac_mutated.pdf")
ggsave(out_file, plt_mis_sensor_data, width=6.5, height=2.3)

out_file = file.path(fig_dir, "msisensor_frac_mutated.rds")
saveRDS(plt_mis_sensor_data, out_file)


plt_storted_n_vars = # plot in numeric order
  somatic_variants_per_patient %>% 
  dplyr::filter(filter_status == TRUE) %>% 
  dplyr::filter(status == "subclonal") %>% 
  dplyr::mutate(case=factor(case, names(msisensor_avg_mut_per_pat))) %>%
  reshape2::melt(measure.vars=c("SNV","MNV","InDel")) %>%
  ggplot(aes(x=case, y=value, color=variable, shape=variable)) + 
  plot_structure + 
  geom_hline(yintercept=0.25, linetype=2, alpha=0.5) + 
  scale_y_log10(limits=c(NA,NA)) + 
  ylab("N variants") + 
  labs(shape="Variant type", color="Variant type") 

out_file = file.path(fig_dir, "n_vars.pdf")
ggsave(out_file, plt_storted_n_vars, width=6.5, height=2.3)

out_file = file.path(fig_dir, "n_vars.rds")
saveRDS(plt_storted_n_vars, out_file)

