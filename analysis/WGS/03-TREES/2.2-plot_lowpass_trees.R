source("setup_environment/0-source.R")

#---------------------------------------------
# Options
#---------------------------------------------

idat_dir = file.path("analysis/WGS/03-TREES/datasets/", c("lp_trees","lp_and_atac_trees"))
fig_dir = "analysis/WGS/03-TREES/plots/"

#---------------------------------------------

ifiles = list.files(idat_dir, "[.]rds$", full.names=1, recursive=1)
ifiles = ifiles[grepl("ml_sample_tree", ifiles)]
names(ifiles) = gsub("[.]rds", "", basename(ifiles))

#---------------------------------------------

plot_cn = function(d) {
  table(cn=d$cn, mm=d$mm) %>% 
    reshape2::melt() %>% 
    dplyr::mutate(lab = paste0(value, " (", round(value/sum(value)*100, 2), "%)")) %>% 
    ggplot(aes(x=cn, y=mm, fill=value, label=lab)) + 
    scale_fill_distiller(palette = 4) + 
    cowplot::theme_cowplot() + 
    geom_tile() + 
    geom_text() + 
    xlab("CN") + 
    ylab("Multiplicity") + 
    labs(fill="N")
}

labeller_function = 
  function(x) gsub("EPICC_", "", x)

for (i in seq_along(ifiles)) {
  
  # construct output file:
  analyte_set = basename(dirname(dirname(dirname(ifiles[i]))))
  file_set = basename(dirname(dirname(ifiles[i])))
  out_name = gsub("[.]rds", ".pdf", basename(ifiles[i]))
  out_file = file.path(fig_dir, file_set, out_name)
  pat = gsub("_.*", "", basename(ifiles[i]))
  
  # load data
  res_data = readRDS(ifiles[i])
  cn_file = file.path(dirname(dirname(ifiles[i])), "cn_states", basename(ifiles[i]))
  cn_states = readRDS(cn_file)
  
  # plot cn state
  cn_plot = plot_cn(cn_states)
  out_file = file.path(fig_dir, analyte_set, file_set, "cn_states", out_name)
  dir.create(dirname(out_file), FALSE, TRUE)
  ggsave(out_file, cn_plot, height = 4.5, width = 6)
  
  # plot trees
  ggplot_tree = res_data$tree %>% MLLPT:::remove_root_tip("GL") %>% plot_tree()
  out_file = file.path(fig_dir, analyte_set, file_set, "lp_tree", out_name)
  dir.create(dirname(out_file), FALSE, TRUE)
  ggsave(out_file, ggplot_tree, height = 3, width = 4)
  
  ggplot_tree = res_data$inital_values$tree %>% MLLPT:::remove_root_tip("GL") %>% plot_tree()
  out_file = file.path(fig_dir, analyte_set, file_set, "original_tree", out_name)
  dir.create(dirname(out_file), FALSE, TRUE)
  ggsave(out_file, ggplot_tree, height = 3, width = 4)
  
  # plot log-lik data
  plot_loglik = 
    MLLPT::plot_lp_loglik(res_data, labeller_function) + 
    ggtitle(paste0(pat, " - LP assignment")) + 
    xlab("")
  
  out_file = file.path(fig_dir, analyte_set, file_set, "log_lik_data", out_name)
  dir.create(dirname(out_file), FALSE, TRUE)
  ggsave(out_file, plot_loglik, height = 6, width = 6)
  
  
  # plot log-lik data
  plot_loglik_edge = 
    MLLPT::plot_lp_loglik_edge(res_data, labeller_function) + 
    ggtitle(paste0(pat, " - LP assignment")) + 
    theme(strip.text.y = element_text(angle=0)) + 
    scale_y_continuous(n.breaks = 3)
  
  old_level = levels(plot_loglik_edge$layers[[2]]$data$edge)
  new_level =  gsub(" [(].*", "", old_level)
  
  plot_loglik_edge$data$edge = 
    factor(plot_loglik_edge$data$edge, old_level, new_level)
  
  plot_loglik_edge$layers[[2]]$data$edge = 
    factor(plot_loglik_edge$layers[[2]]$data$edge , old_level, new_level)
  
  out_file = file.path(fig_dir, analyte_set, file_set, "log_lik_data_edge", out_name)
  dir.create(dirname(out_file), FALSE, TRUE)
  ggsave(out_file, plot_loglik_edge, height=12, width=10)
  
  
  # plot params
  optimisation_plot =
    MLLPT::plot_sample_data(
      res_data,
      labeller_function = labeller_function,
      external_purity_estimate = get_purity(res_data$per_sample_results$sample)
    )
  out_file = file.path(fig_dir, analyte_set, file_set, "parameter_optimisation", out_name)
  dir.create(dirname(out_file), FALSE, TRUE)
  ggsave(out_file, optimisation_plot, height = 4, width = 9)
  
  # tree plus edge labels
  out_file = file.path(fig_dir, analyte_set, file_set, "tree_plus_edge_label", out_name)
  dir.create(dirname(out_file), FALSE, TRUE)
  
  pdf(out_file)
    plot(res_data$inital_values$tree)
    ape::edgelabels()
  dev.off()
  
  
  # tree plus cn data
  cn_plot = MLLPT::plot_tree_and_cn_data(res_data, .cna_data, pat)
  out_file = file.path(fig_dir, analyte_set, file_set, "cna_plus_tree", out_name)
  dir.create(dirname(out_file))
  ggsave(out_file, cn_plot, width = 14, height = 6)
  
}
