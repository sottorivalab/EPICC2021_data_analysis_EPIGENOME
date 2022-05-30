source("setup_environment/0-source.R")

#---------------------------------------------
# Options
#---------------------------------------------

dat_dir = "analysis/WGS/03-TREES/datasets"
fig_dir = "analysis/WGS/03-TREES/plots/trees/"

#---------------------------------------------

ifiles = list.files(dat_dir, "trees[.]rds$", full.names=TRUE, recursive=TRUE)
names(ifiles) = gsub("[.]rds", "", basename(ifiles))

#---------------------------------------------

data = 
  lapply(ifiles, readRDS) %>% 
  (function(x) x[sapply(x, function(y) length(y$data)) != 0])


for (set in names(data)) {
  for (remove_clonal in c(FALSE, TRUE)) {
    
    tree_plots = plot_treeset(data[[set]], remove_clonal_muts=remove_clonal)
    grop_list = do.call(grid::gList, lapply(tree_plots, ggplotGrob))
    
    # plot as panel of all trees
    panel = gridExtra::arrangeGrob(grobs=grop_list, nrow=6, ncol=5)
    suffix = ifelse(remove_clonal, "_no_clonal", "")
    out_file = file.path(fig_dir, paste0(set, suffix, ".pdf"))
    dir.create(dirname(out_file), FALSE, TRUE)
    ggsave(out_file, panel, width=8, height=11)
    
    # plot each case separately:
    for (pat_id in names(tree_plots)) {
      out_file = file.path(fig_dir, paste0(set, suffix), paste0(pat_id, ".pdf"))
      dir.create(dirname(out_file), showWarnings=FALSE, recursive=TRUE)
      ggsave(out_file, tree_plots[[pat_id]], width=4, height=2.8)
    }
  }
}

