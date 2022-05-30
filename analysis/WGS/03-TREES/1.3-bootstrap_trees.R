source("setup_environment/0-source.R")
library(phangorn)

#---------------------------------------------
# Options
#---------------------------------------------

idat_dir = "analysis/WGS/03-TREES/datasets/mp_trees"
dat_dir = "analysis/WGS/03-TREES/datasets/bootstrap_support"
fig_dir = "analysis/WGS/03-TREES/plots/bootstrap_support"
n_cores = 6

#---------------------------------------------

dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

ifiles = list.files(idat_dir, "[.]rds$", full.names = TRUE, recursive = TRUE)
names(ifiles) = gsub("[.]rds$", "", basename(ifiles))
pyData = lapply(ifiles, readRDS)

#---------------------------------------------

set.seed(123)

.get_bs_tree = function(x) {
  bootstrap.phyDat(x, 
    phyDat_to_tree, 
    bs = 100, 
    mc.cores = n_cores
  )
}

bs_data = 
  lapply(pyData, {
    function(x) {
      pbapply::pblapply(x$data, .get_bs_tree)
    }
  })


#---------------------------------------------

for (i in seq_along(bs_data)) {
  plot_list = list()
  for (j in seq_along(bs_data[[i]])) {
    
    id = names(bs_data[[i]])[j]
    tree = pyData[[i]]$trees[[id]]
    tree_ = tree
    bs_dat = bs_data[[i]][[id]]
    
    try({
      tl = gsub("EPICC_C[0-9]+_", "", gsub("_[LD][0-9]$", "", tree$tip.label))
      bs_val = signif(prop.clades(tree, bs_dat) / NROW(bs_dat), 2)
      tree_$node.label = ifelse(is.na(bs_val), "", paste0(bs_val, ""))
      
      wh1 = getRoot(tree) - length(tree$tip.label)
      wh_gl_tip = which(tree$tip.label == "GL")
      wh2 = which(tree$edge[,1] == getRoot(tree) & tree$edge[,2] == wh_gl_tip)
      
      val2 = bs_val[wh2 - length(tree_$tip.label) + 1]
      bs_val[wh2 - length(tree_$tip.label) + 1] = bs_val[wh1]
      bs_val[wh1] = val2
      
      tree_$node.label = ifelse(is.na(bs_val), "", paste0(bs_val, ""))
      tree_$tip.label = tl
    })
    
    plot_list[[j]] = 
      MLLPT:::remove_root_tip(tree_, "GL") %>% 
      MLLPT::plot_tree() +
      ggtree::geom_nodelab(size=2, color="red", nudge_x=0.1) + 
      ggtitle(id)
  }
  
  # save all plots to one file:
  gr_list = do.call(grid::gList, lapply(plot_list, ggplotGrob))
  plots = gridExtra::arrangeGrob(grobs=gr_list, nrow=5, ncol=6)
  out_file = file.path(fig_dir, paste0(names(ifiles)[i], ".pdf"))
  ggsave(out_file, plots, width=10, height=7.5)
}

