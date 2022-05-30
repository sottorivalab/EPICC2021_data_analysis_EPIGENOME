source("setup_environment/0-source.R")

library(ggtree)
library(dplyr)
library(ggplot2)
library(phangorn)
library(tidytree)
library(reshape)
library(ggridges)
library(gridExtra)
library(grid)
library(VariantAnnotation)
library(pheatmap)
library(magrittr)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

tree_file_original = "analysis/WGS/03-TREES/datasets/mp_trees/mp_trees_filtered_including_bulks_trees.rds"
tree_file_dir = "analysis/WGS/03-TREES/datasets/lp_trees/mp_trees_filtered_including_bulks_trees/ml_sample_trees"
fig_dir = "analysis/WGS/03-TREES/plots/driver_annotated_trees"
driver_data_file = "analysis/WGS/02-SNV/datasets/driver_gene_data.rds" 
status = NULL

gene_sets = c("coad_driver_genes_tcga_2012","IntOGen-DriverGenes_COREAD")
gene_lists = THmisc::gene_lists[gene_sets]
gene_lists[["most_frequent"]] = c("TP53","APC","KRAS","PIK3CA")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# load purity estimates from genotpying:
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

driver_data = 
  readRDS(driver_data_file) %>% 
  (function(x) x[!sapply(x, is.null)]) %>%
  reshape2::melt(measure.vars=c()) %>%
  (function(d) lapply(gene_lists, function(x) d[d$gene %in% x,])) %>% 
  lapply(function(x) x$variant %>% set_names(x$mutation)) %>% 
  lapply(sort)

# load the file containing the original trees
tree_and_pydata = readRDS(tree_file_original)
tree_data = tree_and_pydata$trees %>% lapply(MLLPT:::remove_root_tip, "GL")
py_data = tree_and_pydata$data

# 
for (gene_set in names(driver_data)) {
  
  for (i in seq_along(tree_data)) {
    
    tree = tree_data[[i]]
    data_tree = py_data[[i]]

    anc = ancestral.pars(tree, data_tree)
    tree$edge.label = rep("", length(tree$edge.length))
    
    var_idx = attr(anc, "index")
    mt_drivers = driver_data[[gene_set]][attr(anc, "id")]
    var_idx = var_idx[order(mt_drivers)]
    mt_drivers = mt_drivers[order(mt_drivers)]
    
    for (j in which(!is.na(mt_drivers))) {
      var_idx_j = var_idx[j] 
      lab_j = mt_drivers[j]
      
      state_nodes = sapply(anc, function(x) if (NROW(x)) x[var_idx_j,2] else NA)
      mutated_nodes = names(which(state_nodes == 1))
      if (length(mutated_nodes) == 0) next()
      
      first_mutated_node = unique(unlist(lapply(mutated_nodes, function(n) {
        na = get_nodes_above(n, tree)
        tail(na[na %in% mutated_nodes], 1)
      })))
      
      stopifnot(length(first_mutated_node) >= 1)
      if(length(first_mutated_node) > 1) warning(paste0(lab_j, " annotated ", length(first_mutated_node), " times in case ", names(tree_data)[i], ".\n"))
      
      for (mn in first_mutated_node) {
        edge_j = c(which(tree$edge[,"node"] == mn), which(tree$tip.label == mn))
        stopifnot(length(edge_j) == 1)
        tree$edge.label[edge_j] = paste0(tree$edge.label[edge_j], "\n", lab_j)
        
        if (length(first_mutated_node) > 1)
          tree$edge.label[edge_j] = paste0(tree$edge.label[edge_j], "!")
      }       
      
    }
    tree$edge.label = gsub("^\n", "", tree$edge.label)
    tree_data[[i]] = tree
  }
  
  
  trees = tree_data %>% lapply(plot_tree, color_by=status) 
  legend = cowplot::get_legend(trees[[1]])  # keep color legend and remove them from all trees
  for (i in seq_along(trees)) { trees[[i]] = trees[[i]] + guides(color=FALSE) }
  grobs = trees %>% lapply(ggplotGrob) %>% do.call(what=gList)     # convert to grobs
  plot = arrangeGrob(grobs=grobs, nrow=5, ncol=6, bottom=legend)   # then arrage grobs into a matrix and save the final file
  out_file = file.path(fig_dir, paste0(gene_set, "_annotated_trees.pdf"))
  ggsave(out_file, plot, width=20, height=10)
  
  if (out_file == "results/wgs/02-trees/plots/annotated//most_frequent_annotated_trees.pdf") {
    out_file = gsub("[.]pdf", "_small.pdf", out_file)
    ggsave(out_file, plot, width=10, height=8)
  }
}


# overwrite the old trees with the trees including lowpass samples
all_new_files = list.files(tree_file_dir, "_trees.rds", full.names=TRUE)
for (new_file in all_new_files) {
  id = gsub("_.*", "", basename(new_file))
  tree_data[[id]] = readRDS(new_file)$tree %>% MLLPT:::remove_root_tip("GL")
}


for (gene_set in names(driver_data)) {
  
  for (i in seq_along(tree_data)) {
    
    tree = tree_data[[i]]
    data_tree = py_data[[i]]
    
    for (j in which(!tree$tip.label %in% names(data_tree))) {
      id = tree$tip.label[j]
      pnode = tree$edge[tree$edge[,"node"] == j, "parent"]
      snodes = tree$edge[tree$edge[,"parent"] == pnode,"node"]
      sid = na.omit(tree$tip.label[snodes])
      
      sdata = unique(do.call(rbind, data_tree[sid]))
      if (NROW(sdata) > 1) next()
      
      data_tree[[id]] = c(sdata)
    }
    
    
    anc = ancestral.pars(tree, data_tree)
    tree$edge.label = rep("", length(tree$edge.length))
    
    #mt_drivers = match(attr(anc, "id"), names(driver_data$driver_genes))
    var_idx = attr(anc, "index")
    mt_drivers = driver_data[[gene_set]][attr(anc, "id")]
    var_idx = var_idx[order(mt_drivers)]
    mt_drivers = mt_drivers[order(mt_drivers)]
    
    for (j in which(!is.na(mt_drivers))) {
      var_idx_j = var_idx[j] 
      lab_j = mt_drivers[j]
      
      state_nodes = sapply(anc, function(x) if (NROW(x)) x[var_idx_j,2] else NA)
      mutated_nodes = names(which(state_nodes == 1))
      if (length(mutated_nodes) == 0) next()
      
      first_mutated_node = unique(unlist(lapply(mutated_nodes, function(n) {
        na = get_nodes_above(n, tree)
        tail(na[na %in% mutated_nodes], 1)
      })))
      
      stopifnot(length(first_mutated_node) >= 1)
      if(length(first_mutated_node) > 1) warning(paste0(lab_j, " annotated ", length(first_mutated_node), " times in case ", names(tree_data)[i], ".\n"))
      
      for (mn in first_mutated_node) {
        edge_j = c(which(tree$edge[,"node"] == mn), which(tree$tip.label == mn))
        stopifnot(length(edge_j) == 1)
        tree$edge.label[edge_j] = paste0(tree$edge.label[edge_j], "\n", lab_j)
        
        if (length(first_mutated_node) > 1)
          tree$edge.label[edge_j] = paste0(tree$edge.label[edge_j], "!")
      }       
      
    }
    tree$edge.label = gsub("^\n", "", tree$edge.label)
    tree_data[[i]] = tree
  }
  
  trees = tree_data %>% lapply(plot_tree, color_by=status) 
  legend = cowplot::get_legend(trees[[1]])  # keep color legend and remove them from all trees
  for (i in seq_along(trees)) { trees[[i]] = trees[[i]] + guides(color=FALSE) }
  grobs = trees %>% lapply(ggplotGrob) %>% do.call(what=gList)     # convert to grobs
  plot = arrangeGrob(grobs=grobs, nrow=5, ncol=6, bottom=legend)   # then arrage grobs into a matrix and save the final file
  out_file = file.path(fig_dir, paste0(gene_set, "_including_lp_annotated_trees.pdf"))
  ggsave(out_file, plot, width=20, height=15)
}
