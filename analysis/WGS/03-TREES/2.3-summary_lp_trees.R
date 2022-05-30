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

source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# options
dataset_dir = "analysis/WGS/03-TREES/datasets/lp_trees/mp_trees_unfiltered_including_bulks_trees/ml_sample_trees/"
fig_dir = "analysis/WGS/03-TREES/plots/summary_lp_trees_plus_wgs"
tree_file = "analysis/WGS/03-TREES/datasets/mp_trees/mp_trees_unfiltered_including_bulks_trees.rds"
min_purity = 0.05

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Helper functions ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

estimate_lp_correctness = function(x) {
  
  # function that drops lp samples for which deep WGS is available, too
  sample_ids = unlist(x$tip.label[x$tip.label != "GL"])
  annot = suppressWarnings(THmisc::annotation_from_barcode(sample_ids, TRUE))
  node_dists = dist.nodes(x); diag(node_dists) = NA
  tree_height = max(node_dists[which(x$tip.label == "GL"),], na.rm=TRUE)
  
  wgs_ids = annot$tissue_barcode[annot$analyte == "D"]
  wh_lp_with_wgs = annot$analyte == "L" & annot$tissue_barcode %in% wgs_ids
  
  # internal function
  .is_correct = function(i, exact=TRUE) {
    node_id = which(x$tip.label == sample_ids[i])
    node_to_check = phangorn::Ancestors(x, node_id, type="parent")
    desc = phangorn::Descendants(x, node_to_check, type="tips")[[1]]
    str_desc = x$tip.label[desc]
    annot_desc = suppressWarnings(annotation_from_barcode( str_desc[str_desc!="GL"], 1))
    ids_wgs_desc = with(annot_desc, tissue_barcode[analyte_name == "WGS"])
    if (exact) {
      return(all(ids_wgs_desc == annot$tissue_barcode[i]))
    } else {
      if (!any(ids_wgs_desc == annot$tissue_barcode[i])) {
        return(NA)
      } else {
        return(length(ids_wgs_desc))
      }
    }
  }
  
  
  .dist_correct_node = function(i, exact=TRUE) {
    sample_ref = 
      annot %>% 
      filter(tissue_barcode == annot$tissue_barcode[i]) %>% 
      filter(analyte_name == "WGS") %>% 
      dplyr::select(sample_barcode) %>% 
      unlist()
    
    node_id = which(x$tip.label == sample_ids[i])
    parent = phangorn::Ancestors(x, node_id, type="parent")
    target_id = which(x$tip.label == sample_ref)

    node_dists[parent,target_id] / tree_height
  }
  
  
  idx_lp_with_wgs = which(wh_lp_with_wgs)
  correct_edge = sapply(idx_lp_with_wgs, .is_correct)
  correct_linage = sapply(idx_lp_with_wgs, .is_correct, exact = FALSE)
  dist_correct_node = sapply(idx_lp_with_wgs, .dist_correct_node, exact = FALSE)
  dist_correct_node = sapply(dist_correct_node, function(x) ifelse(length(x) == 0, NA, x))
  conf_val = as.numeric(gsub(")", "", gsub(".*Added p = ", "", sample_ids[wh_lp_with_wgs])))
  
  data.frame(
    sample = annot$sample_barcode[wh_lp_with_wgs],
    confidence = conf_val,
    correctly_assigned_edge = correct_edge,
    correctly_assigned_linage = correct_linage,
    dist_correct_node = dist_correct_node,
    number_of_wgs_samples = rep(length(wgs_ids), length(correct_edge))
  )
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load datatasets and functions  ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# load purity estimates from genotpying:
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

trees_original = readRDS(tree_file)$tree
get_annot = function(x) THmisc::annotation_from_barcode(x$tip.label[1])$patient
names(trees_original) = sapply(trees_original, get_annot)

data_lp_assignment=
  list.files(dataset_dir, full.names = TRUE) %>% 
  magrittr::set_names(., gsub("[.]rds","",basename(.))) %>%
  lapply(readRDS)

per_sample_merged = 
  data_lp_assignment %>% 
  lapply(function(x) { x$per_sample_results$init_purity = x$inital_values$purity; return(x)}) %>%
  lapply("[[", "per_sample_results") %>% 
  do.call(what=rbind) %>% 
  cbind(THmisc::annotation_from_barcode(as.character(.$sample)), .) %>% 
  mutate(genotyping_purity = get_purity(sample))

initial_state_vaf = 
  data_lp_assignment %>% 
  sapply(function(x) x$inital_values$vaf_bkgr) %>% 
  unique() %>% 
  data.frame() %>% 
  set_colnames("vaf")

samples_to_drop = per_sample_merged$sample_barcode[per_sample_merged$purity < min_purity]

readr::write_csv(per_sample_merged, file.path(fig_dir, "per_sample_data.csv"))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Summary of purity estimates ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

plot_purity = 
  per_sample_merged %>% 
  ggplot(aes(x=genotyping_purity, y=purity)) + 
  geom_abline(slope=1, intercept = 0, linetype=3) + 
  geom_point(alpha=0.75)  +
  xlim(0, 1) + 
  ylim(0, 1) + 
  theme(legend.position = "bottom") + 
  labs(linetype="") +
  xlab("External purity estimate") + 
  ylab("Purity estimate") 

out_file = file.path(fig_dir, "purity_estimates.pdf")
ggsave(out_file, plot_purity, width=3.3, height=3)


plot_purity_plus_init = 
  plot_purity +  
  geom_point(alpha=0.2, aes(y=init_purity), size=0.75)  + 
  geom_segment(alpha=0.3, aes(xend=genotyping_purity, yend=init_purity), linetype=3, size=0.3)
  
out_file = file.path(fig_dir, "purity_estimates_including_init.pdf")
ggsave(out_file, plot_purity_plus_init, width=3.3, height=3)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Summary of background rates ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

pat_order = 
  with(per_sample_merged, names(sort(tapply(log(background_vaf+1e-4), patient, mean))))

plot_vaf_bkr = 
  per_sample_merged %>% 
  filter(!sample_barcode %in% samples_to_drop) %>% 
  mutate(patient = factor(patient, rev(unique(sort(patient))), ordered = TRUE)) %>% 
  mutate(patient = factor(patient, pat_order, ordered = TRUE)) %>% 
  dplyr::filter(background_vaf > 1e-4) %>% 
  ggplot(aes(x=patient, y=background_vaf+1e-4)) + 
  #geom_violin() + 
  geom_jitter(width=0.2, height=0) + 
  coord_flip() +  
  scale_y_log10(limits=c(min(c(initial_state_vaf$vaf, per_sample_merged$background_vaf+5e-4)), 0.05)) + 
  theme(legend.position = "bottom") + 
  labs(linetype="") +
  ylab("Background VAF") + 
  xlab("") + 
  background_grid() +
  scale_linetype_manual(values=3) + 
  geom_hline(data=initial_state_vaf, aes(yintercept = vaf, linetype="Initial value"))
  
out_file = file.path(fig_dir, "background_vaf.pdf")
ggsave(out_file, plot_vaf_bkr, width=3.3, height=5)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Summary of background rates ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

avg_data = 
  with(per_sample_merged, 
    data.frame(
      purity=tapply(genotyping_purity, patient, mean, na.rm=TRUE),
      vaf=tapply(background_vaf, patient, mean, na.rm=TRUE)
    )
  )


plot_vaf_bkr = 
  avg_data %>% 
  ggplot(aes(x=vaf+1e-4, y=purity)) + 
  geom_point()
  
cor.test(avg_data$purity, avg_data$vaf)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#  Create plot containing all trees ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# get trees and add missing ones
trees_plus_lp = 
  file.path(dataset_dir) %>% 
  list.files(full.names = TRUE) %>% 
  magrittr::set_names(., gsub("_.*", "", basename(.))) %>% 
  lapply(readRDS) %>% 
  lapply("[[", "tree")

wh_missing = sapply(trees_plus_lp[names(trees_original)], is.null)
trees_plus_lp[names(trees_original)[wh_missing]] = trees_original[wh_missing]

.samples_to_drop = function(x) {
  wh = sapply(x$tip.label, function(y) any(sapply(samples_to_drop, grepl, x=y)))
  unlist(x$tip.label)[wh]
}

drop_lp_duplicates = function(x) {
  shrt_bc = gsub("_[LDC][0-9]+( .*$|$)", "", x$tip.label)
  is_lp = grepl("Added", x$tip.label)
  wh_drop = is_lp & shrt_bc %in% shrt_bc[!is_lp]
  drop.tip(x, x$tip.label[wh_drop])
}

# Create plot containing all trees:
tree_plots =
  trees_plus_lp[names(trees_original)] %>% 
  lapply(MLLPT:::remove_root_tip, "GL") %>% 
  lapply(plot_tree, tip_label_size=2, lty_lp = 1, alpha_lp=1) %>% 
  lapply(ggplotGrob)

plots = arrangeGrob(grobs=do.call(gList, tree_plots), nrow=6, ncol=5)
ggsave(file.path(fig_dir, "tree_panel_all.pdf"), plots,  width=9, height=12)


# Create plot containing unique sample trees:
tree_plots =
  trees_plus_lp[names(trees_original)] %>% 
  lapply(function(x) drop.tip(x, .samples_to_drop(x))) %>% 
  lapply(drop_lp_duplicates) %>%
  lapply(MLLPT:::remove_root_tip, "GL") %>% 
  lapply(plot_tree, tip_label_size=2, lty_lp = 1, alpha_lp=1) %>% 
  lapply(ggplotGrob) 
 
plots = arrangeGrob(grobs=do.call(gList, tree_plots), nrow=6, ncol=5)
ggsave(file.path(fig_dir, "tree_panel_unique.pdf"), plots, width=9, height=12)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plot LL ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

plot_nll =
  per_sample_merged %>%
  filter(!sample_barcode %in% samples_to_drop) %>% 
  mutate(patient = factor(patient, rev(unique(sort(patient))), ordered = TRUE)) %>% 
  mutate(patient = factor(patient, names(sort(tapply(log(-mll), patient, mean))), ordered = TRUE)) %>% 
  ggplot(aes(x=patient, y=-mll)) + 
  #geom_violin() + 
  scale_y_continuous(limits=c(0, 8000), oob = scales::squish) +   
  scale_y_log10() +
  geom_jitter(width=0.1, height=0) + 
  coord_flip() +  
  theme(legend.position = "bottom") + 
  labs(linetype="") +
  ylab("Negative Log-Likelihood (NLL)") + 
  xlab("") + 
  background_grid() +
  theme(plot.margin=unit(c(0.2,0.6,0.2,0.2), "cm"))

ggsave(file.path(fig_dir, "nll_all_cases.pdf"), plot_nll, width=3.7, height=5)

