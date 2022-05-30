library(cowplot)
library(dplyr)
library(cowplot)
library(ggpubr)
library(grid)
library(THmisc)

source("setup_environment/0-source.R")

#---------------------------------------------
# Options
#---------------------------------------------

fig_dir = "analysis/WGS/01-CNA/plots"
dataset_dir = "analysis/WGS/01-CNA/datasets"
plot_logr = FALSE
include_low_pass = TRUE
include_blacklisted_samples = FALSE

#---------------------------------------------


# prepare output dirs:
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(dataset_dir, showWarnings=FALSE, recursive=TRUE)


# wgs and lp cn datasets: 
if (plot_logr) {
  columns_keep = c("chromosome","start.pos","end.pos","depth.ratio","L1")
} else {
  columns_keep = c("chromosome","start.pos","end.pos","CNt","L1")
}
  

# load wgs data
wgs_data = reshape::melt(.sequenza_cna_data, measure.vars=c())[,columns_keep] # short to long


# load low-pass data:
if (include_low_pass) {
  if (plot_logr) {
    low_pass_data = reshape::melt(.lowpass_segment_data, measure.vars = c()) # short to long
  } else {
    low_pass_data = reshape::melt(.lowpass_cn_calls, measure.vars = c()) # short to long
  }
  low_pass_data = low_pass_data[, columns_keep]
} else {
  low_pass_data = NULL
}


# bind:
cna_data = rbind(low_pass_data, wgs_data)
ord_chrs = paste0("chr", c(1:22, "X", "Y"))
cna_data$chromosome = factor(cna_data$chromosome, ord_chrs, ordered = TRUE) # insert chromsome order


# annotation of all samples:
annotation = annotation_from_barcode(unique(as.character(cna_data$L1)), TRUE)
rownames(annotation) = unique(as.character(cna_data$L1))
cna_data = cbind(cna_data, annotation[as.character(cna_data$L1), ])
cna_data$L1 = NULL

# remove normals:
cna_data = filter(cna_data, tissue_type != "normal")


# plotting:
for (small in c(FALSE,TRUE)){
  
  plot_samples =
    lapply(sort(unique(cna_data$patient)), function(pat) {
      
      data_plot = 
        filter(cna_data, patient == pat) %>%
        filter(chromosome %in% c(paste0("chr", 1:22))) %>%
        mutate(sample=gsub("EPICC_C[[:digit:]]+_", "", sample_barcode)) 
      
      if (plot_logr) {
        data_plot$fill_value = data_plot$depth.ratio
      } else {
        data_plot$fill_value = 
          factor(
            data_plot$CNt, 
            names(THmisc:::color_cnas), 
            ordered =TRUE
          )
      }
      
      p =
        data_plot %>%
        ggplot(aes(x=(start.pos+end.pos)/2, y=sample, width=(end.pos-start.pos), fill=fill_value)) +
        facet_grid(region~chromosome, space="free", scales="free", as.table=TRUE) +
        geom_tile() + 
        xlab("Genomic position") +
        ylab(pat) +
        theme(
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.text.y = element_text(angle = 0, colour = "gray90"),
          strip.background.y = element_rect(fill = "gray5"),
          axis.text.y = element_text(size = 8),
          strip.text.x=element_blank(), 
          strip.background.x=element_blank()
        )
      
      if (plot_logr) {
        p = p + 
          scale_fill_gradient2(
            low = "#2166ac",
            high = "#b2182b",
            midpoint = (1),
            limits = (c(0.5, 1.5)),
            oob = scales::squish
          ) + 
          guides(fill=FALSE)
      } else {
        p = p + scale_fill_manual(values=THmisc:::color_cnas) + guides(fill=FALSE)
      }
      
      
      if (small) {
        p = p + 
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          theme(strip.text.y=element_text(angle=0, colour="gray5", size=8)) +
          theme(panel.spacing.y=unit(0.5, "pt")) + 
          theme(axis.title.y=element_text(size=10, angle=0)) + 
          theme(strip.text.y=element_blank(), strip.background.y=element_blank()) +
          theme(plot.margin = unit(c(0,2,-1,2), "mm"))
      }
      return(p)
    })
  
  
  cytobands =
    D3GB::GRCh38.bands %>%
    mutate(chr=paste0("chr", chr)) %>%
    filter(chr %in% c(paste0("chr", 1:22))) %>%
    mutate(chr=factor(chr, ord_chrs, ordered=TRUE)) %>%
    ggplot(aes((start+end)/2, y=1, height=3, width=(end-start), fill=score)) +
    facet_grid(1~chr, space="free", scales="free", as.table=TRUE) +
    geom_tile() +  
    scale_fill_manual(
      values = THmisc:::cytoband_colors,
      breaks = names(THmisc:::cytoband_colors)
    ) +
    guides(fill=FALSE) +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      strip.text.x=element_text(angle=90, size=10)
    )
  
  if (small) {
    cytobands = cytobands + theme(strip.text.x=element_blank()) + 
    theme(plot.margin = unit(c(3,2,0,2), "mm"))
  }
  
  
  # arrange and save plots:
  gs = lapply(plot_samples, ggplotGrob) 
  gc = ggplotGrob(cytobands)
  gc$widths = gs[[1]]$widths
  
  if (small) {
    g = gridExtra::arrangeGrob(grobs=c(list(gc), gs), ncol=1, heights=c(1,rep(1,length(gs))))
  } else {
    g = gridExtra::arrangeGrob(grobs=c(list(gc), gs), ncol=1, heights=c(0.3,rep(1,length(gs))))
  }

  if (small) {
    if (plot_logr) {
      fig = file.path(fig_dir, "wgs_and_lowpass_cna_plot_small.png")
    } else {
      fig = file.path(fig_dir, "wgs_and_lowpass_cna_plot_small_cn_states.png")
    }
    ggsave(fig, g, width=6, height=5, limitsize = FALSE)
  } else {
    if (plot_logr) {
      fig = file.path(fig_dir, "wgs_and_lowpass_cna_plot.png")
    } else {
      fig = file.path(fig_dir, "wgs_and_lowpass_cna_plot_cn_states.png")
    }
    ggsave(fig, g, width=10, height=84, limitsize = FALSE)
  }
}
 

if (plot_logr) {
  out_file = file.path(dataset_dir, "logr_ratios.rds")
  saveRDS(cna_data, out_file)
} else {
  out_file = file.path(dataset_dir, "cn_calls.rds")
  saveRDS(cna_data, out_file)
}

