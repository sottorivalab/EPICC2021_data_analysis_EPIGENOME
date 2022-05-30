# This script gets the per-sample read counts for the recurrent peaks
library(GenomicFeatures)
library(GenomicRanges)
library(stringr)
source("setup_environment/0-source.R")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

fig_dir = "analysis/ATAC/CLAIRE/Epigenetic/plots"
dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"
idat_dir = "analysis/ATAC/CLAIRE/Epigenetic/input_data"
regex_cutsites = ".*[1-9]_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"
labels = c(clonal = "Putatively Clonal", subclonal = "Putatively Subclonal")
pat_ids = unique(gsub("_.*", "", list.files(file.path(dat_dir, "fits"))))
col = RColorBrewer::brewer.pal(name="Set1", n=9)
coldata_region = readRDS(file.path(idat_dir, "coldata_region.rds"))
combos = rownames(coldata_region)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

source("analysis/ATAC/CLAIRE/Epigenetic/functions/horizon_mass_plot_function_commonE_CL6.R")
source("analysis/ATAC/CLAIRE/Epigenetic/functions/horizon_plot_functions_commonE_CL6.R")

peaks_original = readRDS("analysis/ATAC/MEGABULKS/datasets/edger_summary/final_reccurence_data.rds")$summary
peak_tracks = readRDS(file.path(dat_dir, "peak_tracks.rds")) 

results = readRDS(file.path(dat_dir, "all_deseq_results.rds"))
peaks_original = peaks_original %>% dplyr::filter(peak %in% results$peak)

rec_data = readRDS("analysis/ATAC/MEGABULKS/datasets/edger_summary/final_reccurence_data.rds")
sig_filter = rec_data$sig_matrix[results$peak,]
fc = rec_data$fc[results$peak,]
p = rec_data$p[results$peak,]
p_adj = rec_data$p_adj[results$peak,]

cov = file.path(dat_dir, "cpm_data") %>% 
  list.files("[.]rds$", full.names = TRUE) %>% 
  magrittr::set_names(gsub("[.]rds", "", basename(.))) %>% 
  lapply(readRDS)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

plot_set_mass = function(file, peaks, patients, ...) {
  
  pdf(file, height = 5 * NROW(peaks), width = 7 * length(patients))
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(NROW(peaks), length(patients))))
  
  for (pat in 1:length(patients)) {
    for (b in 1:NROW(peaks)) {
      grid::pushViewport(grid::viewport(layout.pos.col = pat, layout.pos.row = b) )
      plt = plot_horizon_mass(cov, patients[pat], peaks[b, ]$peak, type = "l", ...)
      print(plt, vp = vplayout(x = pat, y = b))
      grid::popViewport()
    }
  }
  
  dev.off()
}

plot_set_mass = function(file, peaks, patients, ...) {
  
  pdf(file, height = 5 * NROW(peaks), width = 7 * length(patients))
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(NROW(peaks), length(patients))))
  
  for (pat in 1:length(patients)) {
    for (b in 1:NROW(peaks)) {
      grid::pushViewport(grid::viewport(layout.pos.col = pat, layout.pos.row = b) )
      plt = plot_horizon(cov, patients[pat], peaks[b, ]$peak, type = "l", ...)
      print(plt, vp = vplayout(x = pat, y = b))
      grid::popViewport()
    }
  }
  
  dev.off()
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(fig_dir, FALSE, TRUE)

for (event in c("gain", "loss")) {
  for (type in c("promoter", "enhancer")) {
    wh = peaks_original$event_type==event & peaks_original$type==type
    peaks = peaks_original[wh,]
    fname = paste0(type, "_", event, "_suppl.pdf")
    ofile = file.path(fig_dir, fname)
    #plot_set(ofile, peaks, pat_ids)
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Plot subclonal
for(i in 1:length(patients)){
  peakss <- peaks
  patient <- patients[i]
  peaks_with_peaks<-mcols(peaks[peaks$peak %in% unique(unlist(peak_matches[grepl(patient,names(peak_matches))])),])$peak
  for(t in 1:length(peakss)){
    peak <- peakss$peak[t]
    present <-ifelse(peak %in% peaks_with_peaks, "gotpeaks","nopeaks")
    logfc <-round(fc[peak,paste0(patient,".pure")],digits = 2)
    sig <- sig_filter[peak,paste0(patient,".pure")]
    status<-get_status2(patient=patient,
                        peak=peakss$peak[t])
    try(plot_info(pdf_name=paste0("~/subclonal_redo/plots2/",patients[i],"_",
                                  peakss$id[t],"_",
                                  peakss$type[t],"_",
                                  peakss$event_type_chr[t],"_",
                                  status,"_",logfc,"_",sig,"_",present,"_line.pdf"),
                  patient=patient,
                  peak=peakss$peak[t],
                  type="l"))
    graphics.off()
  }}


plot_horizon(cov, "C518", peaks_original$peak)
