# This script selects the peaks to test/plot between regions, based on the top 20 most recurrent:
# gained promoter, lost promoter, gained enhancer, lost enhancer
# Plus driver mutations with recurrence >=5
source("setup_environment/0-source.R")
library(stringr)
library(tidyr)
library(dplyr)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

fig_dir = "analysis/ATAC/MEGABULKS/plots/megabulk_summary2/"
dat_dir = "analysis/ATAC/CLAIRE/Epigenetic/datasets"

drivers = c(
  THmisc::gene_lists$driver_genes_martincorena_2017, 
  THmisc::gene_lists$`IntOGen-DriverGenes_COREAD`,
  readr::read_tsv("external_datasets/IntOGen-Drivers-20200213.tsv")$SYMBOL
) %>% unique() %>% sort()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

peak_to_gr = function(x) {
  stopifnot(all(grepl("(chr)?:[0-9]+-[0-9]+", x)))
  
  id_spl = strsplit(x, "[:-]")
  
  chr = sapply(id_spl, "[", 1)
  start = as.numeric(sapply(id_spl, "[", 2))
  end = as.numeric(sapply(id_spl, "[", 3))
  
  df = data.frame(chr, start, end, peaks_to_plot=x)
  row.names(df) = NULL
  
  GenomicRanges::makeGRangesFromDataFrame(
      df,
      keep.extra.columns = TRUE,
      ignore.strand = FALSE,
      seqinfo = NULL,
      starts.in.df.are.0based = FALSE
    )
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(dat_dir, FALSE, TRUE)

# get relevant peak data
rec_data = readRDS("analysis/ATAC/MEGABULKS/datasets/edger_summary/final_reccurence_data.rds")
driver_peaks = get_driver_peaks(summary_data = rec_data, genes = drivers, min_rec = 1)
peaks = unique(c(driver_peaks, with(rec_data$summary, peak[recurrence>5])))
peaks_all = with(rec_data$summary, peak)

#making the peaks for subclonal analysis
peaks_gr = peaks %>%  unique() %>% sort() %>% peak_to_gr()
saveRDS(peaks_gr, file.path(dat_dir, "recurrent_peaks.rds"))
peals_all_gr = peaks_all %>%  unique() %>% sort() %>% peak_to_gr()
saveRDS(peals_all_gr, file.path(dat_dir, "all_peaks.rds"))
