source("setup_environment/0-source.R")
source("functions/atac.R")
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(magrittr)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

idx = as.numeric(commandArgs(trailingOnly=TRUE)[1])

fig_dir = "analysis/ATAC/MEGABULKS/plots/adenoma_power_analysis"
outdat_dir = "analysis/ATAC/MEGABULKS/datasets/poweranalysis"

dat_dir = "analysis/ATAC/MEGABULKS/datasets/edger_results_plus_cn"
dat_dir_sum = "analysis/ATAC/MEGABULKS/datasets/edger_summary"
dat_dir_norm = "analysis/ATAC/MEGABULKS/datasets/edger_summary_normal_glands"

input_data_p = file.path(dat_dir, "input_data_proximal_nucleosome_free.rds") %>% readRDS()
input_data_d = file.path(dat_dir, "input_data_distal_nucleosome_free.rds") %>% readRDS()
dge_lists = file.path(dat_dir, "dge_lists_proximal_nucleosome_free.rds") %>% readRDS()
summary_table = file.path(dat_dir_sum, "final_reccurence_data.rds") %>% readRDS()
summary_table_n = file.path(dat_dir_norm, "final_reccurence_data_no_cn.rds") %>% readRDS()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(outdat_dir , FALSE, TRUE)

input_data =
  list(
    fc = rbind(input_data_p$fc, input_data_d$fc),
    cn = rbind(input_data_p$cn, input_data_d$cn),
    counts = rbind(input_data_p$counts, input_data_d$counts),
    groups = input_data_p$groups
  )

# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

n_scaas = apply(cbind(summary_table$sig_matrix, summary_table_n$sig_matrix), 2, sum, na.rm=TRUE)
names(n_scaas) = gsub("normal_vs_", "", names(n_scaas))

rec_scaas = with(summary_table$summary, peak[recurrence > 10])
n_scaas_rec = apply(cbind(summary_table$sig_matrix, summary_table_n$sig_matrix)[rec_scaas,], 2, sum, na.rm=TRUE)
names(n_scaas_rec) = gsub("normal_vs_", "", names(n_scaas_rec))
 

summary_of_input_data =
  with(input_data, data.frame(group = colnames(counts), reads = apply(counts, 2, sum))) %>%
  dplyr::mutate(type=gsub("[0-9]+$", "", gsub(" .*", "", gsub(".*-", "", group)))) %>%
  dplyr::mutate(n_scaa = n_scaas[as.character(group)]) %>%
  dplyr::mutate(n_scaa_rec = n_scaas_rec[as.character(group)])

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

n_reads_per_type = with(summary_of_input_data, tapply(reads, type, median))
subsample_to = n_reads_per_type["adenoma"]
sids_test = with(summary_of_input_data, group[type %in% c("adenoma","pure","normal_gland")])

pltE = readRDS("created_datasets/heatmap_scaas_in_driver_genes_martincorena_and_intogen_all.rds")
pltF = readRDS("created_datasets/heatmap_recurrent_epigenetic_changes_few_incl_drivers.rds")
wh_dupE = duplicated(c(pltE$data$peak))
wh_dupF = duplicated(c(pltF$data$peak))

rec_peaks = 
  split(
    c(as.character(pltE$data$peak)[!wh_dupE], as.character(pltF$data$peak)[!wh_dupF]), 
    paste0(
      c(as.character(pltE$data$event_type_chr)[!wh_dupE], as.character(pltF$data$event_type_chr)[!wh_dupF]), ".",
      c(rep("driver", NROW(pltE$data))[!wh_dupE], rep("rec", NROW(pltF$data))[!wh_dupF]), ".",
      c(as.character(pltE$data$type)[!wh_dupE], as.character(pltF$data$type)[!wh_dupF])
    )
  )

n_reps = 50
if (!is.na(idx))
  sids_test = sids_test[idx]

for (sid in sids_test) {
  
  out_file = file.path(outdat_dir, "per_case", paste0(sid, "_power_analysis_results.rds"))
  if (file.exists(out_file)) next()
  
  subsampling_results = NULL
  set.seed(123)
  
  if (!sid %in% colnames(dge_lists$unadj$counts)) next()
  if (sid %in% subsampling_results$sample) next()
  cat(sid)
  res = do_tests(dge_lists, input_data, sid) 
  
  # input data
  if (!sid %in% colnames(dge_lists$unadj$counts)) next()
  subsample_from = sum(dge_lists$unadj$counts[,sid]) 
  counts = dge_lists$unadj$counts[,colnames(dge_lists$cn_adj)]
  idx = rep_each_by(seq_along(counts[,sid]), counts[,sid])
  
  for (n in seq_len(n_reps)) {
    
    cat(".")

    # subsample data
    sc_to = subsample_to #min(c(subsample_from, round(subsample_to)))
    idx_s = table(sample(idx, sc_to, TRUE))[as.character(seq_len(NROW(counts)))]
    idx_s[is.na(idx_s)] = 0

    # update content of dge list:
    dge_lists_ = dge_lists
    dge_lists_$unadj$counts[,sid] = as.numeric(idx_s)
    dge_lists_$unadj$samples[sid,"lib.size"] = sum(idx_s)
    dge_lists_$unadj %<>% edgeR::calcNormFactors(method="TMMwsp") 

    dge_lists_$cn_adj$counts[,sid] = as.numeric(idx_s)
    dge_lists_$cn_adj$samples[sid,"lib.size"] = sum(idx_s)
    dge_lists_$cn_adj$offset[,sid] = log(input_data$fc[rownames(dge_lists_$cn_adj),sid] * subsample_to)
    res_ = do_tests(dge_lists_, input_data, sid)
   
    n_rec = sapply(rec_peaks, function(x) {
      sum(res[x,"PValue"] < 0.01, na.rm=TRUE)
    }) 
    
    n_rec_sc = sapply(rec_peaks, function(x) {
      sum(res_[x,"PValue"] < 0.01, na.rm=TRUE)
    })  %>% magrittr::set_names(paste0(names(.), "_sc"))
    
    summary_res = cbind(
      data.frame(
        sample = sid,
        i = n,
        type = gsub(" .*", "", gsub(".*-", "", sid)),
        n_fp =  sum(res_$PValue < 0.01 & res$PValue >= 0.01),
        n_fn = sum(res_$PValue >= 0.01 & res$PValue < 0.01),
        n = sum(res[,"PValue"] < 0.01, na.rm=TRUE),
        n_sc = sum(res_[,"PValue"] < 0.01, na.rm=TRUE)
      ), 
      t(data.frame(n_rec)), 
      t(data.frame(n_rec_sc))
    )
    
    subsampling_results = rbind(subsampling_results, summary_res)
  }
  
  dir.create(dirname(out_file), FALSE, TRUE)
  saveRDS(subsampling_results, out_file)
  cat("\n")
}

