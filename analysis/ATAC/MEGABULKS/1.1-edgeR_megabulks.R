source("setup_environment/0-source.R")
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Options ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dat_dir = "analysis/ATAC/MEGABULKS/datasets/edger_results"
dis_source = "common" # source of dispersion estimate

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
annodb = "org.Hs.eg.db"

dataset_file = "analysis/ATAC/MEGABULKS/datasets/counts/data_matrix.rds"
atac_purity_data = "created_datasets/genotyping_estimates_per_sample.rds"
set = "nucleosome_free"

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

data = dataset_file %>% readRDS()
annot_d = THmisc::annotation_from_barcode(colnames(data), TRUE)

gr = with(annot_d, paste0(patient, "_nucleosome_free-", tissue_type))
idx = split(seq_along(gr), gr)
idx[["pan_patient-region_E"]] = annot_d$tissue_type == "normal" & annot_d$sample_type == "G"
data_gr = do.call(cbind, lapply(idx, function(x) apply(data[,x,drop=FALSE], 1, sum)))
data = cbind(data, data_gr)

atacseq_purity = atac_purity_data %>% readRDS()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

cn_table = 
  rownames(data) %>% 
  unlist %>% 
  as("GRanges") %>% "-"(250) %>%
  THmisc::get_cnas(.cna_data, names(.cna_data)) %>% 
  magrittr::set_rownames(rownames(data))

annot = THmisc::annotation_from_barcode(colnames(cn_table))
wh_tumour = annot$tissue_type == "cancer"

avg_cn_peaks = 
  apply(
    X = cn_table[, wh_tumour],
    MARGIN = 1,
    FUN = tapply,
    INDEX = annot$patient[wh_tumour],
    mean, na.rm = TRUE
  ) %>% t()


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# groupings to plot:
pairs_to_plot = list(
  NvA = c("normal", "adenoma"),
  NvC = c("normal", "cancer"),
  AvC = c("adenoma", "cancer"),
  NvCp = c("normal","pure")
)

test_group_labels = c(
  "cancer" = "Cancer",
  "impure" = "Impure cancer",
  "pure" = "Pure cancer",
  "NA" = "Unknown purity",
  "00" = "Purity < 20%",
  "20" = "Purity > 20%",
  "40" = "Purity > 40%",
  "60" = "Purity > 60%"
)


# annotate data:
peak_annotation  = 
  data %>%
  rownames() %>% 
  as("GRanges") %>%
  ChIPseeker::annotatePeak(TxDb=txdb, annoDb=annodb)
    
tss_relationship = 
  peak_annotation@detailGenomicAnnotation$Promoter %>% 
  ifelse("proximal", "distal")
  
for (peak_group in unique(tss_relationship)) {
    
  # create dge list:``
  peak_chr = as.character(seqnames(peak_annotation@anno))
  wh_use = tss_relationship == peak_group & peak_chr %in% paste0("chr", 1:22)
  matrix = data[wh_use,]
  colnames(matrix) = gsub("_nucleosome_free", "", colnames(matrix))
  
  # split into pure and impure regions
  wh_sample = grepl("EPICC", colnames(matrix))
  annot = THmisc::annotation_from_barcode(colnames(matrix)[wh_sample], T)
  purity_per_sample = atacseq_purity[annot$sample_barcode,"estimated_purity"]
  
  cut1 = # first cut, pure and impure
    purity_per_sample %>% 
    cut(c(0,0.4,1), include.lowest = TRUE) %>% 
    factor(labels = c("impure","pure")) %>% 
    as.character() %>% 
    magrittr::set_names(colnames(matrix)[wh_sample])
  
  cut2 = # second cut, 4 groups in 0.25 intervals
    purity_per_sample %>% 
    cut(c(0,0.2,0.4,0.6,1), include.lowest = TRUE) %>% 
    factor(labels = c("00","20","40","60")) %>% 
    as.character() %>% 
    magrittr::set_names(colnames(matrix)[wh_sample])
  
  wh_na = annot$tissue_type != "cancer" #is.na(purity_per_sample) |
  cuts = c(cut1[!wh_na], cut2[!wh_na])
  split_label = paste0(rep(annot$patient[!wh_na], 2), ".", cuts)
  purity_groups = split(cuts, split_label)
  
  matrix_purity_groups = 
    do.call(cbind, lapply(purity_groups, function(gr) {
      if (length(gr) == 0) return(NULL)
      apply(matrix[,names(gr),drop=FALSE], 1, sum)
    }))
  
  # construct a per adenoma column
  wh_sample = grepl("EPICC", colnames(matrix))
  annot = THmisc::annotation_from_barcode(colnames(matrix)[wh_sample], T)
  wh_adeno = annot$tissue_type == "adenoma"
  tmp = paste0(annot$tumour_id[wh_adeno], "-per_tumour")
  names(tmp) = colnames(matrix)[wh_sample][wh_adeno]
  adenoma_groups = split(tmp, tmp)
  
  matrix_adenoma_groups = 
    do.call(cbind, lapply(adenoma_groups, function(gr) {
      if (length(gr) == 0) return(NULL)
      apply(matrix[,names(gr),drop=FALSE], 1, sum)
    }))

  # exclude per region measurements
  wh_cols = grepl("(tumour|cancer|adenoma|pan_patient-region_E)", colnames(matrix))
  matrix_megabulks = matrix[,wh_cols]
  colnames(matrix_megabulks) = gsub("tumour", "cancer", colnames(matrix_megabulks))
  
  # normalise and estimate dispersion:
  wh = grepl("C5.*-cancer", colnames(matrix_megabulks))
  gr = rep(1, ncol(matrix_megabulks[,wh]))
  dge_list_megabulks = edgeR::DGEList(matrix_megabulks[,wh], group=gr)
  
  dge_list_megabulks = dge_list_megabulks %>% 
    edgeR::calcNormFactors(method="TMMwsp") %>% 
    edgeR::estimateDisp()  %>% 
    (function(x) {x$samples$group = colnames(matrix_megabulks[,wh]); x})
  
  # get cpm values:
  matrix_relevant = cbind(matrix_megabulks, matrix_purity_groups, matrix_adenoma_groups)
  dge_list = edgeR::DGEList(cbind(matrix_relevant), group=rep(1, ncol(matrix_relevant)))
  dge_list = dge_list %>% edgeR::calcNormFactors(method="TMMwsp") 
  dge_list$samples$group=rownames(dge_list$samples)
  
  for (cn in names(dge_list_megabulks)) {
    if (cn %in% c("counts","samples")) next()
    dge_list[[cn]] = dge_list_megabulks[[cn]]
  }
  
  rm(dge_list_megabulks)
  
  cpm_values = dge_list %>% 
    edgeR::cpm(normalized.lib.sizes=TRUE, prior.count=0.5) %>% 
    data.frame()
  
  max_cpm_per_peak = apply(cpm_values, 1, max)
  
  # grouping of samples
  cases = sapply(strsplit(colnames(cpm_values), "[.]"), "[", 1) 
  type = sapply(strsplit(colnames(cpm_values), "[.]"), "[", 2) 
  type[colnames(cpm_values) == "pan_patient.region_E"] = "normal"
  
  for (case in grep("C5", unique(cases), value = TRUE)) {
    
    print(case)
    
    if (case %in% c("C001","pan_patient")) next() 
    
    # find original column names:
    ids_not_n = unlist(pairs_to_plot) %>% unique() %>% (function(x) {x[x!="normal"]})
    wh_cols = c(colnames(cpm_values)[type == "normal"], paste0(case, ".", ids_not_n))
    names(wh_cols) = c("normal", ids_not_n)
    
    # subset data to current set:
    wh = wh_cols %in% colnames(cpm_values)
    data_currente_case = data.frame(cpm_values[,wh_cols[wh]])
    colnames(data_currente_case) = names(wh_cols)[wh]
    data_currente_case$cn_cancer = avg_cn_peaks[rownames(data_currente_case),case]
    
    # do statistical tests:
    test_results = 
      lapply(pairs_to_plot, function(pair) {
        pair_name = gsub("[.]","-", wh_cols[pair])
        #pair_name = wh_cols[pair]
        if (!all(pair_name %in% colnames(dge_list))) return(NULL)
        edgeR::exactTest(dge_list, pair_name, dispersion=dis_source)$table %>% 
          mutate(p_adj = p.adjust(PValue, "fdr")) %>% 
          mutate(sig=p_adj<0.1) %>% 
          magrittr::set_rownames(rownames(dge_list)) %>% 
          cbind(cpm_values[rownames(.),wh_cols[pair]], .)
      })
    
    data_currente_case$highlight = label_data(test_results[c("NvA","NvCp")])
    
    
    # save results:
    result_set = 
      list(
        cpm_data = data_currente_case, 
        test_results = test_results
      )
    
    out_file = 
      file.path(
        dat_dir,
        set,
        peak_group,
        paste0(set, "-", case, "-", peak_group, ".rds")
      )
    
    dir.create(dirname(out_file), FALSE, TRUE)
    saveRDS(result_set, out_file)
    
  }
}
