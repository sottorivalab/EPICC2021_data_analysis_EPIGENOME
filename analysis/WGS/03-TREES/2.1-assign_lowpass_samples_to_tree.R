source("setup_environment/0-source.R")

#---------------------------------------------
# Options
#---------------------------------------------

idat_dir = "analysis/WGS/03-TREES/datasets/mp_trees/"
odat_dir = "analysis/WGS/03-TREES/datasets"
idx = suppressWarnings(as.numeric(commandArgs(trailingOnly=TRUE)[1]))

min_confidence = 0
vaf_bkgr = 1e-3 # background error rate if mutation is absent 
analytes_to_add = list(lp_trees="L", lp_and_atac_trees=c("L","C"))

samples_to_add = c(# Just low purity
  "EPICC_C536_B1_G4_D1","EPICC_C536_B1_G6_D1",
  "EPICC_C538_A1_G10_D1","EPICC_C539_B1_G6_D1",
  "EPICC_C539_B1_G8_D1","EPICC_C539_D1_G3_D1",
  "EPICC_C539_D1_G8_D1","EPICC_C539_D1_G9_D1",  
  "EPICC_C542_D1_G8_D1","EPICC_C547_B3_L1_D1",
  "EPICC_C549_A1_G8_D1","EPICC_C555_A234_L1_D1",
  "EPICC_C562_B1_G2_D1","EPICC_C538_A1_G10_D1",
  # Moderate contamination with normal crypts
  "EPICC_C528_A1_G9_D1","EPICC_C531_D1_G2_D1",
  "EPICC_C549_C1_G6_D1","EPICC_C550_C1_G7_D1",
  "EPICC_C550_B1_G5_D1","EPICC_C560_A1_G8_D1")

vcf_file_str = file.path(
  "analysis",
  "WGS",
  "02-SNV",
  "datasets",
  "vcf_datasets",
  "%s_GRCh38_platypus_callmode0_mutect2_prior_all_samples_%s_annotated.rds"
)

#---------------------------------------------

ifiles = list.files(idat_dir, "trees[.]rds$", full.names=TRUE, recursive=TRUE)
names(ifiles) = gsub("[.]rds", "", basename(ifiles))

#---------------------------------------------

get_cn_states_of_mutations = function(vcf, smps=NULL) {
  
  m_data = readRDS(vcf) 
  annot = THmisc::annotation_from_barcode(colnames(m_data))     
  if (!is.null(smps)) m_data = m_data[,colnames(m_data) %in% smps]
  smp_pur = get_purity(colnames(m_data))
  
  # determine sites with equal cn in all tumour samples:
  cn_states = geno(m_data)$CN
  wh_equal_cn = apply(cn_states, 1, function(x) length(unique(x)) == 1)
  
  # tabulate estimates:
  cn_and_mm_data =
    data.frame(
      id = rownames(m_data)[wh_equal_cn],
      cn = unlist(c(apply(cn_states[wh_equal_cn,], 1, unique))),
      mm = NA
    ) %>% dplyr::filter(cn >= 1 & cn <= 4 & !is.na(cn))
  
  # vectorized estimation of mm
  for (cn in unique(cn_and_mm_data$cn)) {
    
    # view on data
    wh_cols = cn_and_mm_data$cn == cn
    m_data_ = m_data[cn_and_mm_data$id[wh_cols],]
    
    # get log-lik for each state
    d_array =
      lapply(seq(0, cn), function(mm) {
        # expected vaf for each mutation
        vaf_m = (mm * smp_pur) / (smp_pur * cn + (1 - smp_pur) * 2)
        if (mm == 0) vaf_m = rep(0.01, length(vaf_m))
        lapply(seq_len(NCOL(m_data_)), function(i) {
          N = c(unlist(geno(m_data_)$NR[, i]))
          n = c(unlist(geno(m_data_)$NV[, i]))
          d = dbinom(n, N, vaf_m[i], log = TRUE)
          return(d)
        }) %>% do.call(what = cbind)
      }) %>% abind::abind(along = 3)
    
    
    # drop states where most lik is not mutated at all
    wh_drop = abind::abind(rep(list(apply(d_array, c(1,2), which.max)), cn+1), along = 3)
    d_array[wh_drop == 1] = NA
    d_array = d_array[,,-1,drop=FALSE]
    
    # get most likely estimate of mm:
    mm_est = apply(apply(d_array, c(1, 3), sum, na.rm=TRUE), 1, which.max)
    cn_and_mm_data$mm[wh_cols] = mm_est[cn_and_mm_data$id[wh_cols]]
  }
  
  return(cn_and_mm_data)
}

find_lp_input_files = function(pat, analyte_to_use="L", to_add="") {
  
  sub_dirs = c(
    file.path("SNVS", "atac"),
    file.path("SNVs", "wgs_and_lp")
  )
  
  ifiles = 
    file.path(.pipeline_ddir, sub_dirs) %>%
    list.files("[.]table[.]gz$", full.names=TRUE) %>% 
    grep(pattern=pat, value=TRUE)
  
  if (length(ifiles) == 0) 
    return(NULL)
  
  ifiles %>%
    cbind(file=., THmisc::annotation_from_barcode(basename(.), TRUE)) %>%
    dplyr::filter(patient == pat) %>% 
    magrittr::set_rownames(.$sample_barcode) %>%
    dplyr::filter((!sample_barcode %in% .excluded_samples) | sample_barcode %in% to_add) %>%
    dplyr::filter(analyte %in% analyte_to_use | sample_barcode %in% to_add) %>%
    dplyr::filter(!sample_barcode %in% tree$tip.label) %>%
    dplyr::filter(tissue_type != "normal")
}

get_initial_purity = function(smp_ids_wgs, smp_ids_lp) {
  
  smp_ids_wgs = smp_ids_wgs[smp_ids_wgs != "GL"]
  
  # get annotation of wgs
  annot_smps_wgs = THmisc::annotation_from_barcode(smp_ids_wgs)
  annot_smps_wgs$purity = get_purity(annot_smps_wgs$sample_barcode)
  
  # get annot of lps:
  annot_smps_lp = THmisc::annotation_from_barcode(smp_ids_lp)
  annot_smps_lp$purity = NA
  
  # get purity estimates
  for (i in seq_len(NROW(annot_smps_lp))) {
    
    wh_same_reg = annot_smps_lp$region[i] == annot_smps_wgs$region
    wh_not_same_smp = annot_smps_lp$tissue_barcode[i] == annot_smps_wgs$tissue_barcode
    
    wh_use = wh_not_same_smp & wh_same_reg # same region, but not same sample
    if (sum(wh_use) == 0) wh_use = wh_same_reg # otherwise all other samples
    
    if (sum(wh_use) != 0)
      annot_smps_lp$purity[i] = mean(annot_smps_wgs$purity[wh_use], na.rm=TRUE)
    else 
      annot_smps_lp$purity[i] = 1
  }
  
  magrittr::set_names(annot_smps_lp$purity, annot_smps_lp$sample_barcode)
}

load_mut_data = function(df, ...) {
  
  cat("Loading data:\n")
  
  sample_data = 
    magrittr::set_names(df$file, df$sample_barcode) %>%
    pbapply::pblapply(MLLPT:::load_genotyping_file, ...)
  
  for (i in seq_along(sample_data)) {
    if (!df$analyte_name[i] %in% c("ATAC-seq","WGS")) {
      # for (bad) wgs samples use clonal cn data from wgs (default)
      # otherwise use the actual CN calls
      sid = names(sample_data)[i]
      sample_data[[i]] = sample_data[[i]] %>% 
        dplyr::mutate(copy_number = as.numeric(THmisc::get_cnas(id, .cna_data, sid)))
    }
    sample_data[[i]] = sample_data[[i]] %>% 
      dplyr::filter(!is.na(copy_number)) %>% 
      dplyr::filter(copy_number >= mm)
  }
  
  return(sample_data)
}

#---------------------------------------------

input_combinations = 
  expand.grid(
    i = seq_along(ifiles),
    j = seq_len(30),
    k = seq_along(analytes_to_add)
  )

if (!is.na(idx)) 
  input_combinations = input_combinations[idx, , drop=FALSE]


for (n in seq_len(NROW(input_combinations))) {
  
  i = input_combinations[n,"i"]
  j = input_combinations[n,"j"]
  k = input_combinations[n,"k"]
  
  treeset = readRDS(ifiles[i])
  
  # basic variables
  if(length(treeset$trees) < j) next()
  pat = names(treeset$trees)[j]
  tree = treeset$trees[[j]]
  tree_data = treeset$data[[j]]
  vcf_set = gsub("_.*", "", gsub("mp_trees_", "", names(ifiles)[i]))
  vcf_file = sprintf(vcf_file_str, pat, vcf_set)
  if (!file.exists(vcf_file)) next()
  
  print(pat)
  
  # output dir
  ofile_name = paste0(pat, "_", vcf_set, "_", names(ifiles)[i], ".rds")
  subdir = file.path(odat_dir, names(analytes_to_add)[k], names(ifiles)[i]) 
  out_file =  file.path(subdir, "ml_sample_trees", ofile_name)
  out_file_cn = file.path(subdir, "cn_states", ofile_name)
  if (file.exists(out_file)) next()
      
  dir.create(dirname(out_file), FALSE, TRUE)
  dir.create(dirname(out_file_cn), FALSE, TRUE)
      
  # prepare input data
  cn_state_data = get_cn_states_of_mutations(vcf_file, tree$tip.label)
  saveRDS(cn_state_data, out_file_cn, version = 2)      
  cat("Marginal mutation states:")
  print(table(mm=cn_state_data$mm, cn=cn_state_data$cn))
  
  # lp count file that can be used for patient
  dfiles = find_lp_input_files(pat, analytes_to_add[[k]], samples_to_add)
  if (is.null(dfiles) | NROW(dfiles) == 0) next()
  
  mut_data = load_mut_data(dfiles, cn_data = cn_state_data)
  inital_purity = get_initial_purity(tree$tip.label,  dfiles$sample_barcode)
      
  # add the samples
  tree_plus_lp =
    MLLPT::add_lowpass_sampled(
      tree = tree,
      phydata = tree_data,
      sample_data = mut_data,
      min_confidence = min_confidence,
      vaf_bkgr = vaf_bkgr,
      max_vaf_bkgr = 0.05,
      purity_estimates = inital_purity,
      split_states = TRUE,
      return_details = TRUE, 
      min_edge_length = 0, 
      rescale_tree = TRUE
    )
  
  saveRDS(tree_plus_lp, out_file, version = 2)
}

  
