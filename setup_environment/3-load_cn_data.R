.load_raw_or_dump = function(files, fun, dump_file, ...) {
  if (length(files)) checkmate::assertFileExists(files, access = "r")
  checkmate::assertFunction(fun)
  checkmate::checkPathForOutput(dump_file, overwrite = TRUE)
  
  if (file.exists(dump_file)){
    data = readRDS(dump_file)
  } else {
    if (length(files) > 0) {
      data = fun(files, ...)
      dir.create(dirname(dump_file), showWarnings = FALSE, recursive = TRUE)
      saveRDS(data, dump_file, version = 2)
    } else {
      stop("Couldn't load raw data or dump file.\n")
    }
  }
  
  return(data)
}

.excluded_samples = 
  "created_datasets/excluded_samples.csv" %>% 
  readr::read_csv(progress = FALSE, col_types = readr::cols()) %>% 
  unlist() %>% magrittr::set_names(NULL)

cat("Loading all CNA data.\n")

# wgs ##########################################################################

cat(" -> Loading sequenza fits.\n")
.sequenza_purity_ploidy_fits = 
  .load_raw_or_dump(
    files = .sequenza_fit_files,
    fun = THmisc::load_sequenza_purity_files,
    dump_file = file.path("created_datasets", "sequenza_purity_ploidy_fits.rds"),
    verbose = TRUE
  ) %>% (function(x) x[!names(x) %in% .excluded_samples])


cat(" -> Loading sequenza cn segments.\n")
.sequenza_cna_data = 
  .load_raw_or_dump(
    files = .sequenza_cna_files,
    fun = THmisc::load_sequenza_segment_files,
    dump_file = file.path("created_datasets", "sequenza_cna_data.rds"),
    verbose = TRUE
  ) %>% (function(x) x[!names(x) %in% .excluded_samples])

# low pass #####################################################################

cat(" -> Loading lowpass fits.\n")
.lowpass_ploidy_fits = 
  .load_raw_or_dump(
    files = .lp_metric_files,
    fun = function(x) {
      pbapply::pblapply(x, read.delim) %>% 
        do.call(what = rbind) %>% 
        dplyr::mutate(file=basename(x)) %>% 
        dplyr::mutate(sample=gsub("_500kb_.*","",gsub("^C[[:digit:]]+_","",file)))
    },
    dump_file = file.path("created_datasets", "lowpass_ploidy_fits.rds")
  ) %>% dplyr::filter(!sample %in% .excluded_samples)


.lowpass_ploidy_fits = 
  .lowpass_ploidy_fits %>% 
  dplyr::select(sample, Purity, psit, Solution, PGA) %>%
  magrittr::set_colnames(c("sample","purity","ploidy","fitted","percentage_changed")) %>%
  dplyr::mutate(purity=ifelse(fitted, purity, NA)) %>%
  magrittr::set_rownames(.$sample) 


cat(" -> Loading lp cn segments.\n")
.lowpass_cn_calls = 
  .load_raw_or_dump(
    files = magrittr::set_names(.lp_cn_files, gsub("^C[[:digit:]]*_", "", names(.lp_cn_files))),
    fun = function(x) pbapply::pblapply(x, load_low_pass_custom_call_files),
    dump_file = file.path("created_datasets", "lowpass_cn_calls.rds")
  ) %>% (function(x) x[!names(x) %in% .excluded_samples])


# low pass genotyping data #####################################################

tryCatch({
  cat(" -> Loading genotyping purity estimates.\n")
  if (file.exists(.gt_purity_esimates)) {
    .sample_genotyping_data = 
      readRDS(.gt_purity_esimates) %>% 
      dplyr::select(sample=sample_barcode, estimated_purity) %>% 
      cbind(., THmisc::annotation_from_barcode(.$sample)) %>% 
      dplyr::select(sample_barcode, purity=estimated_purity) %>% 
      magrittr::set_rownames(.$sample_barcode)
  }
}, error=function(e) {print("Couldn't load genotyping purity estimates. Reason is:"); print(e)})

# getter functions #############################################################

.cna_data = c(.sequenza_cna_data, .lowpass_cn_calls)

get_sample_annotation__ = function(samples, attribute, silent=FALSE) {
  
  if (any(duplicated(samples))) {
    
    unique_samples = unique(samples)
    
    annot = get_sample_annotation__(unique_samples, attribute, silent)
    names(annot) = unique_samples
    
    res = annot[samples]
    names(res) = NULL
    
    return(res)
  }
  
  
  results = 
    do.call(rbind, lapply(samples, function(sample) {
      tryCatch({
        
        annotation = tryCatch({THmisc::annotation_from_barcode(sample, TRUE)}, 
                              error=function(e) {return(NULL)})
        
        # invalid annotation
        if (is.null(annotation)) {
          return(data.frame(value = NA, reason = "Invalid barcode."))
        }
        
        # genotyping data only purity
        if (attribute == "purity_gt" ) {
          if (annotation$sample_barcode %in% rownames(.sample_genotyping_data)) {
            value = .sample_genotyping_data[annotation$sample_barcode, "purity"]
            return(data.frame(value = value, reason = NA))
          }
        }
        
        # sequenza fits (WGS)
        if (exists(".sequenza_purity_ploidy_fits")) {
          if (annotation$sample_barcode %in% rownames(.sequenza_purity_ploidy_fits)) {
            value = .sequenza_purity_ploidy_fits[annotation$sample_barcode, attribute]
            return(data.frame(value = value, reason = NA))
          }
        }
        
        # lowpass fits (LP)
        if (annotation$sample_barcode %in% rownames(.lowpass_ploidy_fits)) {
          value = .lowpass_ploidy_fits[annotation$sample_barcode, attribute]
          return(data.frame(value = value, reason = NA))
        }
        
        # throw error, caught below.
        stop("")
        
      }, error=function(e) {
        return(data.frame(value = NA, reason = "Missing data."))
      })
    }))
  
  results = cbind(sample=samples, results)
  
  wh_missing = is.na(results$value)
  
  if (any(wh_missing) & !silent) {
    cat(sprintf("Failed to get %s estimates for the following samples:\n\n", attribute))  
    print(results[wh_missing,c("sample","reason")])
  }
  
  return(results$value)
}

get_purity = function(samples, silent=FALSE) {
  get_sample_annotation__(samples, "purity", silent)
}

get_purity_gt = function(samples, silent=FALSE) {
  get_sample_annotation__(samples, "purity_gt", silent)
}

get_ploidy = function(samples, silent=FALSE) {
  get_sample_annotation__(samples, "ploidy", silent)
}
