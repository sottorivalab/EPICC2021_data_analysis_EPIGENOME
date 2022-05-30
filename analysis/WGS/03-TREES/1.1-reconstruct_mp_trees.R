source("setup_environment/0-source.R")
library(VariantAnnotation)

#---------------------------------------------
# Options
#---------------------------------------------

idat_dir = "analysis/WGS/02-SNV/datasets/vcf_datasets"
dat_dir = "analysis/WGS/03-TREES/datasets/mp_trees"

#---------------------------------------------

dir.create(dat_dir, showWarnings=FALSE, recursive=TRUE)

ifiles = list.files(idat_dir, "[.]rds$", full.names = TRUE, recursive = TRUE)
names(ifiles) = gsub("[.]rds$", "", basename(ifiles))

#---------------------------------------------

drop_bulks =
  function(x) {
    stopifnot("phyDat" %in% class(x))
    not_gl = names(x) != "GL"
    annot = THmisc::annotation_from_barcode(names(x)[not_gl], extract = TRUE)
    wh_keep = !not_gl
    wh_keep[not_gl] = annot$sample_type != "B"
    return(x[wh_keep])
  }

#---------------------------------------------

# find all input rds files:
cat("Loading input data:\n")
pyData =
  pbapply::pblapply(
    ifiles,
    rds_input_to_phyDat,
    ccf_cutoff = 0.25,
    equal_cna_states = FALSE,
    exclude_normals = TRUE
  ) %>% (function(x) x[!sapply(x, is.null)]) # drop cases with no samples/remaining variants


# estimate MP trees
cat("Reconstructing trees:\n")
for (removed_bulks in c(TRUE, FALSE)) {
  for (sample_filter in c("_filtered_", "_unfiltered_")) {
    
    pyData_ = pyData[grepl(sample_filter, names(pyData))]
    names(pyData_) = gsub("_.*", "", names(pyData_))
    if (removed_bulks) pyData_ = lapply(pyData_, drop_bulks)
    trees = pbapply::pblapply(pyData_, phyDat_to_tree, k=100, maxit=1e6)
    
    # save output
    file_name = 
      paste0(
        "mp_trees",
        sample_filter,
        ifelse(removed_bulks, "excluding_", "including_"),
        "bulks_trees.rds"
      )
    
    results = list(trees=trees, data=pyData_)
    saveRDS(results, file.path(dat_dir, file_name))
  }
}

