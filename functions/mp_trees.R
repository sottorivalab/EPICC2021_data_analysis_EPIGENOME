
rds_input_to_phyDat = 
  function(rds, ccf_cutoff = 0.25, excluded_samples=c(), exclude_bulks=FALSE, equal_cna_states=FALSE, exclude_normals=FALSE, tissue=NULL, cn_equal_tissue=NULL) {
    
    # Load and annotate data:
    vcf_data = readRDS(rds)
    annotation = THmisc::annotation_from_barcode(colnames(vcf_data))
    wh_drop = rep(FALSE, NCOL(vcf_data))

    # exclude bulks:
    if (exclude_bulks) {
      wh_drop = wh_drop | annotation$sample_type == "B"
    }
    
    # exclude normals:
    if (exclude_normals) {
      wh_drop = wh_drop | annotation$tissue_type == "normal"
    }
    
    # keep only selected tissues.
    if (!is.null(tissue)) {
      wh_drop = wh_drop | !annotation$tissue_type %in% tissue
    }
    
    # exclude selected samples:
    if (length(excluded_samples) > 0) {
      wh_drop = wh_drop | colnames(vcf_data) %in% excluded_samples
    }
    
    
    # create mutation matrix
    mutated = data.frame(geno(vcf_data[,!wh_drop,drop=FALSE])$CCF) > ccf_cutoff
    storage.mode(mutated) = "numeric"
    
    # find and remove missing entries: 
    whC = apply(is.na(mutated), 2, all)
    whR = apply(is.na(mutated[,!whC, drop=FALSE]), 1, any)
    mutated = mutated[!whR,!whC, drop=FALSE]

    if (NCOL(mutated) == 0 | NROW(mutated) == 0) {
      return(NULL)
    }
    
    # keep only those without subclonal cn changes:
    if (equal_cna_states) {
      annotation = THmisc::annotation_from_barcode(colnames(mutated))
      ab_states = geno(vcf_data)$AB[rownames(mutated),colnames(mutated)]
      
      if (!is.null(cn_equal_tissue)) {
        wh_use = annotation$tissue_type %in% cn_equal_tissue
        ab_states = ab_states[,wh_use]
      }
      
      equal_states = apply(ab_states, 1, function(x) length(unique(x)) == 1)
      cat("Equal CN states:", mean(equal_states), "\n")
      samples = paste0(colnames(ab_states), collapse=", ")
      
      if (sum(equal_states) < 3) return(NULL)
      mutated = mutated[equal_states,]
    }
    
    
    # convert to phyDat object:
    pyData = 
      mutated %>%
      data.frame() %>%
      mutate(GL=0) %>% 
      t() %>% as.matrix() %>%
      phangorn::phyDat(type="USER", levels=c(0,1))
    
    attr(pyData, "id") = rownames(mutated)
    
    return(pyData)
  }

phyDat_to_tree = function(phyDat, all_trees=FALSE, ...) {
  phyDat %>%
    phangorn::pratchet(trace=FALSE, all=all_trees, ...) %>%
    ape::unroot() %>%
    ape::root(out="GL", resolve.root=TRUE) %>%
    phangorn::acctran(phyDat)
}
