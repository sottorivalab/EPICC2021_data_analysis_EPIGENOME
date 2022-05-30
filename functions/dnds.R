parse_rds_file_for_dndscv = function(file) {
  
  tryCatch({
    # load vcf data:
    print(file)
    vcf_data = readRDS(file)
    annot = THmisc::annotation_from_barcode(colnames(vcf_data))
    
    # remove excluded and normal samples:
    is_excluded = annot$sample_barcode %in% .excluded_samples
    is_normal = annot$tissue_type == "normal"
    wh_keep = !(is_excluded | is_normal)
    vcf_data = vcf_data[, wh_keep,drop=FALSE]
    annot = annot[wh_keep, ,drop=FALSE]
    
    # exclude incomplete mutations (i.e. NA entries) and non snvs:
    CCF = data.frame(geno(vcf_data)$CCF)
    wh_keep = apply(!is.na(CCF), 1, all) #& isSNV(vcf_data)
    vcf_data = vcf_data[wh_keep,,drop=FALSE]
    CCF = CCF[wh_keep,,drop=FALSE]
    
    if (NROW(CCF) == 0) return(NULL)
    
    tissue_types = as.character(unique(annot$tumour_id))
    
    data = 
      lapply(tissue_types, function(tissue) {
        
        wh_type = annot$tumour_id == tissue
        
        if (grepl("cancer", tissue)) { # let cancer absorb all variants that are shared between tumour and adenomas/normals
          wh_vars = rep(TRUE, NROW(CCF))
        } else {
          wh_cancer = grepl("cancer", annot$tumour_id)
          wh_vars = apply(CCF[,wh_cancer,drop=FALSE] > ccf_cutoff, 1, sum) == 0
        }
        
        # clonal status:
        n_cases = sum(wh_type)
        n_mut = apply(CCF[wh_vars, wh_type,drop=FALSE] > ccf_cutoff, 1, sum)
        
        status = case_when(
          n_mut == sum(wh_type) ~ "clonal",
          n_mut > 1 ~ "subclonal",
          n_mut == 1 ~ "subclonal"
        )
        
        # data for dnds
        data.frame(
          sampleID = gsub("[.]rds$", "", basename(file)),
          tissue = tissue,
          chr = as.character(seqnames(vcf_data)[wh_vars]),
          pos = start(vcf_data)[wh_vars],
          ref = as.character(ref(vcf_data))[wh_vars],
          alt = as.character(unlist(alt(vcf_data)))[wh_vars],
          status = status,
          samples_mutated=n_mut,
          samples_total=n_cases,
          stringsAsFactors = FALSE
        )
      })
    
    data %>%
      do.call(what=rbind) %>% 
      dplyr::filter(!is.na(status)) %>%
      THmisc::lift_hg38_to_hg19()
    
  }, error=function(e) {print(e); return(NULL)})
}

dndscv_wrapper = 
  function(data, ...) {
    tryCatch({
      dplyr::select(data,sampleID,chr,pos,ref,alt) %>%
        dplyr::mutate(chr=gsub("chr","",chr)) %>%
        dndscv(...)
    }, error=function(e) { warning(e); return(NULL) })
  }



load_n_variant_data_dnds = function(data, genes) {
    lapply(data, function(d) {
    if (is.character(d) & file.exists(d)) d = readRDS(d)
    lapply(genes, function(g) {
      if (sum(g %in% d$annotmuts$gene) == 0) 
        return(data.frame(type="Total", Freq=0))
      d$annotmuts %>% 
        dplyr::filter(gene %in% g) %>% 
        dplyr::select(impact) %>% unlist() %>% 
        table() %>% data.frame() %>% 
        magrittr::set_names(c("type","Freq")) %>% 
        (function(x) rbind(x, data.frame(type="Total", Freq=sum(x$Freq))))
    })
  })
}

plot_n_variants_per_type_group = function(n_variants) {
  
  reshape2::melt(n_variants) %>% 
    (function(x) cbind(x, data.frame(do.call(rbind, strsplit(x$L1, "[.]"))))) %>% 
    #dplyr::filter(gsub("[.]MS[IS]$", "", L2) %in% c("adenoma.all","cancer.clonal","cancer.subclonal","cancer.all")) %>% 
    dplyr::mutate(L2=factor(L2, rev(names(gene_set_labels)), rev(gene_set_labels), ordered=1)) %>%
    ggplot(aes(x=type, y=L2, fill=value, label=value)) + 
      geom_tile() + 
      geom_text() + 
      facet_wrap(X1 ~ paste0(X2, " - ", X3, " - ", X4), ncol=2) + 
      scale_fill_distiller(trans="log", direction = 1, breaks=10^(0:6), limits=c(1,1e5)) + 
      xlab("") + ylab("") + labs(fill="N") + 
      theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  
}


curate_data_for_dndscv = function(x) {
    ids = with(x, paste0(tissue, "-", chr, ":", pos, "_", ref, "/", alt))
    wh_not = ids %in% .variants_drop
    x = x[!wh_not,]
    ids = ids[!wh_not]
    x$status[ids %in% .variants_clonal] = "clonal" 
    return(x)
}

load_genset_result_file = function(f) {
  
  cns = c(
    "geneset",
    "filter_stat",
    "tissue",
    "clonal_status",
    "msi_status")
  
  annot = basename(f) %>% 
    gsub(pattern = "[.]rds", replacement = "") %>% 
    strsplit("[.]") %>% 
    do.call(what=rbind) %>% 
    data.frame %>% 
    magrittr::set_colnames(cns)
  
  
  d = readRDS(f) %>% 
    dplyr::mutate(name = gsub("_geneset", "", rownames(.))) %>% 
    dplyr::mutate(type_name = factor(name, names(dnds_type_names), dnds_type_names)) %>% 
    cbind(., annot)
}


plot_dnds_results = function(d, title="dndncv", width_tick=0.15) {
  
  d_mod = d %>%
    dplyr::filter(!is.na(variant_group)) %>% 
    dplyr::mutate(variant_group = droplevels(variant_group)) %>% 
    dplyr::filter(!is.na(gene_set_label))  %>% 
    dplyr::mutate(gene_set_label = droplevels(gene_set_label)) %>%
    dplyr::mutate(x_offset = as.numeric(gene_set_label)) %>%
    dplyr::mutate(x_offset = ((x_offset - 1) - (length(levels(gene_set_label)) - 1) / 2) / (length(levels(gene_set_label))) / 2) %>%
    dplyr::mutate(x_pos = as.numeric(variant_group)) %>%
    dplyr::mutate(alpha_bar=ifelse(cilow > 0 & !is.infinite(cihigh), 1, 0.4)) %>% 
    dplyr::mutate(width_bar=ifelse(cilow > 0 & !is.infinite(cihigh), width_tick, 0)) %>% 
    dplyr::mutate(cihigh = ifelse(cihigh > 200, Inf, cihigh)) %>%
    split(., paste0(.$L1, ".", .$type_name)) %>% 
    do.call(what=rbind) 
  
  
  plot_dnds = 
    d_mod %>%
    ggplot(aes(y=mle, ymin=cilow, ymax=cihigh, group=gene_set_label, color=gene_set_label, x=x_pos+x_offset))+
    geom_hline(yintercept=1, alpha=0.5) +
    geom_errorbar(aes(width=width_bar), size=1) +
    #geom_point(shape=16, size=2) +  
    geom_point(color="black", shape=16, size=1.5) +  
    ylab("dN/dS") + xlab("") + #xlab("Consequence") +
    labs(color="Legend") +
    scale_color_brewer(palette="Set1") + 
    scale_x_continuous(breaks=seq_along(levels(d_mod$variant_group)), labels=levels(d_mod$variant_group)) + 
    labs(color="Genes", alpha="Genes") + 
    cowplot::background_grid() + 
    scale_y_log10(limits=c(0.01, 500)) + 
    ggtitle(title) + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
    theme(strip.text.y=element_text(angle=0)) + 
    theme(strip.background.y=element_rect(fill="gray99", color="white")) + 
    coord_cartesian(clip = 'off') + 
    scale_alpha_identity()
  
  if ("msi_status" %in% colnames(d)) {
    plot_dnds = plot_dnds + facet_grid(msi_status~type_name, scales="free_y")
  } else {
    plot_dnds = plot_dnds + facet_grid(.~type_name, scales="free_y")
  }
  
  return(plot_dnds)
}

plot_variant_heatmap_dndscv = function(data, genes, removed_patients=character()) {
  
  impact_labels = 
    c(Essential_Splice = "Nonsense & Splice Site",
      Nonsense = "Nonsense & Splice Site",
      Missense = "Missense")
  
  tissue_ids = sort(unique(data$tissue))
  data = data %>% dplyr::filter(gene %in% genes)
  
  wh = data$impact %in% names(impact_labels)
  n_mut_gene = with(data[wh,], tapply(samples_mutated > 0, gene, sum))
  order_genes = rev(unique(c(names(sort(n_mut_gene, decreasing = TRUE)), genes)))
  
  mut_mat = with(data[wh,], tapply(samples_mutated/samples_total, list(gene, tissue), max))
  mut_mat = mut_mat[rev(order_genes[order_genes %in% rownames(mut_mat)]),]
  mut_mat[is.na(mut_mat)] = 0
  
  ord = rev(do.call(order, split(mut_mat, seq_len(NROW(mut_mat)))))
  order_samples = unique(c(colnames(mut_mat)[ord], tissue_ids))
  
  # create group data
  tissues = unique(tissue_ids)
  pats = gsub("[.].*", "", tissues) 
  is_msi_cancer = pats %in% .msi_positiv & grepl("cancer", tissues)
  is_msi_adenoma = pats %in% .msi_positiv_adenoma & grepl("adenoma", tissues)
  is_msi = is_msi_cancer & is_msi_adenoma
  groups = paste0(ifelse(is_msi, "MSI", "MSS"), "\n", Hmisc::capitalize(gsub(" .*", "", gsub(".*[.]", "", tissues))))
  names(groups) = tissues
  

  # modify data 
  data_plot = data %>%
    dplyr::mutate(patient = gsub("[.].*", "", tissue)) %>%
    dplyr::mutate(fraction_mutated = samples_mutated / samples_total) %>%
    dplyr::filter(!patient %in% removed_patients) %>%
    dplyr::mutate(group = groups[tissue]) %>%
    dplyr::select(gene, tissue, impact, fraction_mutated, group) %>%
    mutate(gene = factor(gene, order_genes, ordered = TRUE)) %>%
    mutate(tissue = factor(tissue, order_samples, ordered = TRUE)) %>%
    mutate(impact = factor(impact, names(impact_labels), impact_labels)) %>% 
    dplyr::arrange(-fraction_mutated * as.numeric(!is.na(impact))) %>% 
    dplyr::filter(!duplicated(paste0(tissue, gene)))
  
  
  for (ctissue in levels(data_plot$tissue)) {
    
    wh_tissue = data_plot$tissue == ctissue
    wh = wh_tissue & !is.na(data_plot$impact)
    wh_add = !data_plot$gene %in% data_plot$gene[wh]
    
    data_add = 
      data_plot[wh_add,] %>% 
      mutate(tissue = ctissue) %>% 
      mutate(impact = "Missense") %>% 
      mutate(fraction_mutated = NA) %>%
      mutate(group = groups[ctissue]) %>%
      unique()
    
    data_plot = rbind(data_plot, data_add)
  }
  
  data_plot = data_plot %>% 
    dplyr::filter(!is.na(impact))
  
  plot_heatmap =
    data_plot %>%
    ggplot(aes(x=tissue, y=gene, fill=fraction_mutated))
  
  for (ip in na.omit(unique(data_plot$impact))) {
    
    dat = data_plot %>% dplyr::filter(impact == ip)
    palette = c("Missense"="#335498", "Nonsense & Splice Site"="#B13430")[ip]
    
    plot_heatmap =
      plot_heatmap +
      ggnewscale::new_scale_fill() +
      geom_tile(
        data = dat,
        aes(fill = fraction_mutated),
        color = "gray0",
        size = 0.2
      ) +
      scale_fill_gradient( 
        na.value = "white",
        low="white",
        high = palette,
        #values = c(-0.2, 1.1),
        breaks = c(0.25, 0.5, 0.75, 1.0),
        limits = c(0, 1)
      ) + 
      #scale_fill_distiller(
      #  na.value = "white",
      #  direction = 1,
      #  palette = palette,
      #  values = c(-0.2, 1.1),
      #  breaks = c(0.25, 0.5, 0.75, 1.0),
      #  limits = c(0, 1)
      #) +
      labs(fill = paste0("Fraction mutated (", ip, ")")) + 
      geom_point(
        data = dplyr::filter(dat, fraction_mutated == 1),
        aes(color = "Clonal"),
        size = 0.8
      )
      
  }
  
  fill_guide = guide_legend(
    title.position = "top",
    title.hjust = 0.5,
    override.aes = list(shape = NA),
    breaks = c(0.25, 0.5, 0.75, 1.0)
  )
  
  plot_heatmap = 
    plot_heatmap + 
    guides(fill = fill_guide, fill_new = fill_guide, fill_new_new  = fill_guide)
  
  
  plot_heatmap = 
    plot_heatmap +
    scale_color_manual(values=c("Clonal"="darkorange")) + 
    theme_cowplot() + 
    scale_x_discrete(breaks=order_samples, labels=gsub("[(]", " (", gsub("[.][[:alnum:][:space:]-]*", "", order_samples))) +
    scale_y_discrete() + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
    theme(legend.position="bottom", legend.box="horizontal") +
    theme(strip.text = element_text(size=11), strip.background = element_blank()) + 
    facet_grid(.~group, labeller=labeller(.default=Hmisc::capitalize), scales="free", space="free") + 
    xlab("") + ylab("") + labs(color="") +  
    guides(color = guide_legend(title.position="top", title.hjust = 0.5)) +
    guides(fill = guide_legend(title.position="top", title.hjust = 0.5, override.aes = list(shape=NA))) +
    theme(legend.box = "vertical")
  

  return(plot_heatmap)
  
}

dnds_to_df = function(d) {
  
  if ("globaldnds_knownsites" %in% c(names(d), names(d[[1]]))) {
    # parse standard site dnds results
    
    if ("globaldnds_knownsites" %in% names(d)) {
      d = list("gene_set"=d)
    }
    
    res = 
      lapply(d, "[[", "globaldnds_knownsites") %>% 
      lapply(function(x) as.data.frame(t(x))) %>% 
      reshape2::melt(id.vars=c("obs","exp","dnds","cilow","cihigh")) %>% 
      dplyr::mutate(variant_group = factor(L1, names(d))) %>% 
      dplyr::mutate(type_name = factor("sitednds")) %>% 
      dplyr::mutate(mle = dnds)
    
    
  } else if ("globaldnds" %in% c(names(d), names(d[[1]]))) {
    # parse standard dnds results
    
    if ("globaldnds" %in% names(d)) {
      d = list("site_set"=d)
    }
    
    for (i in seq_along(d)) {
      if (!"names" %in% colnames(d[[i]]$globaldnds)) {
        d[[i]]$globaldnds$name = sub("_.*", "", rownames(d[[i]]$globaldnds))
      }
    }
    
    res = 
      lapply(d, "[[", "globaldnds") %>% 
      reshape2::melt(id.vars=c("name","mle","cilow","cihigh")) %>% 
      dplyr::mutate(variant_group = factor(L1, names(d))) %>% 
      dplyr::mutate(type_name = factor(name))
  }
  
  return(res)
}
