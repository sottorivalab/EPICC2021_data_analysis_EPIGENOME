load_peak_callsets = function(files) {
  files %>% 
    magrittr::set_names(gsub("-GRCh38-filtered_peaks[.]bed", "", basename(.))) %>%
    pbapply::pblapply(readr::read_tsv, col_type=readr::cols(), progress=FALSE) %>% 
    lapply(function(x) arrange(x, -fold_enrichment)) %>% 
    lapply(as, "GRanges") 
}

label_data = function(tr) {
  
  n = sapply(tr, NROW)
  stopifnot(n == max(n) | n == 0)
  n = max(n)
  if (n == 0) return(NA)
  
  label = rep("", n)
  
  for (test in names(tr)) {
    if (is.null(tr[[test]])) next()
    wh = which(tr[[test]]$sig)
    label[wh] = paste0(label[wh], ",", test)
  }
  
  label = gsub("^,", "", label)
  label[label == ""] = NA
  
  return(label)
}

exp_atac_signal_cn = function(purity, cn, ploidy=NULL, sig_n=1) {
  if (is.null(ploidy)) ploidy = mean(cn, na.rm=TRUE)
  site_signal = (sig_n * (1-purity) * 2 + sig_n * cn * purity)
  global_signal = (sig_n * (1-purity) * 2 + sig_n * ploidy * purity)
  site_signal / global_signal
}

recal_fit = function(model) {
  
  ls = model$samples$lib.size
  nf = model$samples$norm.factors
  coef = model$coefficients
  design = model$design
  
  fit = t(round(ls * nf * exp(design %*% t(coef))))
  dimnames(fit) = dimnames(model$fitted.values)
  model$fitted.values = fit
  
  model$deviance = 
    edgeR:::nbinomDeviance(
      model$counts, 
      fit, 
      model$dispersion
    )
  
  return(model)
}

edgeRGlmExp = function(model, coef, lfc=0) {
  
  get_stats = function(x, y, coef) {
    
    if (coef %in% colnames(x$coefficients)) {
      c1 = x$coefficients[, coef]
    } else {
      c1 = 0
    }
    
    if (coef %in% colnames(y$coefficients)) {
      c2 = y$coefficients[, coef]
    } else {
      c2 = 0
    }
    
    logFC = (c1 / log(2)) / (c2 / log(2))
    LR = y$deviance - x$deviance
    df.test = y$df.residual - x$df.residual
    LRT.pvalue = pchisq(LR, df = df.test, lower.tail = FALSE, log.p = FALSE)
    
    data.frame(
      logFC = logFC,
      logCPM = x$AveLogCPM,
      LR = LR,
      PValue = LRT.pvalue,
      row.names = rownames(x)
    )
  }
  
  # test the reduced model
  design = as.matrix(model$design)
  coef_n = match(coef, colnames(design))
  design0 = design[, -coef_n, drop = FALSE]
  design0[,"(Intercept)"] = as.numeric(colnames(model$counts) == "pan_patient-region_E")
  
  model0 =
    edgeR::glmFit(
      model$counts,
      design = design0,
      offset = model$offset,
      weights = model$weights,
      dispersion = model$dispersion,
      prior.count = model$prior.count,
      offsets = -5
    )
  

  # test alternative model
  model_lfc = model
  model_lfc$coefficients[,coef] = lfc * log(2)
  model_lfc = recal_fit(model_lfc)
  model_lfc$df.residual = model0$df.residual
  
  
  # gather stats:
  stats_inc_cn = get_stats(model, model_lfc, coef)
  stats_excl_cn = get_stats(model, model0, coef)
  colnames(stats_excl_cn) = paste0(colnames(stats_excl_cn), "_excl_cn")
  
  wh_clns = c("LR_excl_cn","PValue_excl_cn")
  tab = cbind(stats_inc_cn,  stats_excl_cn[,wh_clns])
  tab$exp_lfc = lfc
  
  return(tab)
}

plot_atac_cn_adjustment = function(peaks, coef, lfc=0, cn=2, title="") {
  
  color_scale = list(scale_color_viridis_d(), scale_fill_viridis_d())
  d_pos = data.frame(as(peaks, "GRanges"))
  d_pl = cbind(d_pos, data.frame(log_fc=lfc, coef=coef, cn=cn))
  
  mean_per_cn = 
    d_pl %>% 
    dplyr::group_by(
      cn = round(cn)
    ) %>% 
    dplyr::summarise(
      coef = mean(coef, na.rm=TRUE), 
      log_fc = mean(log_fc, na.rm=TRUE)
    )
  
  pl1 =
    ggplot(d_pl, aes(
      x = start,
      y = coef,
      color = factor(round(cn))
    )) +
    facet_grid( ~ seqnames, scales = "free_x", space = "free_x") +
    color_scale +
    ylim(quantile(d_pl$coef, 0.001, na.rm = TRUE), NA) +
    geom_point(alpha = 0.8) +
    geom_line(aes(y = log_fc), color = "black") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab("Position") +
    ylab("Observed coefficient") +
    labs(color = "Average CN") +
    theme(legend.position = "bottom") +
    theme(strip.text.x = element_text(angle = 90)) +
    theme(strip.background.x = element_blank())
  
  pl2 =
    ggplot(d_pl, aes(
      x = log_fc,
      y = coef,
      color = factor(round(cn))
    )) +
    geom_point(alpha = 0.8) +
    geom_point(
      data = mean_per_cn,
      shape = "x",
      color = "black",
      size = 4
    ) +
    geom_abline(intercept = 0, slope = 1) +
    ylim(quantile(d_pl$coef, 0.001, na.rm = TRUE), NA) +
    color_scale +
    xlab("Expected coefficient") +
    ylab("Observed coefficient") +
    labs(color = "Average CN")
  
  pl3 =
    ggplot(d_pl, aes(
      x = factor(round(cn)),
      y = coef,
      color = factor(round(cn))
    )) +
    color_scale +
    THmisc::pretty_violin() +
    ylab("Observed coefficient") +
    xlab("CN") +
    guides(color="none")
  
  cowplot::plot_grid(
    ggplot() + cowplot::theme_minimal_grid() + ggtitle(title) + theme(plot.title = element_text(hjust=0.5)), 
    pl1,
    cowplot::plot_grid(pl2 + ggtitle(""), pl3 + ggtitle(""), labels = c("B", "C")),
    nrow = 3,
    labels = c("","A", ""),
    rel_heights = c(0.1, 1, 0.9)
  )
}

test_results_to_analysis_set = function(results) {
  
  el_pairs = list(
    p = "PValue",
    p_no_cn = "PValue_excl_cn",
    max_cpm = "max_cpm",
    fc = "fc",
    cpm_baseline = 1
  )
  
  annot_group =
    strsplit(names(results), "-") %>% 
    do.call(what = rbind) %>% 
    data.frame() %>% 
    dplyr::mutate(X4 = paste0(X4, "-", X5)) %>%
    magrittr::set_colnames(c("set","patient","tss_type","sample"))
  
  stopifnot(length(unique(annot_group[,"set"])) == 1)
  col_ids = unique(paste0(annot_group[,"patient"], "-", annot_group[,"sample"]))
  row_ids = sort(unique(unlist(lapply(results, rownames))))
  
  null_matrix = matrix(
    NA,
    nrow = length(row_ids),
    ncol = length(col_ids),
    dimnames = list(row_ids, col_ids)
  )
  
  analysis_results = list(
    p = null_matrix,
    p_no_cn = null_matrix,
    max_cpm = null_matrix,
    fc = null_matrix,
    sig = null_matrix,
    cpm_baseline = null_matrix,
    ids = row_ids,
    set = magrittr::set_names(rep("", length(row_ids)), row_ids)
  )
    
  for (i in seq_len(NROW(annot_group))){
    
    # identify current subset
    c_results = results[[i]]
    idx_expr = strsplit(annot_group$sample[i], "_vs_", fixed = TRUE)[[1]]
    idx_cn = with(annot_group[i,], paste0(patient, "-", sample))
    idx_rn = rownames(c_results)
    c_results$fc = log2(c_results[,idx_expr[2]] / c_results[,idx_expr[1]])
    c_results$max_cpm = apply(c_results[,idx_expr], 1, max)
    
    # set elements in patient matrix
    for (j in seq_along(el_pairs)) {
      analysis_results[[names(el_pairs)[j]]][idx_rn,idx_cn] = 
        c_results[,el_pairs[[j]]]
    }
    
    # set set id
    tss_type = annot_group$tss_type[i]
    c_tss_type = analysis_results$set[idx_rn]
    stopifnot(c_tss_type == tss_type | c_tss_type == "")
    analysis_results$set[idx_rn] = tss_type
  }
  
  # peak locations
  analysis_results$peak_locations = 
    analysis_results$ids %>% 
    as("GRanges")
  
  return(analysis_results)
}

add_peak_annotation_to_analysis_set = function(analysis_results) {
  
  # add peak annotation 
  analysis_results$peak_annotation = 
    analysis_results$ids %>% 
    as("GRanges") %>% 
    ChIPseeker::annotatePeak(
      TxDb = txdb,
      annoDb = annodb,
      addFlankGeneInfo = TRUE,
      overlap = "TSS",
      flankDistance = 3000
    )
  
  return(analysis_results)
}

add_tss_overlap_to_analysis_set = function(analysis_results) {
  
  #add additional annotations of flanking gene
  overlaps_tss = 
    analysis_results$ids %>% 
    as("GRanges") %>%
    findOverlaps(.tss_pos - 1000 + 3000)
  
  # flanking gene ids by ids
  tg_name = .tss_pos$tg_name[subjectHits(overlaps_tss)]
  ids = analysis_results$ids[queryHits(overlaps_tss)]
  gids = lapply(split(tg_name, ids)[analysis_results$ids], unique)
  names(gids) = analysis_results$ids
  analysis_results$flanking_genes = gids
  
  # flanking gene symbols
  ensgid =
    tg_name %>% unique() %>% 
    AnnotationDbi::select(x=org.Hs.eg.db, keytype = "ENSEMBL", columns="SYMBOL")
  
  symbol_map = ensgid$SYMBOL
  names(symbol_map) = ensgid$ENSEMBL
  
  sids = lapply(gids, function(x) magrittr::set_names(symbol_map[x], NULL))
  analysis_results$flanking_symbol = sids
  
  return(analysis_results)
}

add_genhancer_overlap_to_analysis_set = function(analysis_results) {
  
  gh_loc = as(.genehancer_data$geneHancerRegElements[,1:3], "GRanges")
  gh_overlaps = findOverlaps(analysis_results$peak_locations, gh_loc)
  analysis_results$gh_overlaps = gh_overlaps
  
  # add annotations to result set
  for (el in c("elementType", "name")) {
    
    el_data = .genehancer_data$geneHancerRegElements[[el]]
    gh_annotation = rep(NA, length(analysis_results$ids))
    
    annot_per_query = 
      split(el_data[subjectHits(gh_overlaps)], queryHits(gh_overlaps)) %>%
      lapply(strsplit, "/") %>% lapply(unlist) %>% 
      lapply(unique) %>% lapply(sort) %>% sapply(paste, collapse=",")
    
    gh_annotation[as.numeric(names(annot_per_query))] = annot_per_query
    analysis_results[[paste0("gh_annotation_", el)]] =  gh_annotation
    
  }
  
  return(analysis_results)
}

get_final_reccurence_count = function(analysis_results, p_col="p", adj_method="none", min_fc=2, min_cpm=15) {
  
  final_reccurence_count = list()
  
  for (type in c("promoter","enhancer")) {
    
    # indexes 
    if (type == "promoter") {
      
      # is adjacent to gene tss?
      #wh_feature = analysis_results$set == "proximal"
      
      # is the adjacent gene expressed?
      n_flanking = sapply(analysis_results$flanking_genes, length)
      flanking_genes = unlist(analysis_results$flanking_genes)
      idx_flanking = rep_each_by(seq_along(n_flanking), n_flanking)
      idx_expressed = flanking_genes %in% rownames(.data_rnaseq_normalised)
      wh_feature = seq_along(n_flanking) %in% idx_flanking[idx_expressed]
      
      # preselection
      wh_entries = analysis_results$max_cpm > min_cpm & abs(analysis_results$fc) > log2(min_fc)
      wh_entries = wh_entries[wh_feature,]
      
      # entry ids
      gids = analysis_results$flanking_genes[which(wh_feature)]
      symb = analysis_results$flanking_symbol[which(wh_feature)]
      expr = lapply(gids, "%in%", rownames(.data_rnaseq_normalised))
      entry_ids = sapply(Map("[", symb, expr), paste0, collapse = ",")
      
    } else if (type == "enhancer") {
      n_flanking = sapply(analysis_results$flanking_genes, length)
      wh_feature = n_flanking == 0 & grepl("Enhancer", analysis_results$gh_annotation_elementType)
      
      # preselection
      wh_entries = analysis_results$max_cpm > min_cpm & abs(analysis_results$fc) > log2(min_fc)
      wh_entries = wh_entries[wh_feature,]
      
      # entry ids
      entry_ids = analysis_results$gh_annotation_name[wh_feature]
      
    } else {
      stop()
    }
    
    if (sum(wh_entries, na.rm=TRUE) == 0) next()
    
    # p value
    p_vals = analysis_results[[p_col]][rownames(wh_entries),]
    p_vals[!wh_entries] = NA
    
    #
    fc = analysis_results$fc[rownames(wh_entries),]
    
    stopifnot(all(colnames(fc) == colnames(p_vals)))
    
    # p.adjustment
    p_adj = apply(p_vals, 2, p.adjust, method=adj_method)
    #p_adj = p.adjust(c(unlist(p_vals)), "fdr")
    #dim(p_adj) = dim(p_vals)
    #dimnames(p_adj) = dimnames(p_vals)
    
    # significant matrix
    sig_matrix = p_adj < 0.01
    
    col_for_rec = !grepl("NvA", colnames(p_vals)) 
    
    rec_count = 
      lapply(c(1, -1), function(x) {
        apply(sig_matrix[,col_for_rec] & sign(fc[,col_for_rec]) == x, 1, sum, na.rm = TRUE)
      }) %>% do.call(what=cbind) %>% data.frame()
    
    colnames(rec_count) = c("gain", "loss")
    event_type = colnames(rec_count)[apply(rec_count, 1, which.max)]
    event_type[rec_count$gain == rec_count$loss] = "gain/loss"
    event_type[rec_count$gain == 0 & rec_count$loss == 0] = NA
    rec_count$recurrence = apply(rec_count, 1, max)
    
    summary =
      data.frame(
        peak = rownames(p_vals),
        id = entry_ids,
        rec_count,
        type = type,
        event_type = event_type
      )
    
    el_list = c(
      p = "p_vals",
      p_adj = "p_adj",
      sig_matrix = "sig_matrix",
      summary = "summary",
      fc = "fc"
    )
    
    for (i in seq_along(el_list)) {
      final_reccurence_count[[names(el_list)[i]]] =
        rbind(get(el_list[i]), final_reccurence_count[[names(el_list)[i]]])
    }
  }
  
  final_reccurence_count$flanking_genes = analysis_results$flanking_genes
  final_reccurence_count$flanking_symbol = analysis_results$flanking_symbol
  #final_reccurence_count$geneHancerInteractions = analysis_results$geneHancerInteractions
  
  return(final_reccurence_count)
}

rep_each_by = function(a, n) unlist(sapply(seq_along(a), function(i) rep(a[i], n[i])))

plot_test_results_regions = function(x, sig_cutoff=0.01, p_val="PValue", title="", labels=NULL) {
  
  x$ratio = x[,2] / x[,1]
  x$min_val = apply(x[,1:2], 1, min)
  x$sig = x[,p_val] < sig_cutoff
  x$peak = rownames(x) 
  
  cn = Hmisc::capitalize(colnames(x))
  ylab = paste0("CPM ", cn[2], "/", cn[1])
  
  gr = x %>% 
    rownames() %>% 
    as("GRanges") %>% 
    as.data.frame() %>% 
    cbind(x) %>% 
    dplyr::filter(min_val > 0) %>%
    dplyr::mutate(sig=ifelse(sig, "A", NA))

  plt = ggplot(gr, aes(x=start, y=ratio)) + 
    ggrastr::geom_point_rast(data=gr %>% dplyr::filter(is.na(sig)), alpha=0.8, aes(color=sig)) + 
    geom_point(data=gr %>% dplyr::filter(!is.na(sig)), alpha=0.8, aes(color=sig)) + 
    geom_line(aes(x=start, y=(exp_fc))) + 
    scale_y_log10() + 
    facet_grid(~seqnames, space = "free_x", scales="free_x") + 
    theme(axis.text.x = element_blank()) + 
    theme(strip.text.x = element_text(angle=90)) +
    theme(strip.background.x = element_blank()) +
    xlab("Position") + ylab(ylab) + 
    scale_color_brewer(palette = "Set1", direction = -1, na.value="gray50") + 
    ggtitle(title) + 
    guides(color="none")
  
  if (!is.null(labels)) {
    
    d_lab = gr %>% 
      dplyr::filter(peak %in% names(labels)) %>% 
      dplyr::mutate(label=labels[as.character(peak)])
    
    plt = plt + 
      ggrepel::geom_text_repel(
        data = d_lab %>% dplyr::filter(ratio < 1), 
        aes(label = label), 
        nudge_y = -0.5,
        min.segment.length = 0
      ) + 
      ggrepel::geom_text_repel(
        data = d_lab %>% dplyr::filter(ratio > 1), 
        aes(label = label), 
        nudge_y = 0.5,
        min.segment.length = 0
      )
    
  }

  return(plt)
}

plot_recurrence_heatmap = function(x, subclonal_status=NULL, rnaseq_res=NULL, n_per_group=20, alt_labelling_tissue=NULL, wh_peaks=NULL, drop_peak=c(), add_peak=c(), tissues_for_ordering="Cancer", suborder_cases=TRUE) {
  
  library(ggplot2)
  library(cowplot)
  library(tidyr)
  library(dplyr)
  theme_set(theme_cowplot())
  
  # significant fold change
  s_matrix = x$sig_matrix
  fc_sign = sign(x$fc[rownames(s_matrix),colnames(s_matrix)])
  fc_sign[is.na(s_matrix) | !s_matrix] = 0
  
  # correct direction of change
  n_losses = apply(fc_sign < 0, 1, sum)
  n_gains = apply(fc_sign > 0, 1, sum)
  change_type = case_when(n_losses < n_gains~1, n_losses > n_gains ~ -1, TRUE ~ 0)
  names(change_type) = names(n_losses)
  
  for (i in seq_len(NCOL(fc_sign))) { 
    fc_sign[which(fc_sign[,i] != change_type),i] = 0  
  }
  
  if (!is.null(rnaseq_res)) {
    min_p_rnaseq = with(rnaseq_res, tapply(p_use, peak.peak, min))
  } else {
    peaks_u = unique(rnaseq_res$peak.peak)
    min_p_rnaseq = magrittr::set_names(rep(NA, length(peaks_u)), peaks_u)
  }
  
  # 
  peaks_to_plot = 
    x$summary %>% 
    mutate(peak=as.character(peak)) %>% 
    mutate(event_type = change_type[peak]) %>% 
    mutate(event_type_chr = factor(event_type, c(-1,1), c("loss","gain"))) %>% 
    mutate(p_rna_correlation = min_p_rnaseq[peak]) %>% 
    mutate(used_before=peak %in% rownames(sc_cancer)) %>%
    arrange(-recurrence, -used_before, id) %>% 
    # exclude non consistent changes and keep top 5 gain/losses for enhancers/promotors
    filter(!is.na(event_type_chr))
  
  
  if (is.null(wh_peaks)) {
    peaks_to_plot_head = 
      peaks_to_plot %>%
      # split by event type and select most reccurent ones
      split(., list(.$type, .$event_type_chr)) %>% 
      lapply(head, n=n_per_group) %>% 
      #lapply(function(x) x[unique(c(1:20, which(x$recurrence >=10 | x$peak %in% rownames(sc_cancer)))),]) %>% #head, n=20) %>% 
      #lapply(function(d) d[which(d$p_rna_correlation < 0.01 | seq_len(NROW(d)) <= 10),]) %>% 
      do.call(what=rbind) %>%
      dplyr::filter(!peak %in% drop_peak)
    
    peaks_to_plot_filt = 
      peaks_to_plot %>% 
      dplyr::filter(peak %in% add_peak) %>%
      dplyr::filter(!peak %in% drop_peak)
    
    peaks_to_plot = 
      rbind(peaks_to_plot_head, peaks_to_plot_filt) %>% 
      unique() %>%
      arrange(-recurrence)  %>% 
      mutate(id=as.character(id)) %>% 
      mutate(id=factor(peak, peak[!duplicated(peak)], (id[!duplicated(peak)])))
    
  } else {
    peaks_to_plot = 
      peaks_to_plot %>% 
      dplyr::filter(peak %in% wh_peaks) %>% 
      mutate(id=factor(peak, peak[!duplicated(peak)], (id[!duplicated(peak)])))
  }
 
  #
  effect_group_labels = # final label for plot
    c("Gained Promoters"="Gained Promoters", 
      "Lossed Promoters"="Lost Promoters", 
      "Gained Enhancers"="Gained Enhancers", 
      "Lossed Enhancers"="Lost Enhancers")
  
  plot_data = # actual statistics for the peaks 
    fc_sign[peaks_to_plot$peak,] %>% 
    as.data.frame() %>%
    magrittr::set_colnames(colnames(fc_sign)) %>%
    mutate(peak=rownames(.)) %>% 
    reshape2::melt(id.vars="peak") %>% 
    magrittr::set_colnames(c("peak","patient","effect")) %>% 
    merge(peaks_to_plot, all.x=TRUE, by="peak") %>% 
    dplyr::filter(grepl("normal_vs_", patient)) %>% 
    mutate(effect_group=paste0(Hmisc::capitalize(as.character(event_type_chr)), "ed ", Hmisc::capitalize(type), "s")) %>% 
    mutate(effect_group=factor(effect_group, names(effect_group_labels), effect_group_labels)) %>%
    mutate(tissue=Hmisc::capitalize(gsub(".*normal_vs_", "", patient))) %>%
    mutate(tissue=factor(tissue, unique(c("Adenoma", "Cancer", "Pure", tissue)))) %>%
    mutate(patient=gsub("[-].*", "", patient)) %>% 
    mutate(effect = ifelse(effect == c("gain"=1, "loss"=-1)[as.character(event_type_chr)], 1, NA))
  

  if (!is.null(alt_labelling_tissue)) {
    plot_data = plot_data %>% 
      dplyr::mutate(tissue_label = factor(tissue, names(alt_labelling_tissue), alt_labelling_tissue)) %>% 
      dplyr::filter(!is.na(tissue_label))
  } else {
    plot_data$tissue_label = plot_data$tissue
  }

  # set order of peaks and patients 
  order_peaks = rev(peaks_to_plot$peak[order(peaks_to_plot$recurrence)])
  
  wh = plot_data$tissue_label %in% tissues_for_ordering
  effect_matrix = with(plot_data[wh,], tapply(effect, list(peak, patient), c))[order_peaks,]
  effect_matrix[is.na(effect_matrix)] = 0
  
  effect_list = split(effect_matrix, seq_len(NROW(effect_matrix)))
  ord_pats = rev(colnames(effect_matrix)[do.call(order, effect_list)])

  plot_data = 
    plot_data %>% 
    mutate(peak=factor(peak, rev(order_peaks), ordered = TRUE)) %>% 
    mutate(id=factor(id, rev(unique(rev(id[order(peak)]))), ordered = TRUE)) %>% 
    mutate(patient=factor(patient, ord_pats, ordered = TRUE))
  
  
  # RNA-seq tests to add to plot:
  plot_data_rna = 
    dplyr::filter(plot_data, p_rna_correlation < 0.01 & !is.na(effect))
  
  # create the plot
  plot_heatmap =
    plot_data %>%
    mutate(effect=ifelse(is.na(effect), 0, effect)) %>%
    mutate(x_pos = paste0(patient, "/", tissue)) %>% 
    ggplot(aes(x=x_pos, y=id, fill=event_type_chr, alpha=effect*0.9)) + 
    # geom layers:
    geom_tile(color="gray0", size=0.2) 
    
  if (NROW(plot_data_rna)) {
    plot_heatmap = plot_heatmap + 
      geom_point(data=plot_data_rna, aes(color="Association with RNA-seq"), size=0.8)
  }
  
  plot_heatmap = plot_heatmap +
    facet_grid(effect_group~tissue_label, scales="free", space="free") + 
    # scales and legends
    scale_color_manual(values=c("Association with RNA-seq"="gray10")) + 
    scale_alpha_identity() + 
    scale_fill_manual(na.value="white", values = c("#E62F39", "#4F8DC1"), breaks = c("gain","loss")) + 
    guides(color = guide_legend(title.position="top", title.hjust = 0.5)) + 
    guides(fill = "none") +
    # modify theme:
    theme_cowplot() + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
    theme(legend.position="bottom", legend.box="horizontal") +
    theme(strip.text = element_text(size=11), strip.background = element_blank()) + 
    # labels
    xlab("") + 
    ylab("") + 
    labs(color = "", fill = "Fraction of samples mutated") 
  
  
  # fix order of the axis labels:
  x_axis_breaks = paste0(plot_data$patient, "/", plot_data$tissue)
  wh = !duplicated(x_axis_breaks)
  x_axis_breaks = x_axis_breaks[wh]
  x_axis_labels = plot_data$patient[wh]
  x_axis_ord = order(plot_data$patient[wh])
  
  plot_heatmap$data$x_pos = 
    factor(
      plot_heatmap$data$x_pos, 
      as.character(x_axis_breaks[x_axis_ord])
    )
  
  plot_heatmap = plot_heatmap + 
    scale_x_discrete(
      breaks = as.character(x_axis_breaks), 
      labels = as.character(x_axis_labels)
    )
  
  
  # fix order of y axis labels:
  y_labels = 
    plot_data %>% 
    dplyr::select(id, peak) %>% unique() %>% 
    dplyr::arrange(id) %>% 
    dplyr::mutate(id = as.character(id)) %>% 
    dplyr::mutate(id_l = id)
  
  if (exists(".longer_enh_label")) {
    wh_enh = y_labels$id %in% names(.longer_enh_label)
    y_labels$id_l[wh_enh] = .longer_enh_label[as.character(y_labels$id)[wh_enh]]
  }
  
  plot_heatmap = plot_heatmap + 
    scale_y_discrete(
      breaks = y_labels$id, 
      labels = y_labels$id_l
    )
  
  if (!is.null(subclonal_status)) {
    
    stopifnot("subclonal" %in% colnames(subclonal_status))
    
    wh_missing_pat = !plot_heatmap$data$patient %in% subclonal_status$patient
    if (any(wh_missing_pat)) {
      wrn = paste0(
        "Missing subclonal data for patients: ", 
        paste0(unique(sort(plot_heatmap$data$patient[wh_missing_pat])), collapse=", ")
      )
      warning(wrn)
    }
    
    wh_missing_peak= !plot_heatmap$data$peak %in% subclonal_status$peak
    if (any(wh_missing_peak)) {
      wrn = paste0(
        "Missing subclonal data for peaks: ", 
        paste0(unique(sort(plot_heatmap$data$peak[wh_missing_peak])), collapse=", ")
      )
      warning(wrn)
    }
    
    wh_cols = colnames(subclonal_status)[colnames(subclonal_status) %in% colnames(plot_heatmap$data)]
    is_subclonal = rep(TRUE, NROW(plot_heatmap$data))
    for (i in seq_along(is_subclonal)) {
      ol = subclonal_status$subclonal
      for (j in seq_along(wh_cols)) {
        ol = ol & plot_heatmap$data[i,wh_cols[j]] == subclonal_status[,wh_cols[j]]
        stopifnot(!any(is.na(ol)))
      }
      is_subclonal[i] = any(ol)
    }
    
    plot_heatmap$data$effect[is_subclonal] = 0.4
      #plot_heatmap$data$effect[is_subclonal] * 0.4
    
  }

  return(plot_heatmap)
}


compare_recurrence_counts = function(x, y, wh_ids=NULL, plot_line=TRUE, plot_relative=FALSE) {
  
  id1 = x$summary$peak
  id2 = y$summary$peak
  ids = intersect(id1, id2)
  if (!is.null(wh_ids)) ids = intersect(wh_ids, ids)

  mt1 = match(ids, id1)
  mt2 = match(ids, id2)
  
  r1 = x$summary$recurrence[mt1]
  r2 = y$summary$recurrence[mt2]
  csqt = chisq.test(r1, r2, simulate.p.value = TRUE, B = 1000)
  
  if (plot_relative) {
    r1 = r1 / NCOL(x$p)
    r2 = r2 / NCOL(y$p)
    b1 = seq_len(NCOL(x$p)) /  NCOL(x$p)
    b2 = seq_len(NCOL(y$p)) /  NCOL(y$p)
  } else {
    r1_ = r1
    r2_ = r2
    b1 = seq_len(NCOL(x$p))
    b2 = seq_len(NCOL(y$p))
  }
  
  pt = x$summary$type[mt1]
  et = x$summary$event_type[mt1]
  
  pl = table(with_cn_adj=r2, without_cn_adj=r1, peak_type=pt) %>% 
    reshape2::melt() %>% 
    dplyr::filter(value>0) %>% 
    ggplot(aes(x=without_cn_adj, y=with_cn_adj, fill=log10(value))) + 
    facet_wrap(peak_type~., labeller = as_labeller(Hmisc::capitalize)) +
    geom_tile() + 
    scale_fill_viridis_c(na.value = "white") + 
    xlab("Recurrence without CN adjustment") + 
    ylab("Recurrence with CN adjustment") + 
    scale_x_continuous(minor_breaks = b1) + 
    scale_y_continuous(minor_breaks = b2) + 
    labs(caption=paste0("X2 = ", round(csqt$statistic, 4), ", p = ", signif(csqt$p.value, 3)))
    
  
  
  if (plot_line) {
    pl = pl + geom_abline(intercept = 0, slope=1)
  } 
  
  
  return(pl)
}



plot_loss_gain_numbers = function(x, label_n=4, nudge_x=2, nudge_y=2) {
  
  plt =
    x$summary %>% 
    dplyr::group_by(gain, loss, type) %>% 
    dplyr::summarise(n=length(type)) %>%
    mutate(type=Hmisc::capitalize(paste0(type, "s"))) %>% 
    #dplyr::mutate(n_lab=ifelse(nchar(as.character(n)) > 2, NA, n)) %>% 
    dplyr::mutate(n_lab = as.character(n)) %>%
    ggplot(aes(x=loss, y=gain, fill=n)) + 
    geom_tile() + 
    geom_text(aes(label=n_lab, color=ifelse(n<500, "gray10", "gray95")), size=2.5) + 
    facet_wrap(Hmisc::capitalize(type)~., scales='free') +
    scale_fill_viridis_c(direction = -1, trans="log10") + 
    theme(strip.background = element_blank()) + 
    xlab("Lost [# patients]") + 
    ylab("Gained [# patients]") + 
    labs(fill="N peaks") + 
    scale_color_identity()
  
  
  if (label_n) {
    
    add_label = function(x) {
      d = x %>% 
        dplyr::select(type, loss, gain) %>%
        unique() %>%
        dplyr::mutate(label=paste0(x$label, collapse = "\n"))
    }
    
    annot = 
      x$summary %>% 
      arrange(-recurrence) %>% 
      dplyr::filter(event_type %in% c("gain","loss")) %>%
      split(., paste0(.$event_type, .$type)) %>% 
      lapply(head, n=label_n) %>%
      do.call(what=rbind) %>% 
      mutate(loss=as.numeric(as.character(loss))) %>% 
      mutate(gain=as.numeric(as.character(gain))) %>%
      mutate(type=Hmisc::capitalize(paste0(type, "s"))) %>% 
      mutate(label = ifelse(type == "Promoters", id, .longer_enh_label[id])) %>% 
      split(., list(.$loss, .$gain, .$type)) %>% 
      lapply(add_label) %>% 
      do.call(what=rbind)
    
    
    plt = 
      plt +
      ggrepel::geom_text_repel(
        data = annot,
        aes(label = label, x = loss, y = gain + 0.5),
        lineheight = 0.75,
        inherit.aes = FALSE,
        force = 5,
        segment.alpha = 0.75,
        min.segment.length = 0,
        nudge_x = nudge_x,
        nudge_y = nudge_y,
        size = 3
      )
    
  }
  
  return(plt)
  
}


test_rnaseq_association_peaks = function(peaks, atac_dataset, dds_, n_bs=10000, test_data=TRUE, verbose=TRUE) {
  
  require(progress)
  #require(DESeq2)
  
  results = NULL
  
  if (verbose) pb = progress_bar$new("[:bar] :percent eta: :eta", NROW(peaks))
  
  for (i in seq_len(NROW(peaks))) {
    
    if (verbose) pb$tick()
    
    # retrieve cpm values for current peak
    peak = peaks[i,]
    
    rna_per_group = tryCatch(
      get_rna_data_for_peak(dds_, peak$peak, atac_dataset),
      error = function(e) {
        print(e)
        return(NULL)
      }
    )
    
    if (is.null(rna_per_group)) next()
    
    for (gene in names(rna_per_group)) {
      
      # get ids for current gene
      gene_symbol = gsub(".*[|]", "", gene)
      ensg_id = gsub("[|].*", "", gene)
      
      # otherwise prepare tests:
      d_per_gr = rna_per_group[[gene]]
      contrast = lapply(d_per_gr, names)
      contrast_deseq = lapply(contrast, function(x) paste0("Patient", gsub("[.]cancer", "", x)))
      exp_sig =ifelse(peak$type == "two.sided", "greaterAbs", ifelse(peak$loss > peak$gain, "less", "greater"))
      exp_sig_deseq = ifelse(peak$type == "enhancer", "greaterAbs", ifelse(peak$loss > peak$gain, "less", "greater"))
      
      # do tests:
      if (test_data) {
        p_wilcox = tryCatch(wilcox.test(d_per_gr$sig, d_per_gr$not_sig, exp_sig)$p.value, error=function(e) return(NA))
        p_wilcox_twosided =  tryCatch(wilcox.test(d_per_gr$sig, d_per_gr$not_sig)$p.value, error=function(e) return(NA))
        p_deseq = tryCatch(results(dds_[ensg_id,], contrast_deseq, altHypothesis=exp_sig_deseq)$pvalue, error=function(e) return(NA))
        p_deseq_force = tryCatch(results(dds_[ensg_id,], contrast_deseq, altHypothesis=exp_sig_deseq, cooksCutoff=Inf)$pvalue, error=function(e) return(NA))
        p_deseq_twosided = tryCatch(results(dds_[ensg_id,], contrast_deseq)$pvalue, error=function(e) return(NA))
        p_deseq_twosided_force = tryCatch(results(dds_[ensg_id,], contrast_deseq, cooksCutoff=Inf)$pvalue, error=function(e) return(NA))
        
        if (n_bs) {
          len = sapply(d_per_gr, length)
          gr = rep_each_by(names(d_per_gr), len)
          d = unlist(d_per_gr)
          d_bs = replicate(n_bs, split(d, sample(gr, replace=FALSE)), simplify=FALSE)
          d_bs = lapply(d_bs, "[", c("not_sig","sig"))
          mean_bs = lapply(d_bs, sapply, mean, na.rm=TRUE)
          diff_bs = sapply(mean_bs, diff)
          diff = diff(sapply(d_per_gr, mean, na.rm=TRUE)[c("not_sig","sig")])
          
          if (exp_sig_deseq == "greater") {
            p_bs = 1 - mean(c(0, diff > diff_bs))
          } else if (exp_sig_deseq == "less")  {
            p_bs = 1 - mean(c(0, diff > diff_bs))
          } else if (exp_sig_deseq == "greaterAbs")  {
            p_bs = 1 - mean(c(0, abs(diff) > abs(diff_bs)))
          } else  {
            stop()
          }
        } else {
          p_bs = NA
        }
        
      } else {
        p_wilcox = p_wilcox_twosided = p_deseq = p_deseq_force = p_deseq_twosided = p_deseq_twosided_force = p_bs =  NA
      }
      
      mean_per_group = sapply(d_per_gr, mean)
      fc_cpm = mean_per_group[2]/mean_per_group[1]
      
      # append test results:
      dtest = data.frame(
        peak = peak,
        gene = gene,
        symbol = gene_symbol,
        p_wilcox = p_wilcox,
        p_wilcox_twosided = p_wilcox_twosided,
        p_deseq = p_deseq,
        p_deseq_force = p_deseq_force,
        p_deseq_twosided = p_deseq_twosided,
        p_deseq_twosided_force = p_deseq_twosided_force,
        p_bs = p_bs,
        fc_cpm = fc_cpm,
        n_sig = length(d_per_gr$sig),
        n_not_sig = length(d_per_gr$not_sig),
        hypothesis = exp_sig
      )
      
      results = rbind(results, dtest)
    }
    
  }
  
  if (is.null(results)) 
    return(results)
  
  results = 
    results %>% 
    mutate(p_adj = p.adjust(p_deseq_force, "fdr"))
  
  return(results)
  
}


get_rna_data_for_peak = function(dds_, peak, atac_dataset, add_normal=FALSE) {
  
  require(progress)
  require(DESeq2)
  
  
  # effect type 
  s_matrix = atac_dataset$sig_matrix
  fc_sign = sign(atac_dataset$fc[rownames(s_matrix),colnames(s_matrix)])
  fc_sign[is.na(s_matrix) | !s_matrix] = 0
  n_losses = apply(fc_sign < 0, 1, sum)
  n_gains = apply(fc_sign > 0, 1, sum)
  
  peak = atac_dataset$summary[peak,]
  
  if (peak$type == "promoter") {
    genes = atac_dataset$flanking_genes[[peak$peak]]
    gene_symbols = atac_dataset$flanking_symbol[[peak$peak]]
  } else if (peak$type == "enhancer") {
    library(org.Hs.eg.db)
    wh = .genehancer_data$geneHancerInteractions$geneHancerIdentifier == as.character(peak$id)
    gh_id = gsub("/.*", "", .genehancer_data$geneHancerInteractions$name[wh])
    annot = tryCatch(suppressMessages(AnnotationDbi::select(gh_id, x=org.Hs.eg.db, keytype="SYMBOL", columns="ENSEMBL")), error=function(e) return(NULL))
    if (is.null(annot)) return(NULL)
    gene_symbols = annot$SYMBOL
    genes = annot$ENSEMBL
  } else {
    stop()
  }
  
  
  results = list()
  
  for (gene in genes) {
    
    gene_symbol = gene_symbols[which(genes == gene)]
    
    if (!gene %in% rownames(dds_)) {
      warning(paste0("Missing RNA-seq data for gene ", gene_symbol))
      next()
    }
    
    # atac values
    wh_pats_sig_atacseq = names(which(s_matrix[peak$peak,]))
    wh_not_pats_sig_atacseq = colnames(s_matrix)[!colnames(s_matrix) %in% wh_pats_sig_atacseq]
    
    # tests
    wh_pats_sig_atacseq = wh_pats_sig_atacseq[!grepl("adenoma", wh_pats_sig_atacseq)]
    wh_not_pats_sig_atacseq = wh_not_pats_sig_atacseq[!grepl("adenoma", wh_not_pats_sig_atacseq)]
    contrast = list(sig=wh_pats_sig_atacseq, not_sig=wh_not_pats_sig_atacseq)
    contrast = lapply(contrast, function(x) gsub("[-.].*", "", x))
    
    fpm_values = tapply(assay(normTransform(dds_))[gene,], dds_$Patient, mean)
    fpm_per_group = lapply(contrast, function(x) unlist(fpm_values[x[x %in% names(fpm_values)]]))
    
    if (add_normal) {
      wh_norm = dds_$Patient == "Normal"
      pat_id = gsub("_.*", "", colnames(dds_))[wh_norm]
      fpms =  assay(normTransform(dds_))[gene,wh_norm]
      fpm_per_group[["normal"]] = tapply(fpms, pat_id, mean)
    }
    
    results[[paste0(gene, "|", gene_symbol)]] = fpm_per_group
    
  }
  
  return(results)
}

get_purity_groups = function(d) {
  
  wh_sample = grepl("EPICC", colnames(d))
  annot = THmisc::annotation_from_barcode(colnames(d)[wh_sample], T)
  annot$purity_per_sample = atacseq_purity[annot$sample_barcode,"estimated_purity"]
  
  cut1 = # first cut, pure and impure
    annot$purity_per_sample %>% 
    cut(c(0,0.4,1), include.lowest = TRUE) %>% 
    factor(labels = c("impure","pure")) %>% 
    as.character() %>% 
    magrittr::set_names(colnames(d)[wh_sample])
  
  cut2 = # second cut, 4 groups in 0.25 intervals
    annot$purity_per_sample %>% 
    cut(c(0,0.2,0.4,0.6,1), include.lowest = TRUE) %>% 
    factor(labels = c("00","20","40","60")) %>% 
    as.character() %>% 
    magrittr::set_names(colnames(d)[wh_sample])
  
  wh_na = annot$tissue_type != "cancer" #is.na(purity_per_sample) |
  cuts = c(cut1[!wh_na], cut2[!wh_na])
  split_label = paste0(rep(annot$patient[!wh_na], 2), "-", cuts)
  purity_groups = split(rep(names(cut1)[!wh_na], 2), split_label)
  
  return(purity_groups)
}

which.min_all = function(x) {
  which(x == min(x, na.rm=TRUE))
}

get_closest_non_na = function(x) {
  wh_non_na = which(!is.na(x))
  wh_na = which(is.na(x))
  if (length(wh_na)) {
    .get_d= function(i) { mean(x[wh_non_na[which.min_all(abs(wh_non_na - i))]], na.rm = TRUE) }
    x[wh_na] = sapply(wh_na, .get_d)
  }
  return(x)
}

get_tissue_groups = function(d, tissue="cancer") {
  wh_sample = grepl("EPICC", colnames(d))
  annot = THmisc::annotation_from_barcode(colnames(d)[wh_sample], T)
  wh_adeno = annot$tissue_type == tissue
  ids = gsub("[.]", "-", annot$tumour_id[wh_adeno])
  names(ids) = NULL
  groups = split(colnames(d)[wh_sample][wh_adeno],  ids)
  return(groups)
}

get_normal_groups = function(d) {
  ids_ = grep("E[1-9]_G[0-9]+_C1", colnames(d), value = TRUE)
  ids = split(ids_, ids_)
  
  case = sapply(strsplit(names(ids), "-"), "[[", 1)
  gland = substr(sapply(strsplit(names(ids), "-"), "[[", 2), 16, 16)
  names(ids) = paste0(case, "-normal_gland", gland)
  
  ids2 = split(ids_, sapply(strsplit(ids_, "-"), "[[", 1))
  names(ids2) = paste0(names(ids2), "-normal_glands")
  
  return(c(ids, ids2))
}

get_data_for_groups = function(d, g) {
  matrix_purity_groups = 
    do.call(cbind, lapply(g, function(gr) {
      if (length(gr) == 0) return(NULL)
      apply(d[,unique(gr),drop=FALSE], 1, sum)
    }))
}

get_correction_factor_for_groups = function(d, g, sample_data) {
  
  if (is.list(g)) 
    return(lapply(g, get_correction_factor_for_groups, d=d, sample_data=sample_data))

  #
  smps_annot = THmisc::annotation_from_barcode(gsub("C[0-9]+-", "", g), TRUE)
  pur_per_smp = sample_data$purity[smps_annot$sample_barcode,"estimated_purity"]
  w_per_smp = apply(matrix[,g, drop=FALSE], 2, sum) %>% (function(x) x/sum(x))
  smps_region = with(smps_annot, paste0(patient, region))
  
  if (all(smps_annot$tissue_type == "normal")) {
    return(data.frame(cn=rep(2, NROW(d)), fc=rep(1, NROW(d)), row.names = rownames(d)))
  }
  
  # get cn data
  mt = match(smps_region, colnames(sample_data$cn_per_region))
  cn_per_smp = sample_data$cn_per_region[rownames(d), mt, drop=FALSE]
  mt_ = match(smps_annot$patient[is.na(mt)], colnames(sample_data$cn_per_cancer))
  cn_per_smp[,is.na(mt)] = sample_data$cn_per_cancer[rownames(d),mt_]
  for (i in seq_len(NCOL(cn_per_smp))) 
    cn_per_smp[,i] = get_closest_non_na(cn_per_smp[,i])
  
  # calc expected ratio
  fc_ratio = cn_per_smp; fc_ratio[TRUE] = 1
  for (i in seq_len(NCOL(fc_ratio))) {
    fc_ratio[,i] = 
      exp_atac_signal_cn(
        purity = pur_per_smp[i], 
        cn = cn_per_smp[,i]
      )     
  }
  
  # get weighted means
  cn_w = apply(cn_per_smp, 1, weighted.mean, w_per_smp, na.rm=TRUE)
  fc_w = apply(fc_ratio, 1, weighted.mean, w_per_smp, na.rm=TRUE)
  
  return(data.frame(cn=cn_w, fc=fc_w, row.names = rownames(d)))
}


get_dgelist_objects = function(data, fc_data) {
  
  ids_missing_fc_data = apply(is.na(fc_data), 2, all) %>% which %>% names
  
  # normalise and estimate dispersion:
  wh = (grepl("C5.*-cancer", colnames(data)) | colnames(data) == "all_normals") & !colnames(data) %in% ids_missing_fc_data
  dge_list_megabulks = edgeR::DGEList(data[,wh], group=rep(1, ncol(data[,wh])))
  dge_list_megabulks$offset = log(t(dge_list_megabulks$samples$lib.size * t(fc_data[,colnames(dge_list_megabulks)])))
  dge_list_megabulks %<>% edgeR::calcNormFactors(method="TMMwsp") %>% edgeR::estimateDisp() 
  dge_list_megabulks$samples$group = colnames(data[,wh])
  
  # get cpm values:
  wh = !colnames(data) %in% ids_missing_fc_data
  dge_list = edgeR::DGEList(data[,wh], group=rep(1, ncol(data[,wh])))
  dge_list$samples$group = rownames(dge_list$samples)
  
  dge_list_cn_corrected = dge_list
  dge_list_cn_corrected$offset = log(t(dge_list$samples$lib.size * t(fc_data[,colnames(dge_list)])))
  dge_list_cn_corrected$samples$group = rownames(dge_list$samples)
  
  dge_list$offset = NULL
  dge_list %<>% edgeR::calcNormFactors(method="TMMwsp") 
  
  # copy over dispersion estimates 
  for (cn in names(dge_list_megabulks)) {
    if (cn %in% c("counts","samples","offset")) next()
    dge_list[[cn]] = dge_list_megabulks[[cn]]
    dge_list_cn_corrected[[cn]] = dge_list_megabulks[[cn]]
  }
  
  return(list(cn_adj=dge_list_cn_corrected, unadj=dge_list))
}

get_input_data = function(counts, sample_data) {
  
  # construct a purity, tumour, etc based groups
  data_groups = c(
    get_purity_groups(counts),
    get_tissue_groups(counts, "cancer"),
    get_tissue_groups(counts, "adenoma"),
    get_normal_groups(counts),
    list("all_normals"= unlist(get_normal_groups(counts)))
  ) %>% lapply(magrittr::set_names, NULL)
  
  wh_keep = (grepl("C5", names(data_groups))  |
               names(data_groups) == "all_normals") &
    !grepl("-([0246]0)|(NA)|impure$", names(data_groups))
  data_groups = data_groups[wh_keep]
  
  
  # get expected fc ratios for each sample group:
  correction_factor_data = 
    get_correction_factor_for_groups(
      counts,
      data_groups,
      sample_data
    )
  
  fc_matrix = 
    do.call(
      what=cbind, 
      lapply(correction_factor_data, "[[", "fc")
    ) %>% magrittr::set_rownames(rownames(counts))
  
  cn_matrix = 
    do.call(
      what=cbind, 
      lapply(correction_factor_data, "[[", "cn")
    ) %>% magrittr::set_rownames(rownames(counts))
  
  
  # merge data of each group:
  matrix_megabulks = 
    get_data_for_groups(
      counts, 
      data_groups
    )
  
  
  return(list(
    fc = fc_matrix,
    cn = cn_matrix,
    counts = matrix_megabulks,
    groups = data_groups
  ))
}

do_tests = function(dge_data, input_data, group, dis_source="common") {
  
  # do statistical tests:
  stopifnot(dis_source == "common")
  disp = dge_data$unadj$common.dispersion
  
  # unadjusted
  gr_f = relevel(factor(dge_data$unadj$samples$group), "all_normals")
  design = model.matrix(~gr_f)
  colnames(design) = gsub("gr_f", "", colnames(design))
  glmqlfit = edgeR::glmFit(dge_data$unadj, design, dispersion=disp)
  
  # CN adjusted
  gr_f = relevel(factor(dge_data$cn_adj$samples$group), "all_normals")
  design = model.matrix(~gr_f)
  colnames(design) = gsub("gr_f", "", colnames(design))
  glmqlfit_adj = edgeR::glmFit(dge_data$cn_adj, design, dispersion=disp)
  
  # LRT
  test_res = glmLRT(glmqlfit, group)$table
  #test_res = edgeR::exactTest(dge_data$unadj, c("all_normals", group), dispersion = dis_source)$table
  test_res$coef = coefficients(glmqlfit)[,group]
  colnames(test_res) = paste0(colnames(test_res), "_excl_cn")
   
  test_res_adj = glmLRT(glmqlfit_adj, group)$table
  #test_res_adj = edgeR::exactTest(dge_data$cn_adj, c("all_normals", group), dispersion = dis_source)$table
  test_res_adj$coef = exp(coefficients(glmqlfit_adj)[,group])
  test_res_adj$sig = test_res_adj$PValue < 0.01
  
  # cpm data
  pair = c("all_normals", group)
  cpms = edgeR::cpm(dge_data$unadj[, pair], normalized.lib.sizes = TRUE, prior.count = 0.5) %>% data.frame()
  colnames(cpms) = pair
  
  # merge data
  cn_w = input_data$cn[rownames(glmqlfit),group]
  fc_w = input_data$fc[rownames(glmqlfit),group]
  tr = cbind(cpms, cn=cn_w, exp_fc=fc_w, test_res_adj, test_res)
  colnames(tr)[seq_along(pair)] = pair
  
  return(tr)
}

plot_scatter = function(d, a, b, annot=NULL) {
  
  if (!all(c(a, b) %in% colnames(d))) {
    blank_plot = ggplot() + theme_void()
    return(blank_plot)
  }
  
  if (!is.null(annot)) {
    d$label = as.character(annot[rownames(d),"id"])
    d$label[is.na(d$highlight)] = NA
  }
  
  p =
    d %>% 
    (function(x) x[order(x$highlight, na.last=FALSE),]) %>% 
    ggplot(aes(x=get(a), y=get(b), color=highlight)) + 
    ggrastr::geom_point_rast(alpha=0.9) + 
    xlab(tools::toTitleCase(paste0("CPM ", a))) + 
    ylab(tools::toTitleCase(paste0("CPM ", b))) +
    scale_color_brewer(palette="Set1", na.value="grey75") + 
    geom_abline(slope=1, intercept=0, linetype=2, color="grey30") + 
    geom_abline(slope=2, intercept=0, linetype=3, color="grey30") + 
    geom_abline(slope=4, intercept=0, linetype=3, color="grey30") + 
    geom_abline(slope=1/2, intercept=0, linetype=3, color="grey30") + 
    geom_abline(slope=1/4, intercept=0, linetype=3, color="grey30") + 
    labs(color="")
  
  if (!is.null(annot)) {
    p = p + ggrepel::geom_text_repel(
      aes(label = label),
      color = "black",
      min.segment.length = 0,
      segment.alpha = 0.6,
      size = 3,
      force = 1, 
      max.overlaps = 100
    )
  }
  
  return(p)
} 


get_driver_peaks = function(summary_data, genes = THmisc::gene_lists$`IntOGen-DriverGenes_COREAD`, min_rec = 1) {
  
  ids = 
    with(summary_data$summary %>% dplyr::filter(recurrence >= min_rec), {
      id = as.character(id)
      ifelse(id %in% names(.longer_enh_label), .longer_enh_label[id], id) %>% 
        magrittr::set_names(peak) %>% 
        gsub(pattern = " .*", replacement = "") %>% strsplit(",")
    })
  
  id_to_driver = lapply(ids, function(x) x[x %in% genes])
  id_is_driver = sapply(id_to_driver, length) > 0
  id_to_driver = id_to_driver[id_is_driver]
  len = sapply(id_to_driver, length)
  
  magrittr::set_names(rep_each_by(names(id_to_driver), len), unlist(id_to_driver))
}
