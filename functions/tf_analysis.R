read_macs2_peak_file = function(x) {
  stopifnot(file.exists(as.character(x)))
  cn = c("chr","start","end","length","abs_summit","pileup","neg_log10_p","fold_enrichment","neg_log10_q","name")
  readr::read_tsv(x, comment = "#", col_names = cn, skip = 25, col_types = readr::cols(), progress = FALSE)
}

read_filtered_peak_file = function(x) {
  stopifnot(file.exists(as.character(x)))
  cn = c("chr","start","end","length","abs_summit","pileup","neg_log10_p","fold_enrichment","neg_log10_q","name","rank","quantile")
  readr::read_tsv(x, comment = "#", col_names = cn, skip = 1, col_types = readr::cols(), progress = FALSE)
}


read_bed_file = function(f) {
  require(GenomicRanges)
  
  stopifnot(file.exists(f))
  
  fl = readr::read_tsv(f, col_names=FALSE, col_types = readr::cols(), n_max=2)
  if (!suppressWarnings(is.na(as.numeric(fl[1,2])))) { # asssume no header
    if (NROW(fl) == 0) return(NULL)
    bed_data = readr::read_tsv(f, col_names = FALSE, col_types = "cii")
  } else { # assume header
    if (NROW(fl) == 1) return(NULL)
    stopifnot(!suppressWarnings(is.na(as.numeric(fl[2,2])))) 
    bed_data = readr::read_tsv(f, col_names = TRUE, col_types = "cii")
  }
  colnames(bed_data) = c("chr","start","end")
  bed_data = as(bed_data, "GRanges")
  return(bed_data)
}


shorten_tf_name = function(x) {
  bn = basename(x)
  gr = gsub("-.*", "", bn)
  
  ngrn = c( # new group names
    "unclustered" = "uc",
    "tss_close" = "cTSS",
    "tss_distal" = "dTSS",
    "tss_proximal" = "pTSS",
    "not_overlapping_peak" = "nPEAK",
    "overlapping_peak" = "oPEAK"
  )
  
  for (i in seq_along(ngrn)) {
    gr = gsub(names(ngrn)[i], ngrn[i], gr)
  }
  
  tf = gsub("_.*", "", gsub(".*LINE[[:digit:]]+_", "", bn))
  
  paste0(tf, "_", gr)
}


read_shifted_count_file = function(x) {
  stopifnot(file.exists(x))
  d = readr::read_tsv(x, col_names = FALSE, col_types = "ciic", progress = FALSE)
  colnames(d) = c("chr","start","end","strand")
  as(d, "GRanges")
}


read_unshifted_count_file = function(x) {
  
  stopifnot(file.exists(x))
  
  # load file
  d = readr::read_tsv(x, col_names = FALSE, col_types = "ciicii--cc", progress = FALSE)
  colnames(d) = c("chr1","start1","end1","chr2","start2","end2","strand1","strand2")
  
  # minimal checks, might have to be changed
  stopifnot(all(d$chr1 == d$chr2))
  stopifnot(all(d$strand1 != d$strand2))
  
  # shift and return as granges with metadata
  chr = d$chr1
  start = ifelse(d$strand1 == "+", d$start1, d$start2) + 5
  end = ifelse(d$strand2 == "-", d$end2, d$end1) - 4
  
  pos = c(start, end)
  #strand = rep(c("+","-"), each=length(start))
  isize = rep(end - start, 2)
  gr = GRanges(rep(chr, 2), IRanges(pos, pos), insert_size=isize)
  
  sort(gr)
}


findOverlapsTF = function(reads, sites, include_offset=TRUE) {
  
  stopifnot("GRanges" %in% class(d))
  stopifnot("GRanges" %in% class(sites))
  stopifnot(all(width(reads) == 1))
  
  # find overlaps
  ol = findOverlaps(reads, sites, type="within", select="all", ignore.strand=TRUE)
  site_center = start(sites-floor(width(sites)/2))
  
  # summary of overlaps in data.frame
  if (include_offset) {
    offset = site_center[subjectHits(ol)] - start(reads[queryHits(ol)])
  } else {
    offset = NA
  }
  strand = as.character(strand(reads[queryHits(ol)]))
  insert_size = elementMetadata(reads[queryHits(ol)])$insert_size
  insert_size = cut(insert_size, seq(0, 10000, by=5), ordered_result = TRUE)
  tf = sites$tf[subjectHits(ol)]
  
  # count 
  df = data.frame(offset, insert_size, tf)
  df_c = plyr::count(df, vars=colnames(df))
  
  return(df_c)
}


plot2d_tf_grid = function() {
  
  plot_ctcf_2d_grid = 
    data_collection_ctcf %>% 
    ggplot(aes(x=offset, y=insert_size)) + 
    stat_bin_2d(binwidth=c(10, 5), geom="raster", drop=FALSE, size=0, linetype=0) + 
    ylim(0, 500) + 
    scale_fill_viridis_c(na.value = "#440154FF") + 
    xlab("Distance from motif center") + 
    ylab("Insert size") + 
    ggtitle("CTCF motifs") + 
    labs(fill="N insertions")
  
}

plot2d_tf_grid = function(d) {
  
  gr_to_numeric = function(x) {
    
    if (is.factor(x)) {
      lis = levels(x)
      nvars = strsplit(gsub("[(]", "", gsub("[]]", "", lis)), ",") %>% lapply(as.numeric)
      nval = sapply(nvars, max)
      names(nval) = lis
      x_return = nval[as.character(x)]
    } else {
      x_return = x
    }
    
    return(x_return)
  }
  
  
  p = d %>% 
    reshape::melt() %>% 
    magrittr::set_colnames(c("insert_size","offset","freq_plot")) %>% 
    mutate(offset=as.numeric(as.character(offset))) %>% 
    mutate(insert_size=gr_to_numeric(insert_size)) %>%
    ggplot(aes(x=offset, y=insert_size, z=freq_plot)) + 
    stat_summary_2d(binwidth=c(10, 5), geom="raster", drop=FALSE, size=0, linetype=0, fun="mean", size=0) + 
    ylim(0, 500) + 
    scale_fill_viridis_c(na.value = 0) + 
    xlab("Distance from motif center") + 
    ylab("Insert size") + 
    labs(fill="N insertions")
  
  if (isTRUE(all(d == round(d)))) {
    p = p + labs(fill="N insertions")
  } else {
    p = p + labs(fill="CPM")
  }
  
  return(p)
}


summary_tf_reads = function(d_tf, width_center=150) {
  
  if (is.factor(d_tf$insert_siz)) {
    lis = levels(d_tf$insert_size)
    nvars = strsplit(gsub("[(]", "", gsub("[]]", "", lis)), ",") %>% lapply(as.numeric)
    nval = sapply(nvars, max)
    names(nval) = lis
    d_tf$insert_size_numeric = nval[as.character(d_tf$insert_size)]
  } else {
    d_tf$insert_size_numeric = d_tf$insert_size
  }
  
  wh = d_tf$insert_size_numeric <= 100
  
  d_tf$abs_os = abs(d_tf$offset)
  
  d_tf$gr = case_when(d_tf$abs_os >= 800 ~ "background", d_tf$abs_os <= width_center ~ "signal")
  
  with(d_tf[wh,], tapply(freq, list(group=gr), sum)) %>% reshape::melt()
}


plot2d_tf_grid_comparison = function(d, idx1, idx2, tf, r, labels=c(), plot_data=TRUE, prior_count=0, tsse_values=c(), plot_ratio=FALSE) {
  
  collapse_d = function(d, idx) {
    
    if (length(idx) == 1) {
      res = d[[idx]][[tf]]
    } else {
      res = d[idx] %>%
        lapply(function(x) x[[tf]]) %>% 
        do.call(what=rbind) %>%
        plyr::count(c("offset","insert_size"), wt_var="freq")
    }
    
    if (is.factor(res$insert_size)) {
      res$insert_size_gr = res$insert_size
      res$insert_size = as.numeric(res$insert_size) * 5 - 2.5
    } else {
      cuts = seq(0, 10000, by=5)
      res$insert_size_gr = cut(res$insert_size, cuts, ordered_result=TRUE)
    }
    
    return(res)
  }
  
  collapse_r = function(r, idx) {
    
    if (length(idx) == 1) {
      return(r[[idx]])  
    }
    
    r[idx] %>%
      do.call(what=rbind) %>%
      plyr::count(c("insert_size","isize_group"), wt_var="freq")
  }
  
  get_d_matrix= function(d) {
    res = d %>% 
      dplyr::mutate(insert_size_gr=factor(insert_size_gr, vals_y)) %>%
      dplyr::mutate(offset=factor(offset, vals_x)) %>%
      dplyr::filter(!is.na(insert_size_gr)) %>% 
      reshape2::dcast(insert_size_gr~offset, value.var="freq", drop=FALSE) %>% 
      magrittr::set_rownames(., .$insert_size_gr) %>% 
      dplyr::select(-insert_size_gr) %>% 
      as.matrix()
    
    res[is.na(res)] = 0
    return(res)
  }
  
  adjust_by_reads = function(d, r, prior_count = 0, per_insert_size=TRUE) {
    
    if (!per_insert_size) {
      if (is.data.frame(r)) {
        r = sum(r$freq)
      }
    } else {
      stopifnot(is.data.frame(r))
      stopifnot(all(c("isize_group","freq") %in% colnames(r)))
      r = tapply(r$freq, r$isize_group, sum, na.rm = TRUE)
    }
    
    if (prior_count) {
      d = d + prior_count
      r = r + prior_count
    }
    
    if (per_insert_size) {
      d = apply(d, 2, "/", r[rownames(d)])
    } else {
      d = d / r
    }
    
    return(d)
  }
  
  
  
  # calculate ratios of data in plots
  vals_x = as.character(seq(-1000, 1000))
  cuts_y = seq(0, 10000, by=5)
  vals_y = as.character(unique(cut(25:500, cuts_y)))
  
  g1 = collapse_d(d, idx1)
  g2 = collapse_d(d, idx2)
  
  r1 = collapse_r(r, idx1)
  r2 = collapse_r(r, idx2)
  
  dm1_uadj = get_d_matrix(g1)
  dm2_uadj = get_d_matrix(g2)
  
  dm1 = adjust_by_reads(dm1_uadj, r1, prior_count, per_insert_size = TRUE)
  dm2 = adjust_by_reads(dm2_uadj, r2, prior_count, per_insert_size = TRUE)
  
  dm1o = adjust_by_reads(dm1_uadj, r1, prior_count, per_insert_size = FALSE)
  dm2o = adjust_by_reads(dm2_uadj, r2, prior_count, per_insert_size = FALSE)
  
  # remove background and center_data
  idx_background = abs(seq(-1000, 1000)) > 750
  
  bkg1 = apply(dm1[,idx_background], 1, mean)
  bkg2 = apply(dm2[,idx_background], 1, mean)
  
  rel_background = FALSE
  
  if (rel_background) { # signal less then background
    dm1b = dm1 / bkg1 #/ sf1 + max(c(max(bkg1) / sf1, max(bkg2) / sf2))
    dm2b = dm2 / bkg2 #/ sf2 + 
  } else { # signal above background, substract signal
    dm1b = dm1 - bkg1 #/ sf1 + max(c(max(bkg1) / sf1, max(bkg2) / sf2))
    dm2b = dm2 - bkg2 #/ sf2 + max(c(max(bkg1) / sf1, max(bkg2) / sf2))
  }
  
  #plot(apply(dm1b[1:10,], 2, mean))
  #plot(apply(dm2b[1:10,], 2, mean))
  #plot(apply(dm1b[1:10,], 2, mean) / apply(dm2b[1:10,], 2, mean), ylim=c(0, 10))
  #plot(apply(dm1b[1:10,], 2, mean) - apply(dm2b[1:10,], 2, mean))

  if (!plot_data) {
    return(list(d1_all_reads=dm1o, d2_all_reads=dm2o, d1_is_reads=dm1, d2_is_reads=dm2, d1_unadjusted=dm1_uadj, d2_unadjusted=dm2_uadj, dm1b=dm1b, dm2b=dm2b, reads1=r1, reads2=r2, bkg1=bkg1, bkg2=bkg2))
  }
  
  # plot data  
  p1 = plot2d_tf_grid(dm1) + ylim(25, 500) 
  p2 = plot2d_tf_grid(dm2) + ylim(25, 500)
  
  p1o = plot2d_tf_grid(dm1o) + ylim(25, 500) + ggtitle(labels[1])
  p2o = plot2d_tf_grid(dm2o) + ylim(25, 500) + ggtitle(labels[2])
  
  p1b = plot2d_tf_grid(dm1b) + ylim(25, 500) + ggtitle(labels[1])
  p2b = plot2d_tf_grid(dm2b) + ylim(25, 500) + ggtitle(labels[2])
  
  
  # squish data in the plots:
  wh1 = between(p1o$data$insert_size, 25, 300)
  wh2 = between(p2o$data$insert_size, 25, 300)
  max_cpm = max(c(quantile(p1o$data$freq_plot[wh1], 0.95), quantile(p2o$data$freq_plot[wh2], 0.95)))
  p1o = p1o + scale_fill_viridis_c(limits=c(0, max_cpm), oob=scales::squish) 
  p2o = p2o + scale_fill_viridis_c(limits=c(0, max_cpm), oob=scales::squish)
  
  wh1 = between(p1$data$insert_size, 25, 300)
  wh2 = between(p2$data$insert_size, 25, 300)
  max_cpm = max(c(quantile(p1$data$freq_plot[wh1], 0.95), quantile(p1$data$freq_plot[wh2], 0.95)))
  p1 = p1 + scale_fill_viridis_c(limits=c(0, max_cpm), oob=scales::squish) 
  p2 = p2 + scale_fill_viridis_c(limits=c(0, max_cpm), oob=scales::squish)
  
  wh1 = between(p1b$data$insert_size, 25, 300)
  wh2 = between(p2b$data$insert_size, 25, 300)
  max_cpm = quantile(c(p1b$data$freq_plot[wh1], p2b$data$freq_plot[wh2]), c(0, 0.95))
  p1b = p1b + scale_fill_viridis_c(limits=max_cpm, oob=scales::squish) + labs("rel. CPM")
  p2b = p2b + scale_fill_viridis_c(limits=max_cpm, oob=scales::squish) + labs("rel. CPM")
  
  
  p4 = r1 %>% 
    dplyr::filter(between(insert_size, 25, 500)) %>% 
    mutate(freq=freq/sum(freq)) %>% 
    ggplot(aes(x=insert_size, y=freq, color="Sample")) + 
    geom_line() + 
    geom_line(data=plyr::count(p1o$data, c("insert_size"), wt="freq_plot") %>% dplyr::filter(between(insert_size, 25, 500)) %>% mutate(freq=freq/sum(freq)), aes(x=insert_size, y=freq/5, color="TF site")) + 
    xlim(25, 500) +
    coord_flip() + 
    scale_color_brewer(palette = "Set1") + 
    labs(color="") + 
    xlab("Insert size") +
    ylab("Frequency") +
    scale_y_continuous(n.breaks = 2) #+ 
  #theme(legend.position = "bottom")
  
  
  p5 = ggplot(r2 %>% dplyr::filter(between(insert_size, 25, 500)) %>% mutate(freq=freq/sum(freq)), aes(x=insert_size, y=freq, color="Sample")) + 
    geom_line() + 
    geom_line(data=plyr::count(p2o$data, c("insert_size"), wt="freq_plot") %>% dplyr::filter(between(insert_size, 25, 500)) %>% mutate(freq=freq/sum(freq)), aes(x=insert_size, y=freq/5, color="TF site")) + 
    xlim(25, 500) +
    coord_flip() + 
    scale_color_brewer(palette = "Set1") + 
    labs(color="") + 
    xlab("Insert size") +
    ylab("Frequency") + 
    scale_y_continuous(n.breaks = 2) #+ 
  #theme(legend.position = "bottom")
  
  if (plot_ratio) {
    rel_data1 = log2(dm1 / dm2)
    rel_data2 = log2(dm1o / dm2o)
    limits1 = limits2 = c(-1, 1)
    label = "log2(A/B)"
  } else {
    rel_data1 = dm1 - dm2
    rel_data2 = dm1o - dm2o
    
    wh1 = between(p1$data$insert_size, 25, 300)
    wh2 = between(p2$data$insert_size, 25, 300)
    max_cpm = max(c(quantile(p1$data$freq_plot[wh1], 0.95), quantile(p2$data$freq_plot[wh2], 0.95))) * 0.25
    max_cpm2 = quantile(abs(rel_data1), 0.95)
    limits1 = min(c(max_cpm, max_cpm2)) * c(-1, 1)
    
    wh1 = between(p1o$data$insert_size, 25, 300)
    wh2 = between(p2o$data$insert_size, 25, 300)
    max_cpm = max(c(quantile(p1o$data$freq_plot[wh1], 0.95), quantile(p2o$data$freq_plot[wh2], 0.95))) * 0.25
    max_cpm2 = quantile(abs(rel_data2), 0.95)
    limits2 = min(c(max_cpm, max_cpm2)) * c(-1, 1)
    
    label = "A - B"
  }
  
  p3 =
    plot2d_tf_grid(rel_data1) +
    scale_fill_viridis_c(limits=limits1, oob=scales::squish, na.value = 0) + 
    labs(fill=label) + 
    ylim(25, 500) + 
    ggtitle(labels[3])
  
  p3o = 
    plot2d_tf_grid(rel_data2) +
    scale_fill_viridis_c(limits=limits2, oob=scales::squish, na.value = 0) + 
    labs(fill=label) + 
    ylim(25, 500) + 
    ggtitle(labels[3])


  if (rel_background) { # signal less then background
    
    p3b = 
      plot2d_tf_grid(log2(dm1b / dm2b)) +
      scale_fill_viridis_c(limits=c(-1, 1), oob=scales::squish, na.value = 0) + 
      labs(fill="log2(A/B)") + 
      ylim(25, 500) + 
      ggtitle(labels[3])
    
  } else { # signal above background, substract signal
    
    rel_val = (dm1b - dm2b)
    
    wh1 = between(p1b$data$insert_size, 25, 300)
    wh2 = between(p2b$data$insert_size, 25, 300)
    max_cpm = quantile(c(p1b$data$freq_plot[wh1], p2b$data$freq_plot[wh2]), 0.95) * 0.25
    max_cpm2 = quantile(abs(rel_val), 0.95)
    limits = min(c(max_cpm, max_cpm2)) * c(-1, 1)
    
    p3b = 
      plot2d_tf_grid(rel_val) +
      scale_fill_viridis_c(limits=limits, oob=scales::squish, na.value = 0) + 
      labs(fill="A - B") + 
      ylim(25, 500) + 
      ggtitle(labels[3])
    
  }
  

  
  
  if (length(tsse_values)) {
    
    names(tsse_values) = labels[1:2]
    for (i in seq_along(tsse_values)) {
      if (is.null(tsse_values[[i]]))
        tsse_values[[i]] = list(nf=NA, all=NA)
      for (j in seq_along(tsse_values[[i]])) {
        if (is.null(tsse_values[[i]][[j]])) {
          tsse_values[[i]][[j]] = NA
        }
      }
    }
    
    p6 = 
      reshape::melt(lapply(tsse_values, as.list)) %>%
      mutate(L1=factor(L1, labels[1:2], ordered = TRUE)) %>%
      mutate(L2=factor(L2, c("all","nf"), c("All", "IS <= 150"), ordered = TRUE)) %>%
      ggplot(aes(x=L2, y=value, group=L2)) + 
      geom_bar(stat="identity", width=0.6, fill="gray25") + 
      scale_fill_brewer(palette = "Set1") + 
      facet_wrap(~L1, ncol=1) + 
      ylab("TSS enrichment") + 
      xlab("") + 
      guides(fill=FALSE) + 
      theme(strip.text.x = element_text(size=9)) + 
      coord_flip()
    
  } else {
    
    p6 = ggplot(data.frame())
    
  }
  
  
  
  plot_rel_freq = function(d, r, bkg) {
    
    n_per_gr = tapply(r$freq, r$isize_group, sum, na.rm=TRUE)
    freq_per_gr1 = n_per_gr / sum(n_per_gr, na.rm=TRUE)
    
    r_tf = 
      d %>% 
      plyr::count(c("insert_size"), wt = "freq_plot") %>% 
      dplyr::filter(between(insert_size, 26, 500)) %>%
      mutate(freq=freq/sum(freq)) %>% 
      mutate(isize_gr = cut(insert_size, cuts_y)) %>% 
      mutate(freq_adj=freq/freq_per_gr1[as.character(isize_gr)]) %>% 
      mutate(freq_adj=freq_adj / mean(freq_adj)) %>% 
      mutate(tf=factor("TF", c("Sample","TF"), ordered = TRUE)) %>% 
      mutate(background=bkg[as.character(isize_gr)]) %>% 
      mutate(background=background / mean(background))
    
    ggplot(data=r_tf, aes(x=insert_size, y=freq_adj, color="TF")) + 
      geom_line() + 
      geom_line(aes(y=background, color="Background")) + 
      xlim(25, 500) +
      coord_flip() + 
      scale_color_brewer(palette = "Set1", drop=FALSE) + 
      labs(color="") + 
      xlab("Insert size") +
      ylab("Frequency") +
      geom_hline(yintercept = 1, linetype=2)
  }
  
  p7 = plot_rel_freq(p1o$data, r1, bkg1)
  p8 = plot_rel_freq(p2o$data, r2, bkg2)
  
  r_data = range(c(p7$data$freq_adj, p8$data$freq_adj))
  p7 = p7 + scale_y_log10(limits=r_data, breaks=sort(c(round(r_data, 2), 1)))
  p8 = p8 + scale_y_log10(limits=r_data, breaks=sort(c(round(r_data, 2), 1)))
  
  #theme(legend.position = "bottom")
  
  #plot(p7$data$freq, type="l")
  #abline(v=seq(1, 13) * 10.93 + 40, lty=3)

  title_pl = 
    with(calc_element("plot.title", theme_cowplot()), {
      ggdraw() + 
        draw_label(
          tf, fontfamily = family,
          fontface = face,
          size = size
        )
    })
  

  pe = ggplot() + theme_void()
  pl = list(p1o, p4, p1, p7, p1b, 
            p2o, p5, p2, p8, p2b, 
            p3o, p6, p3, pe, p3b)
  
  pg = cowplot::plot_grid(plotlist=pl, ncol=5, nrow=3, align = "h", axis="bt", rel_widths = c(1.5,1.0,1.5,1.0,1.5))
  
  #pl = list(p1o, p4, p1, p2o, p5, p2, p3o, p6, p3)
  #pg = cowplot::plot_grid(plotlist=pl, ncol=3, align = "h", axis="bt", rel_widths = c(1.5,1.0,1.5))
  pg2 = cowplot::plot_grid(title_pl, cowplot::as_grob(pg), ncol = 1, rel_heights = c(0.1,0.9))
  
  return(pg2)
}
 

calc_tss_enrichtment = 
  function(file, tss, genome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) {
    
    if (is.character(file)) {
      stopifnot(file.exists(file))
      chrs = paste0("chr", 1:22)
      d = read_shifted_count_file(as.character(file))
      d = d[as.character(seqnames(d)) %in% chrs]
    } else {
      stopifnot("GRanges" %in% class(file))
      d = file
    }
    
    frac = sum(width(reduce(tss))) / sum(seqlengths(genome)[chrs])
    ol_any = overlapsAny(d, tss)
    
    mean(ol_any) /  frac
  }


calc_tss_enrichtment = 
  function(file, tss, genome_size = 3257347282) {
    
    if (is.character(file)) {
      stopifnot(file.exists(file))
      chrs = paste0("chr", 1:22)
      d = read_shifted_count_file(as.character(file))
      d = d[as.character(seqnames(d)) %in% chrs]
    } else {
      stopifnot("GRanges" %in% class(file))
      d = file
    }
    
    chrs =  unique(as.character(seqnames(d)))
    frac = sum(width(reduce(tss))) / genome_size
    ol_any = overlapsAny(d, tss)
    
    mean(ol_any) /  frac
  }

dist_to_tss = 
  function(file, tss, genome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) {
    
    if (is.character(file)) {
      stopifnot(file.exists(file))
      chrs = paste0("chr", 1:22)
      d = read_shifted_count_file(as.character(file))
      d = d[as.character(seqnames(d)) %in% chrs]
    } else {
      stopifnot("GRanges" %in% class(file))
      d = file
    }
    
    d_closest = distanceToNearest(d, tss-(width(tss)-1)/2)
    dist = scales::squish(d_closest@elementMetadata$distance, c(0, 1e5))
    dplyr::count(data.frame(d = dist %/% 10 * 10), d)
    
  }

get_mark_dup_data = 
  function(x) {
    d = data.frame(readr::read_tsv(x, skip=6, n_max = 1, col_types=readr::cols()))
    d$SAMPLE = gsub("-GRCh38[.]metrics", "", basename(x))
    return(d)
  }


correct_tf_labels = function(x) {
  
  wh_corrected = grepl("_uc_", x)
  wh_to_correct = grepl("[no]PEAK", x)
  
  wh = grepl("unclustered", x)
  x[wh] = gsub("_unclustered_", "_uc_", x[wh])
  
  if (any(wh)) {
    wh_n_peak = grepl("nPEAK", x)
    wh_o_peak = grepl("oPEAK", x)
    x[wh_n_peak] = gsub("_nPEAK", "_oPEAK", x[wh_n_peak])
    x[wh_o_peak] = gsub("_oPEAK", "_nPEAK", x[wh_o_peak])
  } else {
    if (!any(wh_corrected) & any(wh_to_correct)) {
      stop("Unsure about correction status!")
    }
  }
  
  return(x)
}
