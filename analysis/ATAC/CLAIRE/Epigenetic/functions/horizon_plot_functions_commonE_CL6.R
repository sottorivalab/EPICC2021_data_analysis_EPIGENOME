
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

plot_horizon <- function(cov, patient, peak, smooth, lwd=4, type="histogram", legend=TRUE, win=0, include_adenoma=FALSE, patient_old=patient){
  
  combos <- combos[grepl(paste0(patient, "_[A-D,F-I]"),combos)]
  if (!include_adenoma) combos <- combos[grepl("_A|_B|_C|_D|C542_F",combos)]
  regions <- unlist(lapply(1: length(combos), function(c){unlist(str_split(combos[c],"_"))[2]}))
  regions <- c(regions,"E")
  col <- col[regions]
  
  peak_spl <- strsplit(peak, "[:-]")[[1]]
  chr <- peak_spl[1]
  start <- as.numeric(peak_spl[2])
  end <- as.numeric(peak_spl[3])
  peak_gr <- GRanges(chr, IRanges(start, end))
  
  #Get the order of presentation so that max value is plotted first
  cov <- c(cov[combos],cov["panpatientE"])
  names(cov) <- regions # todo maybe fix this?
  if (any(sapply(cov, is.null))) 
    stop("Missing some cov data for regions: ", paste0(names(cov)[sapply(cov, is.null)], collapse=", "))
  cov<-lapply(cov,subsetByOverlaps,peak_gr+win/2) 
  suppressWarnings(max<-unlist(lapply(1:length(regions),function(n){
    c <- regions[n]
    max(cov[[c]]$score)})))
  
  max[which(!is.finite(max))] <- 0
  names(max) <- regions
  max_order <- order(max, decreasing=T)
  max_order <- names(max[max_order])
  
  # get tracks from coverage and colour them
  tracks <- lapply(1:length(regions), function(c) {
    region <-regions[c]
    mycov <- cov[[c]]
    track <- Gviz::DataTrack(range=mycov,
                       fill.mountain="transparent",
                       col.mountain=col,
      name=combos[c],fill=col,
      col.histogram=col,
      fill.histogram="transparent",
      fill="transparent",
      col=col,
      col.line=col,
      from=start- win/2,
      to = end + win/2,
      chr=chr
      )
    Gviz::displayPars(track) <- list(groups = factor(region, 
                               levels = regions),
                               col=col[region],
                               legend = TRUE,
                               missingAsZero=F,
                               lwd=3)
    return(track)})
  
  names(tracks) <- regions
  
  tracks <- tracks[max_order]
  
  names(tracks[[1]]) <- paste0(as.character(patient)," normalised coverage")
  
  ylims <- range(lapply(tracks,values))
  
  tracks_to_plot<-Gviz::OverlayTrack(tracks,name=as.character(patient),
                               ylim=ylims)
  
  annotracks <- peak_tracks[c(combos,"pan_E")]
  

  Gviz::plotTracks(background.title="white",fontcolor.title="black",col.axis="black",
             showTitle=T, clip=FALSE,
             c(Gviz::AnnotationTrack(peak_gr, name = NULL, fill="grey"),
               annotracks,
               tracks_to_plot),
             legend=legend,
             type="l",
             from=start-win/2,
             to = end+win/2,
             min.distance=10,
             chromosome=chr,
             lwd.baseline=1,
             col.baseline="black",
             windowSize =40,window=-1,
             main=paste(peaks[peaks$peak==peak,"id"],
                        peaks[peaks$peak==peak,"type"] ,
                        peaks[peaks$peak==peak,"event_type"],
                        peak, 
                        sep=" "),
             ylim=ylims,cex.title=1,cex.axis=0.7, cex.main=1, add=TRUE,
             sizes=c(0.02,rep(0.02,length(combos)+1),1),
             detail="coverage")
}



col <- RColorBrewer::brewer.pal(name="Set1", n=9)
names(col) <- c("A","B","C","D","E","F","G","H","I")


plot_purity <- function(patient=patient, include_adenoma=FALSE){
  purity_plot <- coldata[coldata$patient==patient & (coldata$subtissue=="cancer" | include_adenoma),]
  purity_plot$region <- as.factor(as.character(purity_plot$region)) 
  
  ggplot(purity_plot, aes(region,purity, fill=region), ylim=c(0,1)) + 
    geom_boxplot() +
    theme_classic() +
    scale_y_continuous(limits = c(0,1.0), expand=c(0,0),labels = scales::percent) +
    theme(axis.line.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = col)
}


Pval_table <- function(patient=patient,peak=peak, patient_old=patient, include_adenoma=FALSE){
  
  wh_cols <- c("purity_only_nf_pvalue", "purity_only_nf_padj", "region_only_nf_pvalue", "region_only_nf_padj", "region_nf_pvalue", "region_nf_padj")
  results_plot <- unlist(results[results$peak==peak, -1])[paste0(patient, "_", wh_cols)]
  results_plot <- c(results_plot, p[peak, patient_old])
  results_plot <- c(results_plot, p_adj[peak, patient_old])
  results_plot<-as.data.frame(sapply(results_plot, as.numeric))
  mynames <- rownames(results_plot)
  sig <- format_p(results_plot[,1])
  results_plot <- as.data.frame(sapply(sapply(results_plot, as.numeric), scales::scientific))
  results_plot[,1] <- paste0(results_plot[,1],sig)
  rownames(results_plot) <- mynames
  P. <- c(results_plot[c(1,3,5,7),1])
 `Adjusted P.` <- c(results_plot[c(2,4,6,8),1])
  results_plot <- cbind(P.,`Adjusted P.`)
  rownames(results_plot) <- c("~region vs. ~","~purity vs ~","~purity+region vs ~purity","Tumour vs. normal")
  
  title <- textGrob(paste0(ifelse(sig[6]=="","Putatively Clonal","Putatively Subclonal"),
                           " (logFC = ",
                           round(fc[peak,patient_old],digits = 2),
                           ")"),
                    gp=gpar(fontsize=12))
  padding <- unit(5,"mm")
  t1<-gridExtra::tableGrob(results_plot, theme=gridExtra::ttheme_default(base_size=10))
  table <- gtable::gtable_add_rows(
    t1, 
    heights = grobHeight(title) + padding,
    pos = 0)
  table <- gtable::gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))
  grid.draw(table)
}



get_status <- function(patient=patient,peak=peak){
  p_val <- results[results$peak==peak, paste0(patient, "_region_nf_padj")]
  sig <- format_p(p_val)
  title <- ifelse(sig=="","Clonal","Subclonal")
  return(title)
}


plot_info <- function(pdf_name, patient, peak, type, patient_old=patient, include_adenoma=FALSE,...){
  pdf(pdf_name, height=5, width=9.8)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(3, 2, widths = c(6,3.8))))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  pushViewport(viewport(layout.pos.row=1:3,layout.pos.col=1))
  print(
    plot_horizon(
      cov = cov,
      patient = gsub("-.*", "", patient),
      patient_old = patient_old,
      peak = peak,
      smooth = 30,
      lwd = 1.5,
      legend = T,
      type = "l",
      include_adenoma = include_adenoma,
      ...
    ),
    vp = vplayout(x = 1:3, y = 1)
  )
  popViewport()
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=2))
  print(Pval_table(patient=patient,peak=peak, patient_old=patient_old),vp = vplayout(1, 2))
  popViewport()
  pushViewport(viewport(layout.pos.row=2:3,layout.pos.col=2))
  p1<-plot_purity(patient=gsub("-.*", "", patient), include_adenoma = include_adenoma)
  print(p1,vp = vplayout(2:3,2))
  popViewport()
  graphics.off()
}


plot_nopurity <- function(pdf_name, patient, peak, type){
  pdf(pdf_name, height=6.5, width=6)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2, 1, heights = c(4.7,1.8))))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=1))
  print(plot_horizon(cov=cov,
                     patient=patient,
                     peak=peak,
                     smooth=30,
                     lwd=1.5,
                     legend=T,
                     type="l"), vp= vplayout(x=1,y=1))
  popViewport()
  pushViewport(viewport(layout.pos.row=2,layout.pos.col=1))
  print(Pval_table(patient=patient,peak=peak),vp = vplayout(2, 1))
  popViewport()
  graphics.off()
  
}

format_p = function(x){
  case_when(
    x <= 0.0001 ~ "****",
    x <= 0.001 ~ "***",
    x <= 0.01 ~ "**",
    x <= 0.05 ~ "*", 
    TRUE ~ ""
  )
}
