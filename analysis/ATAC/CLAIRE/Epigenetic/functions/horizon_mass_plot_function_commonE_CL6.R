plot_horizon_mass <- function(cov, patient, peak, smooth= 40, lwd=1.5, type="histogram", legend=TRUE, win_size=1000, patient_old=patient){
  library(Gviz)

  combos <- combos[grepl(patient,combos)]
  combos <- combos[grepl("_A|_B|_C|_D|C542_F",combos)]
  regions <- unlist(lapply(1: length(combos), function(c){unlist(str_split(combos[c],"_"))[2]}))
  regions <- c(regions,"E")
  col <- col[regions]
  
  peak_spl <- strsplit(peak, "[:-]")[[1]]
  chr <- peak_spl[1]
  start <- as.numeric(peak_spl[2])
  end <- as.numeric(peak_spl[3])
  peak_gr <- GRanges(chr, IRanges(start, end))
  
  #Get the order of presentation so that max value is plotted first
  cov <- c(cov[grepl(patient,names(cov))],cov["panpatientE"])
  cov <- cov[grepl("_A|_B|_C|_D|C542_F|panpatientE",names(cov))]
  cov<-lapply(cov, IRanges::subsetByOverlaps,peak_gr+win_size/2)
  suppressWarnings(max<-unlist(lapply(1:length(cov), function(n){
    c <- regions[n]
    max(cov[[c]]$score)})))
  max[which(!is.finite(max))] <- 0
  names(max) <- regions
  max_order <- order(max, decreasing=T)
  max_order <- names(max[max_order])
  
  
  
  # get tracks from coverage and colour them
  tracks <- lapply(1:length(regions), function(c) {
  
    mycov <- cov[[c]]
    
    track <- 
      DataTrack(
        genome = "hg38",
        range = mycov, 
        type = type,
        fill.mountain = "transparent",
        col.mountain = col,
        name = regions[c],
        fill = col,
        col.histogram = col,
        fill.histogram = "transparent",
        fill = "transparent",
        col = col,
        col.line = col,
        start = start - win_size/2,
        end = end + win_size/2,
        chr = chr
      )
    
    displayPars(track) <- 
      list(
        groups = factor(regions[c], levels = regions),
        col = col[regions[c]],
        legend = TRUE,
        missingAsZero = F,
        lwd = 3
      )
    
    return(track)})
  
  names(tracks) <- regions
  tracks <- tracks[max_order]
  names(tracks[[1]]) <- paste0(as.character(patient)," normalised coverage")
  ylims <- range(lapply(tracks,values))
  tracks_to_plot<-OverlayTrack(tracks,name=as.character(patient), ylim=ylims)
  annotracks <- peak_tracks[combos]
  
  pat_idx <- patient_old
  logfc <-round(fc[peak,pat_idx],digits = 2)
  sig2 <- sig_filter[peak,pat_idx]
  status <- get_status(patient=patient,   peak=peak)
  
  results_plot <- results[results$peak==peak, grepl(patient,colnames(results))]
  results_plot$CvN_pvalue <- p[peak,pat_idx]
  results_plot$CvN_padj <-  p_adj[peak,pat_idx]
  results_plot <- as.data.frame(sapply(results_plot, as.numeric))
  mynames <- rownames(results_plot)
  sig = sapply(results_plot[,1], format_p)
  
  results_plot <- as.data.frame(sapply(sapply(results_plot, as.numeric), scales::scientific))
  results_plot[,1] <- paste0(results_plot[,1],sig)
  rownames(results_plot) <- mynames
  P. <- c(results_plot[c(1,3,5,7),1])
  `Adjusted P.` <- c(results_plot[c(2,4,6,8),1])
  results_plot <- cbind(P.,`Adjusted P.`)
  rownames(results_plot) <- c("~region vs. ~","~purity vs ~","~purity+region vs ~purity","Carcinoma vs. normal")
  
  if (!is.na(sig2) & as.character(sig2)=="TRUE") {# & get_status(patient, peak)=="Clonal" ){
    
    main = paste0(
      peaks_original[peaks_original$peak==peak,"id"], " ",  peak, "\n",
      "Putatively ",get_status(patient,peak),", ", 
      "PASS=",as.character(sig2),", ",
      "LogFC=",round(fc[peak,pat_idx],digits = 2),", ",
      "~purity+region Padj=", results_plot[3,2]
    )
    
    plotTracks(
      background.title = "white",
      fontcolor.title = "black",
      col.axis = "black",
      showTitle = TRUE,
      clip = FALSE,
      c(
        AnnotationTrack(peak_gr, name = NULL, fill = "grey"),
        annotracks,
        tracks_to_plot
      ),
      legend = T,
      type = type,
      from = start - win_size/2,
      to = end + win_size/2,
      min.distance = 10,
      chromosome = chr,
      lwd.baseline = 1,
      col.baseline = "black",
      windowSize = 40,
      window = -1,
      main = main,
      ylim = ylims,
      cex.title = 1,
      cex.axis = 0.7,
      cex.main = 1,
      add = TRUE,
      sizes = c(0.02, rep(0.02, length(annotracks)), 1),
      detail = "coverage"
    )
}

}
