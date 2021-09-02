low_pass_parse_site_annotation = function(site) {
  do.call(rbind, strsplit(site, "[:-]")) %>% 
    magrittr::set_colnames(c("chromosome","start.pos","end.pos")) %>% 
    data.frame(stringsAsFactors = FALSE) %>% 
    mutate(chromosome=paste0("chr", chromosome)) %>%
    mutate(start.pos=as.numeric(start.pos)) %>% 
    mutate(end.pos=as.numeric(end.pos))
    
}


load_low_pass_segment_files = function(f) {
  
  data =
    readr::read_tsv(
      file = f,
      col_names = c("site", "logR"),
      skip = 1,
      col_types = "cd",
      progress = FALSE
    ) %>% mutate(depth.ratio = 2 ^ logR)
  
  site_annotation = low_pass_parse_site_annotation(data$site)
  invisible(cbind(site_annotation, data))
}


load_low_pass_custom_call_files = function(f) {
  
  data = 
    readr::read_tsv(
      file = f,
      col_names = c("site", "CNt"),
      skip = 1,
      col_types = "cd",
      progress = FALSE
    )
  
  site_annotation = low_pass_parse_site_annotation(data$site)
  invisible(cbind(site_annotation, data))
}

