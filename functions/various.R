turn_off_clipping = function(x) {
  x = ggplot_gtable(ggplot_build(x))
  x$layout$clip = "off"
  
  wh = grep("strip-(t|r)", x$layout$name)
  for (i in wh) {
    x$grobs[[i]]$layout$clip = "off"
  }
  
  return(x)
}

convert_pdf_to_png = function(x) {
  out_dir = tempfile()
  dir.create(out_dir, FALSE, TRUE)
  out_file = file.path(out_dir, gsub("[.]pdf", ".png", basename(x)))
  cmd = paste0("convert -quality 100 -density 300 '", x, "' '", out_file, "'")
  system(cmd)
  return(out_file)
}


convert_png_to_pdf = function(x) {
  out_dir = tempfile()
  dir.create(out_dir, FALSE, TRUE)
  out_file = file.path(out_dir, gsub("[.]png", ".pdf", basename(x)))
  cmd = paste0("convert -quality 100 -density 300 '", x, "' '", out_file, "'")
  system(cmd)
  return(out_file)
}


merge_pdf = function(x, out_file) {
  in_files = paste0("'", x, "'", collapse = " ")
  cmd = paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile='", out_file, "' ", in_files)
  system(cmd)
  return(out_file)
}
