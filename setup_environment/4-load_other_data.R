# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Load TSS data ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tsg_mapper = readr::read_tsv("external_datasets/ENSTandENSG.txt", col_types = readr::cols())
tsg_mapper = tsg_mapper$`Gene stable ID` %>% magrittr::set_names(tsg_mapper$`Transcript stable ID`)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
all_ts = transcripts(txdb)
tss_within = 1000

annodb = "org.Hs.eg.db"

#
all_ts$tx_name = gsub("[.][[:digit:]]*$", "", all_ts$tx_name)
all_ts$tg_name = tsg_mapper[all_ts$tx_name] %>%  magrittr::set_names(NULL)

# find tts sites:
tss_start = ifelse(as.character(strand(all_ts)) == "+", start(all_ts), end(all_ts))
.tss_pos = GRanges(seqnames(all_ts), IRanges(tss_start, tss_start)) + tss_within
.tss_pos$tx_name = all_ts$tx_name 
.tss_pos$tg_name = all_ts$tg_name 


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Genehancer DB ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

tryCatch({
  
  .genehancer_data = 
    "external_datasets/genehancer_tracks" %>%
    list.files("[.]bed[.]gz", full.names=TRUE) %>% 
    (function(x) x[grepl("DoubleElite", x)]) %>%
    magrittr::set_names(., gsub("(DoubleElite)?[.]bed[.]gz","", basename(.))) %>%
    lapply(readr::read_tsv)
  
  .longer_enh_label = 
    with(.genehancer_data$geneHancerInteractions, {
      
      get_gene_list = function(x) {
        x = unique(sort(x))
        if(length(x) < 5) x 
          else c(head(x, 4), "...")
      }
      
      gpe = split(geneName, geneHancerIdentifier) %>% lapply(get_gene_list)
      e2g = sapply(gpe, paste, collapse=", ")
      long_label = paste0(e2g, " (", names(e2g),")")
      names(long_label) = names(e2g)
      return(long_label)
    })
  
}, error=function(e) print("Couldn't load genehancer dataset."))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# List of reference ids ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

.reference_sample_ids = 
  "created_datasets/reference_samples.csv" %>% 
  readr::read_csv() %>% 
  unlist() %>% 
  magrittr::set_names(NULL)

.msi_positiv = c("C536","C548","C516","C518","C552","C562") 
.msi_positiv_adenoma = c("C516") 

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Colors ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

color_data_type = c(
  WGS = "#440154ff",
  ATAC = '#21908dff',
  RNA = '#fde725ff',
  LP = 'gray50',
  'lowpass-WGS' = 'gray50',
  'ATAC-seq' = '#21908dff',
  'RNA-seq' = '#fde725ff'
)

color_tissue_type = c(
  cancer = "#FF1800",
  adenoma = "#008199",
  normal = '#6FB400'
)

region_colors = c(
  A = "#E41A1C",
  B = "#377EB8",
  C = "#4DAF4A",
  D = "#984EA3",
  E = "#FF7F00"
)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if (require("DESeq2")) {
  .data_rnaseq_normalised = readRDS(file.path(.pipeline_ddir, "RNA", "data_rnaseq_normalised.rds"))
  colnames(.data_rnaseq_normalised) = paste0("EPICC_", colnames(.data_rnaseq_normalised), "_R1")
  .annot_rnaseq_normalised = THmisc::annotation_from_barcode(colnames(.data_rnaseq_normalised))
  .dds_data = readRDS(file.path(.pipeline_ddir, "RNA", "dds.expanded.filgenes.ensembl.rds"))
}
