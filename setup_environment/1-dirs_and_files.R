#---------------------------------------------
# Various file locations
#---------------------------------------------

# liftover file locations:
options(chain_hg38_to_hg19="external_datasets/lift_over_chains/hg38ToHg19.over.chain")
options(chain_hg19_to_hg38="external_datasets/lift_over_chains/hg19ToHg38.over.chain")

# location of function source files
.function_source_file_dir = "functions"

#---------------------------------------------
# raw data file locations
#---------------------------------------------
  
.list_dedup_name = function(d, s, pf="") {
  # small helper function to list files and name these 
  list.files(d, s, full.names=T) %>% 
    sort(decreasing = TRUE) %>% (function(x) x[!duplicated(basename(x))]) %>% # this gets rid of duplicated fits (e.g. for CN)
    magrittr::set_names(gsub(paste0("^", pf), "", gsub(paste0(s,"$"), "", basename(.))))
}


# vcf locations
vcf_subdir = c("VCF/filtered", "VCF/unfiltered")
.vcf_data_dir = file.path(.pipeline_ddir, vcf_subdir) 
.vcfs = .list_dedup_name(.vcf_data_dir, "_all_samples_filtered_annotated[.]vcf[.]gz$")
.vcfs_uf = .list_dedup_name(.vcf_data_dir, "_all_samples_unfiltered_annotated[.]vcf.gz$")
checkmate::assertSetEqual(names(.vcfs), names(.vcfs_uf))


# sequenza file locations
sequenza_subdir = "CNA/sequenza"
.sequenza_ddirs = file.path(.pipeline_ddir, sequenza_subdir)
.sequenza_fit_files = .list_dedup_name(.sequenza_ddirs, "_GRCh38_confints_CP[.]txt$")
.sequenza_cna_files = .list_dedup_name(.sequenza_ddirs, "_GRCh38_segments[.]txt$")
checkmate::assertTRUE(all(names(.sequenza_fit_files) %in% names(.sequenza_cna_files)))


# qdnaseq lowpass files
.lp_cn_ddir = file.path(.pipeline_ddir, "CNA/qdnaseq")
.lp_metric_files = .list_dedup_name(.lp_cn_ddir, "_500kb_GRCh38_multiregion_metrics[.]txt", "C[0-9]+_")
.lp_segment_files = .list_dedup_name(.lp_cn_ddir, "_500kb_GRCh38_multiregion_cna_segments[.]txt", "C[0-9]+_")
.lp_cn_files = .list_dedup_name(.lp_cn_ddir, "_500kb_GRCh38_multiregion_segmentation_calls[.]txt", "C[0-9]+_")


# purity estimates:
.gt_purity_esimates = "created_datasets/genotyping_estimates_per_sample.rds" # "analysis/ATAC & WGS/datasets/genotyping_estimates_per_sample.rds"


# genotyping of SNVs
.gt_dir = file.path(.pipeline_ddir, c("SNVs/wgs_and_lp", "SNVs/atac"))
.gt_files = list.files(.gt_dir, "[.]table[.]gz$", full.names=TRUE, recursive=TRUE)
names(.gt_files) = gsub(".*_vs_", "", gsub("-GRCh38.*", "", basename(.gt_files)))


# other dirs/files:
.coverage_data_dir = 
  file.path(.pipeline_ddir, "qualimap")

.msisensor_files = 
  file.path(.pipeline_ddir, "msisensor") %>% 
  .list_dedup_name("(_GRCh38)?", "C[0-9]+_") %>% 
  (function(x) x[file.size(x) > 0]) %>% 
  (function(x) x[!duplicated(names(x))])
  

# atac-seq read files:
.atac_bed_dir = file.path(.pipeline_ddir, "atac", "reads")
.regex_cutsites = "_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz$"
.atac_bed_files = list.files(.atac_bed_dir, .regex_cutsites, full.names=TRUE)
names(.atac_bed_files) = gsub(.regex_cutsites, "", basename(.atac_bed_files))
