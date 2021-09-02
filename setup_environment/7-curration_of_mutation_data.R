# this is a list of variants that are to be treated as 
# clonal (i.e., '.variants_clonal') or subclonal (i.e., '.variants_drop')
# since they were identified based on the lifted dN/dS datasets, they are 
# all hg19! They are lifted to the ids in hg38 below.

.variants_drop = c(
  "C516.cancer-chr12:115112326_C/T",
  "C562.cancer-chrY:19708368_C/T",
  "C516.cancer-chr17:7577139_G/A",
  "C516.cancer-chr19:11143966_G/A",
  "C516.cancer-chr4:126241566_G/A",
  "C516.cancer-chr9:96055271_G/A",
  "C525.cancer-chr7:151945334_T/C",
  "C543.cancer-chr7:151945334_T/C",
  "C551.cancer-chr12:115110005_G/A",
  "C561.cancer-chr13:25016723_C/T",
  "C525.cancer-chr7:151945334_T/C", 
  "C549.cancer-chr7:151932897_C/G", 
  "C562.cancer-chrY:21870254_C/T",
  "C516.cancer-chr11:92534000_C/A", 
  "C549.cancer-chr7:151935831_G/A", 
  "C550.cancer-chr7:151935831_G/A", 
  "C555.cancer-chr7:151945101_G/C", 
  "C560.cancer-chr7:151935831_G/A",
  "C516.cancer-chr12:115112153_GC/G", # adenoma contamination
  "C516.cancer-chr13:49030367_G/GA",
  "C562.cancer-chr20:35663817_T/TG",
  "C551.adenoma-chr11:108160340_GA_G"
)

.variants_clonal = c(
  "C516.cancer-chr16:9984856_C/T",
  "C516.cancer-chr2:141356219_T/C",
  "C516.cancer-chr3:10102088_G/C",
  "C525.cancer-chr13:25052393_G/T",
  "C536.adenoma-chrX:39931622_G/A",
  "C551.cancer-chr3:30686344_G/T",
  "C552.cancer-chr1:27023010_C/A",
  "C554.cancer-chr17:29562775_G/A",
  "C560.cancer-chr5:112128143_C/T",
  "C530.cancer-chr12:25398284_C/T",
  "C532.cancer-chr3:178936091_G/A",
  "C528.cancer-chr3:53902821_T/G", 
  "C531.cancer-chr3:47164925_G/A", 
  "C551.cancer-chr10:70406352_C/T", 
  "C561.adenoma (G)-chr5:176719141_C/G", 
  "C562.cancer-chr8:41790858_G/A",
  "C516.cancer-chr4:153247164_T/C", 
  "C552.cancer-chr2:141055396_G/A",
  "C524.cancer-chr17:12032541_C/T", # subclonal LOH 
  "C524.cancer-chr2:148654004_T/C", # subclonal LOH 
  "C527.cancer-chr5:112175153_G/T", # subclonal LOH
  "C537.cancer-chr3:89528584_C/G",  # subclonal LOH 
  "C549.cancer-chr4:126329733_C/T",  # subclonal LOH 
  "C560.cancer-chr5:112175939_GA/G",  # subclonal LOH
  "C538.cancer-chr5:112175211_TAAAAG/T",  # subclonal LOH
  "C516.cancer-chr14:105932775_G/GC",  # subclonal LOH
  "C525.cancer-chr18:3457667_A/ATGTCCCTTTCTCTCTCTGCCGAC"
)


# lift the variants above to hg38
lift_ids_removed = function(x) {
  type_data = gsub("-.*", "", x)
  mut_data = gsub(".*_", "", x)
  pos_data = gsub("_.*", "", gsub(".*-", "", x))
  pos_data_gr = as(pos_data, "GRanges")
  hg38_pos = THmisc::lift_hg19_to_hg38(pos_data_gr)
  paste0(type_data, "-", as.character(hg38_pos), "_", mut_data) 
}

.variants_drop_hg38 = lift_ids_removed(.variants_drop)
.variants_clonal_hg38 = lift_ids_removed(.variants_clonal)
