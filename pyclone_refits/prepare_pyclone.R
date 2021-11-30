
format_to_pyclone <-  function(xx,purity){

  library(dplyr)
  
  
  # Format data for pyclone, we are using everything as it was clonal
  mutation_id <- paste0("chr",xx$Chromosome, ":", xx$Start_position, ":", xx$Start_position, ":", round(xx$ccf,3), ":", xx$gene, ":", xx$isdriver)
  ref_counts <-  xx$t_ref_count
  alt_counts <- xx$t_alt_count
  sample_id <- xx$Tumor_Sample_Barcode
  major_cn <- xx$major_cn
  minor_cn <-  xx$minor_cn
  normal_cn <-  2
  tumour_content <- purity
  gene_symb <- xx$gene

  df <-  data.frame(
    mutation_id = mutation_id,
    ref_counts = ref_counts,
    alt_counts = alt_counts,
    sample_id = sample_id,
    major_cn = as.integer(major_cn),
    minor_cn = as.integer(minor_cn),
    normal_cn = as.integer(normal_cn),
    tumour_content = tumour_content,
    gene_symb = gene_symb
  )
  
  # Before returning the dataframe we filter NAs in CNV and NV/NR values
  return(df %>%  filter(!is.na(major_cn), !is.na(ref_counts), !is.na(alt_counts)) %>% unique() %>%
           filter(!duplicated(mutation_id)) %>%  unique())

}


# We assume to have a csv with all the PCWAG mutation assignment files binded
# snvs <-  data.table::fread("maf_ccf_withdriver.csv", data.table = F)


# Split the dataset in a list with SNVs values for each sample and add names
snvs_list <- snvs %>% filter(Variant_Type == "SNP") %>%  group_split(Tumor_Sample_Barcode)

names <- sapply(snvs_list, function(x) x$Tumor_Sample_Barcode %>%  unique())

names(snvs_list) <-  names

# We also need purity values for each sample
# They can be donwloaded from https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_cnv/consensus.20170217.purity.ploidy.txt.gz

#meta <- lapply(names, function(x) CDSLabData::load_dataset(type = "WGS", sample = x, cohort = "PCAWG")$metadata$purity )

meta <-  meta %>% as.numeric()

names(meta) <-  meta

# Format data for pyclone-vi

py_clone_VI <-  mapply(snvs_list, meta,FUN = function(x,y) format_to_pyclone(x, y), SIMPLIFY = F)
names(py_clone_VI) <-  names(snvs_list)


# Generate a structure with a directory for each sample
mapply(py_clone_VI, names(py_clone_VI), FUN = function(x,y) {

  dir.create(y, showWarnings = F)
  write.table(x, file = paste0("./", y, "/pyclone_input.tsv"), sep = "\t", row.names = F, quote = F)

})

# save the name of the samples that we used
samples <-   names(py_clone_VI)
save(names,file =  "samples.rda")



