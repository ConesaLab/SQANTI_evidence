library(dplyr)
library(readr)
library(stringr)

# Args
args <- commandArgs(trailingOnly = TRUE)
gff_path <- args[1]
output <- args[2]
thr <- as.numeric(args[3])

gff_cols <- c("chr","source","type","start","end","score","strand","phase","attributes")

df <- read_tsv(gff_path, col_names = gff_cols)

# Print threshold
cat("Threshold given:", thr, "\n")

# Extract gene ids
df %>% 
    mutate(id = gsub(".*=","",str_split(attributes,";", simplify = TRUE)[,1])) -> df_with_id

initial_genes <- length(unique(df_with_id$id))
cat("Initial number of genes:", initial_genes, "\n")

# Summarise and filter
good_genes <- df_with_id %>%
    group_by(id) %>%
    summarise(type = paste0(unique(type),collapse=";"),
              attributes = first(attributes)) %>%
    ungroup() %>%
    mutate(identity = str_extract(attributes, "Identity=([^;]+)") %>% 
                        gsub("Identity=","",.) %>% 
                        as.numeric()) %>%
    filter(identity >= thr & str_detect(type,"stop_codon")) %>% 
    pull(id)

cat("Number of genes passing filter:", length(good_genes), "\n")

# Final genes in output
final_genes <- df_with_id %>%
    filter(id %in% good_genes)

cat("Final number of genes in output:", length(unique(final_genes$id)), "\n")

# Write output
final_genes %>%
    select(-id) %>%
    write_tsv(output,col_names=F)
