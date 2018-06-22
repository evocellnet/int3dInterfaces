library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
## intfcfile <- "../results/interfaces_notmapped.tab"
## mapingsfile <- "../data/pdb_mappings.rds"
intfcfile <- args[1]
mapingsfile <- args[2]
output <- args[3]

mapings <- readRDS(mapingsfile) %>%
    distinct() %>%
    mutate(pdbpos = as.numeric(pdbpos))

intfc <- read_tsv(intfcfile)

out <- intfc %>%
    inner_join(mapings, by = c("PROTEIN" = "acc",
                              "INTERFACE" = "file",
                              "CHAIN" = "pdbchain",
                              "POS_PDB" = "pdbpos",
                              "AA" = "uniprot_res")) %>%
    select(-POS, -pdb_res) %>%
    dplyr::rename(CHAIN_INT3D = chain, UNIPROT_POS = uniprot_pos)

write.table(out,
            file = output,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
    
