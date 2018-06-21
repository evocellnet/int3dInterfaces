library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
intfcfile <- "../results/interfaces_notmapped.tab"
mappingsfile <- "../data/pdb_mappings.rds"
intfcfile <- args[1]
mappings <- args[2]
output <- args[3]

mapppings <- readRDS(mappingsfile) %>% distinct()
intfc <- read_tsv(intfcfile)

intfc %>%
    left_join(mappings, by = c("PROTEIN" = "acc",
                               "INTERFACE" = "file",
                               "CHAIN" = "pdbchain",
                               "POS_PDB" = "pdbpos")) %>%
    select(-POS) %>%
    rename(CHAIN_INT3D = chain, POS = uniprot_pos)
    
