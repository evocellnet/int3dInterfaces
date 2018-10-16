library(tidyverse)
library(bio3d)
library(seqinr)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)
interactionsFile <- args[1]
pdbdir <- args[2]
uniprotfastafile <- args[3]
outfile <- args[4]

getChainSeqs <- function(filename){
    invisible(capture.output(pdb <- suppressWarnings(read.pdb(file = filename))))
    pdbname <- basename(filename)
    message("\t* Reading ",pdbname)
    chains <- unique(pdb$atom$chain)
    seq1 <- pdbseq(bio3d::trim(pdb, chain = chains[1]))
    seq2 <- pdbseq(bio3d::trim(pdb, chain = chains[2]))
    out <- bind_rows(
        data.frame(pdb_aa = seq1,
                   pdb_pos = as.numeric(names(seq1)),
                   chain = chains[1],
                   file = pdbname,
                   PROTEIN_CHAIN_LENGTH = length(seq1),
                   PARTNER_CHAIN_LENGTH = length(seq2),
                   stringsAsFactors = FALSE),
        data.frame(pdb_aa = seq2,
                   pdb_pos = as.numeric(names(seq2)),
                   chain = chains[2],
                   file = pdbname,
                   PROTEIN_CHAIN_LENGTH = length(seq2),
                   PARTNER_CHAIN_LENGTH = length(seq1),
                   stringsAsFactors = FALSE))
    return(out)
}

message("- Reading files...")

#Read interactions file
interactionsFile <- read_tsv(interactionsFile)

##Read uniprot fasta
up_fa <- readAAStringSet(uniprotfastafile)
names(up_fa) = sapply( strsplit(names(up_fa), '\\|'), function(i) i[2] )

#Extract all pdb residues
allpdbresidues <- interactionsFile %>%
    pull(FILENAME) %>%
    map_df(~getChainSeqs(file.path(pdbdir, .))) %>%
    left_join(bind_rows(interactionsFile %>%
                        select(acc = PROT1, file = FILENAME, pdbchain = CHAIN1) %>%
                        mutate(chain = "A"),
                        interactionsFile %>%
                        select(acc = PROT2, file = FILENAME, pdbchain = CHAIN2) %>%
                        mutate(chain = "B")),
              by = c("file", "chain"))

message("- Pairwise alignments...")

#sequences to align
toalign <- allpdbresidues %>%
    group_by(acc, file, chain, pdbchain, PROTEIN_CHAIN_LENGTH, PARTNER_CHAIN_LENGTH) %>%
    arrange(as.numeric(pdb_pos)) %>%
    summarise(pdbseq = paste(pdb_aa, collapse = ""),
              pdbpos = paste(pdb_pos, collapse = ";")) %>%
    ##excluding protein accessions without uniprot fasta
    filter(acc %in% names(up_fa)) %>%
    ungroup()

#alignment
aln <- pairwiseAlignment(up_fa[toalign$acc], toalign$pdbseq)

#Extracting alignment positions
message("- Extracting alignment positions...")
dat <- lapply(1:length(aln), function(i){
  if(i %% 100 == 0) message(i)
  up = strsplit( as.character( pattern(aln)[i] ), '')[[1]]
  pdbaln =  strsplit( as.character( subject(aln)[i] ), '')[[1]]
  df = data.frame(uniprot_res = up, pdb_res = pdbaln, stringsAsFactors = FALSE) %>%
      mutate(uniprot_pos = NA, pdbaln_pos = NA)
  
  ind <- up != '-'
  df$uniprot_pos[ ind ] <- start(pattern(aln)[i]):(start(pattern(aln)[i]) - 1 + sum(ind))
  
  ind <- pdbaln != '-'
  df$pdbaln_pos[ ind ] = start(subject(aln)[i]):(start(subject(aln)[i]) - 1 + sum(ind))

  df$acc = toalign$acc[i]
  df$file = toalign$file[i]
  df$chain = toalign$chain[i]  
  df
})
dat <- bind_rows(dat)

#Extracting positions in pdb
pdb_pos_aln_mapping <- toalign %>%
    dplyr::select(-pdbseq) %>%
    mutate(pdbpos = str_split(pdbpos, ";")) %>%
    mutate(newposmap = map_int(pdbpos, length)) %>%
    mutate(pdbaln_pos = map(newposmap, function(x) 1:x)) %>%
    unnest() %>%
    dplyr::select(-newposmap)

message("- Position mapping...")
#Mapping 
dat <- dat %>%
    inner_join(pdb_pos_aln_mapping, by = c("acc", "file", "chain", "pdbaln_pos")) %>%
    select(-pdbaln_pos)

message("- Printing...")
saveRDS(dat, file = outfile)
