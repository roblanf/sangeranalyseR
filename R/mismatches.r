library(sangerseqR)
library(Biostrings)

# a script to calculate mismatches from ab1 chromatograms, based on a cutoff. Similar to what PHRED does.

process_abi_folder <- function(folder, cutoff){
  abi_files = list.files(folder, pattern = "\\.ab1$", full.names = T, recursive = T)
  lapply(abi_files, process_abi_file, cutoff = cutoff)
}


process_abi_file <- function(inputfile, cutoff = 0.33){
  
  # drop the output in the same folder as the inputfile
  outputfolder = dirname(inputfile)

  name = basename(inputfile)
  
  chromname = paste(name, "_", "chromatogram.pdf", sep='')
  tablename = paste(name, "_", "mismatches.csv", sep='')
  
  s = readsangerseq(inputfile)  
  basecalls = makeBaseCalls(s, ratio = cutoff)
  chromatogram(basecalls, height = 2, showcalls = 'both', filename = file.path(outputfolder, chromname))
  
  # align the primary and secondary base calls
  pa = pairwiseAlignment(primarySeq(basecalls), secondarySeq(basecalls), type = "global-local")
  
  # get number of disagreements
  mismatches = mismatchTable(pa)

  write.csv(mismatches, file = file.path(outputfolder, tablename))
  
}

