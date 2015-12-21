#' Check for secondary peaks in a folder of .ab1 files
#' 
#' This function recursively scans a folder that contains abi files either in itself or its subfolders, and looks for secondary peaks in each sequence. It produces an annotated plot of the secondary peaks above a user-defined cutoff, as well as a table of the mismatches.
#' 
#' It recursively scanes the input folder for all abi files. For each abi file, it creates and returns a data frame and a plot of mismatches, and also writes both to the same file as the .ab1 file (as a .csv and .pdf respectively).
#'
#' @param input.folder the folder in which to search recursively for *.ab1 files.
#' @param cutoff the ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are not. 
#' @param write.files TRUE/FALSE. FALSE (the default): do not write any files to disk. TRUE: write plots and .csv files of secondary peaks to disk, in this case, files are named with the same filename as the .ab1 files, but with .csv and .pdf extensions, and are saved in the same folder as the .ab1 file from which they were created.
#' 
#' @return output a dataframe of secondary peaks
#'
#' @keywords chromatogram, peak, mismatch
#'
#' @export secondary.peaks.abi.folder

secondary.peaks.abi.folder <- function(input.folder, cutoff, mc.cores = 1, write.files = FALSE){
  abi_files = list.files(input.folder, pattern = "\\.ab1$", full.names = T, recursive = T)
  mclapply(abi_files, process.abi.file, cutoff = cutoff, write.files = write.files)
}


