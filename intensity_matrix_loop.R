#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

# Load package
library("Cardinal")

# file name
folder <- args[1] # "tmp/"
filename <- args[2] # 20190523_MS37_MetaboliteMix_Spots_60_900_SDHB_pos_A25_75um_88x70
resolution <- args[3] # 0.001 
units <- args[4] # "mz"
range <- args[5] # 50000 (describing the amount of bins that is processed at a time)
i_threshold <- args[6] # 100 (for not normalized datasets)

#-------------------------------------------------------------

# Import MSi dataset (ImzML and ibd!)
msi <- readImzML(folder,
                 filename,
                 attach.only = FALSE,
                 as="MSImagingExperiment",
                 resolution = resolution, 
                 units = units) 

# define fragmentation
pixel       <- as.character(seq(1,ncol(msi),1))
mz          <- as.vector(mz(msi))

# create range list
range_list <- list()
range_list[[1]] <- 1:range
for (i in 2:floor(length(mz)/range)) {
  range_list[[i]] <- range_list[[i-1]]+range
}
range_list[[length(range_list)+1]] <- (floor(length(mz)/range)*range+1):length(mz)

# create mz list
mz_list <- list()
for (i in 1:length(range_list)) {
  mz_list[[i]] <- mz[range_list[[i]]]
}

# create minor matrices
matrix_list <- list()
for (i in 1:length(mz_list)) {
  matrix_list[[i]]            <- spectra(msi)[range_list[[i]],]
  colnames(matrix_list[[i]])  <- pixel
  rownames(matrix_list[[i]])  <- mz_list[[i]]
  matrix_list[[i]]            <- matrix_list[[i]][rowSums(matrix_list[[i]] > i_threshold) >= 1, ] # strips matrix by deleting peaks (rows) under threshold
  matrix_list[[i]]            <- t(matrix_list[[i]])
  mz_list[[i]]                <- as.numeric(colnames(matrix_list[[i]]))
  matrix_list[[i]]            <- rbind(mz_list[[i]], matrix_list[[i]])
}

# fuse matrices
rm(msi, mz, range_list, mz_list)
intensity_matrix <- do.call(cbind, matrix_list)
rm(matrix_list)

# Save data as *.csv
write.table(intensity_matrix,
            file = paste(filename, ".csv"), # as previous filename
            row.names = TRUE,
            na = "",
            col.names = FALSE,
            sep = ",")




