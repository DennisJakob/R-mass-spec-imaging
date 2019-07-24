# Load package
library("tictoc")
library("Cardinal")

# file name
file <- "tmp/20190523_MS37_MetaboliteMix_Spots_60_900_SDHB_pos_A25_75um_88x70" #args[1]

#-------------------------------------------------------------

# Import MSi dataset (ImzML and ibd!)
msi <- readImzML(file,
                 attach.only = FALSE,
                 as="MSImagingExperiment",
                 resolution = 0.001, #args [2]
                 units = "mz") #args [3]

# define fragmentation
pixel       <- as.character(seq(1,ncol(msi),1))
mz          <- as.vector(mz(msi))
range       <- 50000 #args [4]

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
i_threshold <- 100 #args [5]

matrix_list <- list()

for (i in 1:length(mz_list)) {
  matrix_list[[i]]            <- spectra(msi)[range_list[[i]],]
  colnames(matrix_list[[i]])  <- pixel
  rownames(matrix_list[[i]])  <- mz_list[[i]]
  matrix_list[[i]]            <- matrix_list[[i]][rowSums(matrix_list[[i]] > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
  matrix_list[[i]]            <- t(matrix_list[[i]])
  mz_list[[i]]                <- as.numeric(colnames(matrix_list[[i]]))
  matrix_list[[i]]            <- rbind(mz_list[[i]], matrix_list[[i]])
}

tictoc::tic("fuse")
rm(msi, mz, range_list, mz_list)
intensity_matrix <- do.call(cbind, matrix_list)
rm(matrix_list)
tictoc::toc()

tictoc::toc()
# Save data as *.csv
tictoc::tic("writing")
write.table(intensity_matrix,
            file = paste(file, ".csv"), # as previous filename
            row.names = TRUE,
            na = "",
            col.names = FALSE,
            sep = ",")
tictoc::toc()
tictoc::toc()