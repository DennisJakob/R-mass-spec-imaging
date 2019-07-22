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
pixel       <- as.character(seq(1, ncol(msi),1))
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
  matrix_list[[i]]            <- spectra(msi)[mz_list[[i]],]
  colnames(matrix_list[[i]])  <- pixel
  rownames(matrix_list[[i]])  <- mz_list[[i]]
  matrix_list[[i]]            <- matrix_list[[i]][rowSums(matrix_list[[i]] > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
  matrix_list[[i]]            <- t(matrix_list[[i]])
  mz_list[[i]]                <- as.numeric(colnames(matrix_list[[i]]))
  matrix_list[[i]]            <- rbind(mz_list[[i]], matrix_list[[i]])
}

for (i in 1:length(mz_list)) {
  matrix_list[[i]]             <- spectra(msi)[mz_list[[i]],]
  matrix_list[[i]]             <- t(matrix_list[[i]])
  colnames(matrix_list[[i]])   <- mz_list[[i]]
  matrix_list[[i]]             <- matrix_list[[i]][, colSums(matrix_list[[i]] > i_threshold) >= 1]
  mz_list[[i]]                 <- as.numeric(colnames(matrix_list[[i]]))
  rownames(matrix_list[[i]])   <- pixel
  matrix_list[[i]]             <- rbind(mz_list[[i]], matrix_list[[i]])
}


tictoc::tic("matrix1")
matrix1             <- spectra(msi)[range1,]
colnames(matrix1)   <- pixel
rownames(matrix1)   <- mz1
matrix1             <- matrix1[rowSums(matrix1 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix1             <- t(matrix1)
mz1                 <- as.numeric(colnames(matrix1))
matrix1             <- rbind(mz1, matrix1)
tictoc::toc()

tictoc::tic("matrix2")
matrix2             <- spectra(msi)[range2,]
colnames(matrix2)   <- pixel
rownames(matrix2)   <- mz2
matrix2             <- matrix2[rowSums(matrix2 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix2             <- t(matrix2)
mz2                 <- as.numeric(colnames(matrix2))
matrix2             <- rbind(mz2, matrix2)
tictoc::toc()

tictoc::tic("matrix3")
matrix3             <- spectra(msi)[range3,]
colnames(matrix3)   <- pixel
rownames(matrix3)   <- mz3
matrix3             <- matrix3[rowSums(matrix3 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix3             <- t(matrix3)
mz3                 <- as.numeric(colnames(matrix3))
matrix3             <- rbind(mz3, matrix3)
tictoc::toc()

tictoc::tic("matrix4")
matrix4             <- spectra(msi)[range4,]
colnames(matrix4)   <- pixel
rownames(matrix4)   <- mz4
matrix4             <- matrix4[rowSums(matrix4 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix4             <- t(matrix4)
mz4                 <- as.numeric(colnames(matrix4))
matrix4             <- rbind(mz4, matrix4)
tictoc::toc()

tictoc::tic("matrix5")
matrix5             <- spectra(msi)[range5,]
colnames(matrix5)   <- pixel
rownames(matrix5)   <- mz5
matrix5             <- matrix5[rowSums(matrix5 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix5             <- t(matrix5)
mz5                 <- as.numeric(colnames(matrix5))
matrix5             <- rbind(mz5, matrix5)
tictoc::toc()

tictoc::tic("matrix6")
matrix6             <- spectra(msi)[range6,]
colnames(matrix6)   <- pixel
rownames(matrix6)   <- mz6
matrix6             <- matrix6[rowSums(matrix6 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix6             <- t(matrix6)
mz6                 <- as.numeric(colnames(matrix6))
matrix6             <- rbind(mz6, matrix6)
tictoc::toc()

tictoc::tic("matrix7")
matrix7             <- spectra(msi)[range7,]
colnames(matrix7)   <- pixel
rownames(matrix7)   <- mz7
matrix7             <- matrix7[rowSums(matrix7 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix7             <- t(matrix7)
mz7                 <- as.numeric(colnames(matrix7))
matrix7             <- rbind(mz7, matrix7)
tictoc::toc()

tictoc::tic("matrix8")
matrix8             <- spectra(msi)[range8,]
colnames(matrix8)   <- pixel
rownames(matrix8)   <- mz8
matrix8             <- matrix8[rowSums(matrix8 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix8             <- t(matrix8)
mz8                 <- as.numeric(colnames(matrix8))
matrix8             <- rbind(mz8, matrix8)
tictoc::toc()

tictoc::tic("matrix9")
matrix9             <- spectra(msi)[range9,]
colnames(matrix9)   <- pixel
rownames(matrix9)   <- mz9
matrix9             <- matrix9[rowSums(matrix9 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix9             <- t(matrix9)
mz9                 <- as.numeric(colnames(matrix9))
matrix9             <- rbind(mz9, matrix9)
tictoc::toc()

tictoc::tic("matrix10")
matrix10             <- spectra(msi)[range10,]
colnames(matrix10)   <- pixel
rownames(matrix10)   <- mz10
matrix10             <- matrix10[rowSums(matrix10 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix10             <- t(matrix10)
mz10                 <- as.numeric(colnames(matrix10))
matrix10             <- rbind(mz10, matrix10)
tictoc::toc()

tictoc::tic("matrix11")
matrix11             <- spectra(msi)[range11,]
colnames(matrix11)   <- pixel
rownames(matrix11)   <- mz11
matrix11             <- matrix11[rowSums(matrix11 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix11             <- t(matrix11)
mz11                 <- as.numeric(colnames(matrix11))
matrix11             <- rbind(mz11, matrix11)
tictoc::toc()

tictoc::tic("matrix12")
matrix12             <- spectra(msi)[range12,]
colnames(matrix12)   <- pixel
rownames(matrix12)   <- mz12
matrix12             <- matrix12[rowSums(matrix12 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix12             <- t(matrix12)
mz12                 <- as.numeric(colnames(matrix12))
matrix12             <- rbind(mz12, matrix12)
tictoc::toc()

tictoc::tic("matrix13")
matrix13             <- spectra(msi)[range13,]
colnames(matrix13)   <- pixel
rownames(matrix13)   <- mz13
matrix13             <- matrix13[rowSums(matrix13 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix13             <- t(matrix13)
mz13                 <- as.numeric(colnames(matrix13))
matrix13             <- rbind(mz13, matrix13)
tictoc::toc()

tictoc::tic("matrix14")
matrix14             <- spectra(msi)[range14,]
colnames(matrix14)   <- pixel
rownames(matrix14)   <- mz14
matrix14             <- matrix14[rowSums(matrix14 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix14             <- t(matrix14)
mz14                 <- as.numeric(colnames(matrix14))
matrix14             <- rbind(mz14, matrix14)
tictoc::toc()

tictoc::tic("matrix15")
matrix15             <- spectra(msi)[range15,]
colnames(matrix15)   <- pixel
rownames(matrix15)   <- mz15
matrix15             <- matrix15[rowSums(matrix15 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix15             <- t(matrix15)
mz15                 <- as.numeric(colnames(matrix15))
matrix15             <- rbind(mz15, matrix15)
tictoc::toc()

tictoc::tic("matrix16")
matrix16             <- spectra(msi)[range16,]
colnames(matrix16)   <- pixel
rownames(matrix16)   <- mz16
matrix16             <- matrix16[rowSums(matrix16 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix16             <- t(matrix16)
mz16                 <- as.numeric(colnames(matrix16))
matrix16             <- rbind(mz16, matrix16)
tictoc::toc()

tictoc::tic("matrix17")
matrix17             <- spectra(msi)[range17,]
colnames(matrix17)   <- pixel
rownames(matrix17)   <- mz17
matrix17             <- matrix17[rowSums(matrix17 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix17             <- t(matrix17)
mz17                 <- as.numeric(colnames(matrix17))
matrix17             <- rbind(mz17, matrix17)
tictoc::toc()

tictoc::tic("matrix18")
matrix18             <- spectra(msi)[range18,]
colnames(matrix18)   <- pixel
rownames(matrix18)   <- mz18
matrix18             <- matrix18[rowSums(matrix18 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix18             <- t(matrix18)
mz18                 <- as.numeric(colnames(matrix18))
matrix18             <- rbind(mz18, matrix18)
tictoc::toc()

tictoc::tic("matrix19")
matrix19             <- spectra(msi)[range19,]
colnames(matrix19)   <- pixel
rownames(matrix19)   <- mz19
matrix19             <- matrix19[rowSums(matrix19 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix19             <- t(matrix19)
mz19                 <- as.numeric(colnames(matrix19))
matrix19             <- rbind(mz19, matrix19)
tictoc::toc()

tictoc::tic("matrix20")
matrix20             <- spectra(msi)[range20,]
colnames(matrix20)   <- pixel
rownames(matrix20)   <- mz20
matrix20             <- matrix20[rowSums(matrix20 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix20             <- t(matrix20)
mz20                 <- as.numeric(colnames(matrix20))
matrix20             <- rbind(mz20, matrix20)
tictoc::toc()

tictoc::tic("matrix21")
matrix21             <- spectra(msi)[range21,]
colnames(matrix21)   <- pixel
rownames(matrix21)   <- mz21
matrix21             <- matrix21[rowSums(matrix21 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix21             <- t(matrix21)
mz21                 <- as.numeric(colnames(matrix21))
matrix21             <- rbind(mz21, matrix21)
tictoc::toc()

tictoc::tic("matrix22")
matrix22             <- spectra(msi)[range22,]
colnames(matrix22)   <- pixel
rownames(matrix22)   <- mz22
matrix22             <- matrix22[rowSums(matrix22 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix22             <- t(matrix22)
mz22                 <- as.numeric(colnames(matrix22))
matrix22             <- rbind(mz22, matrix22)
tictoc::toc()

tictoc::tic("matrix23")
matrix23             <- spectra(msi)[range23,]
colnames(matrix23)   <- pixel
rownames(matrix23)   <- mz23
matrix23             <- matrix23[rowSums(matrix23 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix23             <- t(matrix23)
mz23                 <- as.numeric(colnames(matrix23))
matrix23             <- rbind(mz23, matrix23)
tictoc::toc()

tictoc::tic("matrix24")
matrix24             <- spectra(msi)[range24,]
colnames(matrix24)   <- pixel
rownames(matrix24)   <- mz24
matrix24             <- matrix24[rowSums(matrix24 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix24             <- t(matrix24)
mz24                 <- as.numeric(colnames(matrix24))
matrix24             <- rbind(mz24, matrix24)
tictoc::toc()

tictoc::tic("matrix25")
matrix25             <- spectra(msi)[range25,]
colnames(matrix25)   <- pixel
rownames(matrix25)   <- mz25
matrix25             <- matrix25[rowSums(matrix25 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix25             <- t(matrix25)
mz25                 <- as.numeric(colnames(matrix25))
matrix25             <- rbind(mz25, matrix25)
tictoc::toc()

tictoc::tic("matrix26")
matrix26             <- spectra(msi)[range26,]
colnames(matrix26)   <- pixel
rownames(matrix26)   <- mz26
matrix26             <- matrix26[rowSums(matrix26 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix26             <- t(matrix26)
mz26                 <- as.numeric(colnames(matrix26))
matrix26             <- rbind(mz26, matrix26)
tictoc::toc()

tictoc::tic("matrix27")
matrix27             <- spectra(msi)[range27,]
colnames(matrix27)   <- pixel
rownames(matrix27)   <- mz27
matrix27             <- matrix27[rowSums(matrix27 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix27             <- t(matrix27)
mz27                 <- as.numeric(colnames(matrix27))
matrix27             <- rbind(mz27, matrix27)
tictoc::toc()

tictoc::tic("matrix28")
matrix28             <- spectra(msi)[range28,]
colnames(matrix28)   <- pixel
rownames(matrix28)   <- mz28
matrix28             <- matrix28[rowSums(matrix28 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix28             <- t(matrix28)
mz28                 <- as.numeric(colnames(matrix28))
matrix28             <- rbind(mz28, matrix28)
tictoc::toc()

tictoc::tic("matrix29")
matrix29             <- spectra(msi)[range29,]
colnames(matrix29)   <- pixel
rownames(matrix29)   <- mz29
matrix29             <- matrix29[rowSums(matrix29 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix29             <- t(matrix29)
mz29                 <- as.numeric(colnames(matrix29))
matrix29             <- rbind(mz29, matrix29)
tictoc::toc()

tictoc::tic("matrix30")
matrix30             <- spectra(msi)[range30,]
colnames(matrix30)   <- pixel
rownames(matrix30)   <- mz30
matrix30             <- matrix30[rowSums(matrix30 > i_threshold) >= 1, ] # strip matrix by deleting peaks under threshold
matrix30             <- t(matrix30)
mz30                 <- as.numeric(colnames(matrix30))
matrix30             <- rbind(mz30, matrix30)
tictoc::toc()

tictoc::tic("fuse")
rm(msi, mz)
intensity_matrix    <- cbind(matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7, matrix8, matrix9, matrix10, matrix11, matrix12, matrix13, matrix14, matrix15, matrix16, matrix17, matrix18, matrix19, matrix20, matrix21, matrix22, matrix23, matrix24, matrix25, matrix26, matrix27, matrix28, matrix29, matrix30)
rm(matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7, matrix8, matrix9, matrix10, matrix11, matrix12, matrix13, matrix14, matrix15, matrix16, matrix17, matrix18, matrix19, matrix20, matrix21, matrix22, matrix23, matrix24, matrix25, matrix26, matrix27, matrix28, matrix29, matrix30)
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