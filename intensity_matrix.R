# Load package
library("tictoc")
library("Cardinal")

tictoc::tic("total")
tictoc::tic("import")
# Import MSi dataset (ImzML and ibd!)
file <- "tmp/20190523_MS37_MetaboliteMix_Spots_60_900_SDHB_pos_A25_75um_88x70" # file name
msi <- readImzML(file,
                 attach.only=F,
                 as="MSImagingExperiment",
                 resolution=0.001, # 0.0005 results in a bin width of 0.001 Da
                 units="mz")

tictoc::toc()

pixel               <- as.character(seq(1,ncol(msi),1))
mz                  <- as.vector(mz(msi))
range1              <-       1:50000
range2              <-   50001:100000
range3              <-  100001:150000
range4              <-  150001:200000
range5              <-  200001:250000
range6              <-  250001:300000
range7              <-  300001:350000
range8              <-  350001:400000
range9              <-  400001:450000
range10             <-  450001:500000
range11             <-  500001:550000
range12             <-  550001:600000
range13             <-  600001:650000
range14             <-  650001:700000
range15             <-  700001:750000
range16             <-  750001:800000
range17             <-  800001:850000
range18             <-  850001:900000
range19             <-  900001:950000
range20             <-  950001:1000000
range21             <- 1000001:1050000
range22             <- 1050001:1100000
range23             <- 1100001:1150000
range24             <- 1150001:1200000
range25             <- 1200001:1250000
range26             <- 1250001:1300000
range27             <- 1300001:1350000
range28             <- 1350001:1400000
range29             <- 1400001:1450000
range30             <- 1450001:1500001
mz1                <- mz[range1]
mz2                <- mz[range2]
mz3                <- mz[range3]
mz4                <- mz[range4]
mz5                <- mz[range5]
mz6                <- mz[range6]
mz7                <- mz[range7]
mz8                <- mz[range8]
mz9                <- mz[range9]
mz10               <- mz[range10]
mz11               <- mz[range11]
mz12               <- mz[range12]
mz13               <- mz[range13]
mz14               <- mz[range14]
mz15               <- mz[range15]
mz16               <- mz[range16]
mz17               <- mz[range17]
mz18               <- mz[range18]
mz19               <- mz[range19]
mz20               <- mz[range20]
mz21               <- mz[range21]
mz22               <- mz[range22]
mz23               <- mz[range23]
mz24               <- mz[range24]
mz25               <- mz[range25]
mz26               <- mz[range26]
mz27               <- mz[range27]
mz28               <- mz[range28]
mz29               <- mz[range29]
mz30               <- mz[range30]

i_threshold        <- 100

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
