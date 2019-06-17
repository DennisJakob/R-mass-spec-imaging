# Load packages
library(Cardinal) # MSI analysis tool
library(plot3D)   # for external colorkey of MS images
library(viridis)  # color palette

# import datafiles
msi_file <- "20190523_MS37_MetaboliteMix_Spots_60_900_SDHB_pos_A25_75um_88x70"
prediction_file <- "prediction_MetaMix.csv"

msi <- readImzML(msi_file,
                 attach.only = FALSE,
                 #as="MSImagingExperiment",
                 #centroided = FALSE,
                 #attach.only = FALSE,  
                 resolution = 0.005,
                 units = "mz")
prediction <- read.csv(prediction_file,
                       header = TRUE,
                       sep = ",",
                       dec = ".",
                       row.names = NULL)

# normalization
msi_TIC <- normalize(msi,
                     method = "tic",
                     tic = 100)

# save files
save(msi_TIC, "msi.RDATa") # check!!!

# plotting ion images
## white background
for (i in 1:nrow(prediction)) {
  png(filename = paste("MetaMix_SDHBp_", # experiment design for file name
                       prediction$compound[i],
                       ".png",
                       sep = ""),
      width = 2000,
      height = 2000,
      pointsize = 12,
      bg = "white",
      res = 200)
  par(mfrow = c(2,2),
      mar = c(5, 5, 5, 5),
      oma = c(0, 0, 5, 0),
      cex = 0.8)
  for (j in 6:9) {
    image(msi_TIC,
          mz = prediction[i, j],
          plusminus = 0.001, # plotting deviation
          col.regions = viridis(n = 100), # or inferno
          colorkey = FALSE) # TRUE for relative scale bar in ion image
    colkey(col = viridis(n = 100),
           clim = c(0, 100),
           clab = "rel. int.",
           col.clab = "black", # default is same as main title
           cex.clab = 1.1,
           line.clab = 1,
           adj.clab = 0,
           side = 4,
           add = TRUE,
           length = 1,
           width = 0.5)
    mtext(colnames(prediction[j]),
          col = "black",
          side = 3,
          line = 1,
          adj = 0,
          cex = 1.5)
  }
  mtext(prediction$compound[i],
        col = "black",
        side = 3,
        line = 0,
        outer = TRUE,
        cex = 2)
  dev.off()
}
## black background
for (i in 1:nrow(prediction)) {
  png(filename = paste("MetaMix_SDHBp_", # experiment design for file name
                       prediction$compound[i],
                       ".png",
                       sep = ""),
      width = 2000,
      height = 2000,
      pointsize = 12,
      bg = "black", # image scales are not visible
      res = 200)
  par(mfrow = c(2,2),
      mar = c(5, 5, 5, 5),
      oma = c(0, 0, 5, 0),
      cex = 0.8)
  for (j in 6:9) {
    image(msi_TIC,
          mz = prediction[i, j],
          plusminus = 0.001, # plotting deviation
          col.regions = viridis(n = 100), # or inferno
          colorkey = FALSE) # TRUE for relative intensity scale bar in ion image
    colkey(col = viridis(n = 100),
           clim = c(0, 100),
           clab = "rel. int.",
           col.clab = "white", # default is same as main title
           cex.clab = 1.1,
           line.clab = 1,
           adj.clab = 0,
           side = 4,
           add = TRUE,
           length = 1,
           width = 0.5)
    mtext(colnames(prediction[j]),
          col = "white",
          side = 3,
          line = 1,
          adj = 0,
          cex = 1.5)
  }
  mtext(prediction$compound[i],
        col = "white",
        side = 3,
        line = 0,
        outer = TRUE,
        cex = 2)
  dev.off()
}