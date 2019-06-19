# Load packages
library(Cardinal) # MSI analysis tool

# Import datafiles
msi_file <- "20190523_MS37_MetaboliteMix_Spots_60_900_SDHB_pos_A25_75um_88x70"
prediction_file <- "prediction_MetaMix.csv"

msi <- readImzML(msi_file,
                 as = "MSImagingExperiment",
                 resolution = 0.01,
                 units = "mz")
prediction <- read.csv(prediction_file,
                       header = TRUE,
                       sep = ",",
                       dec = ".",
                       row.names = NULL)

# Normalization
msi_TIC <- process(normalize(msi,
                             method = "tic"))

# Save files for later use
save(msi_TIC, "msi_TIC.Rdata")

# Plotting ion images 
options(Cardinal.dark = FALSE) # For dark mode set Cardinal.dark=TRUE
text_col <- "black" # For dark mode use "white"

for (i in 1:nrow(prediction)) {
  png(filename = paste("MetaMix_SDHBp_", # Experimental design for file name
                       prediction$compound[i],
                       ".png",
                       sep = ""),
      width = 2000,
      height = 2000,
      pointsize = 12,
      res = 200)
  par(mfrow = c(2,2),
      mar = c(5, 5, 5, 5),
      oma = c(0, 0, 5, 0),
      cex = 0.8)
  for (j in 6:9) {
    print(image(msi,
                mz = prediction[i, j],
                plusminus = 0.01, # Plotting deviation
                colorscale = magma,
                colorkey = TRUE,
                layout = FALSE))
    mtext(colnames(prediction[j]), # Minor title --> Adduct name
          col = text_col,
          side = 3,
          line = 1,
          adj = 0,
          cex = 1.5)
  }
  mtext(prediction$compound[i], # Main title --> Compound name
        col = text_col, 
        side = 3,
        line = 0,
        outer = TRUE,
        cex = 2)
  dev.off()
}
