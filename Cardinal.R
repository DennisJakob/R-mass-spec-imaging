# Load packages
library(Cardinal) # MSI analysis tool

# Import datafiles
msi_file <- "20190625_MS37_ArrayPrint_Test_100_1100_SDHB_pos_A15_50um_87x112"
prediction_file <- "predictions/prediction_MetaMix.csv"

msi <- readImzML(name = msi_file,
                 folder = "tmp",
                 as = "MSImagingExperiment",
                 resolution = 0.001, # binning width
                 units = "mz")

prediction <- read.csv(prediction_file,
                       header = TRUE,
                       sep = ",",
                       dec = ".",
                       row.names = NULL)
colnames(prediction) <- c("ID", "compound", "formula", "mass", "M+H", "M+Na", "M+K", "M+DHB", "M-H", "M+Cl")

# Normalization
register(SerialParam())
msi_TIC <- process(normalize(msi, method = "tic", bpparam = SerialParam()))
save(msi_TIC, file = paste(msi_file, "_.RData", sep = ""))

# Plotting ion images 
options(Cardinal.dark = FALSE) # For dark mode set Cardinal.dark=TRUE
text_col <- "black" # For dark mode use "white"

for (i in 1:nrow(prediction)) {
  png(filename = paste("MetaMix_array_SDHB_", # Experimental design for file name
                       prediction$compound[i],
                       ".png",
                       sep = ""),
      width = 2000,
      height = 2000,
      pointsize = 12,
      res = 200)
  par(mfrow = c(2, 2),
      mar = c(5, 5, 5, 5),
      oma = c(0, 0, 5, 0),
      cex = 0.8)
  for (j in 5:8) {
    print(image(msi,
                mz = prediction[i, j], # Plotting deviation
                plusminus = 0.005,
                colorscale = viridis,
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
