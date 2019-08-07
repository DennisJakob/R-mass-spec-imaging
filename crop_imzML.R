# edit imzML 
library(Cardinal)

ImzMLfile <- "tmp/20190124_MS37_MetaboliteMix_Spots_60_900_HCCA_pos_A28_100um_75x89"
binwidth <- 0.001
pixelrange <- 1692:4361 # to keep
outfile <- "tmp/20190510_MetaboliteMix_HCCAp_60_900_A28_100um_25x89"


# import old dataset
msidata <- readImzML(ImzMLfile,
                     attach.only = FALSE,
                     as="MSImagingExperiment",
                     resolution = binwidth, # binning width
                     units = "mz")
# crop dataset
msidata2 <- msidata[,pixelrange]
# write new dataset
writeImzML(msidata2,
           outfile)
# import new dataset
msidata2 <- readImzML(outfile,
                      attach.only = FALSE,
                      as="MSImagingExperiment",
                      resolution = binwidth, # binning width
                      units = "mz")