# edit imzML 

# crop ms image
msidata <- readImzML("datasets/20190124_MS37_MetaboliteMix_Spots_60_900_HCCA_pos_A28_100um_75x89",
                     attach.only = FALSE,
                     as="MSImagingExperiment",
                     resolution = 0.0005, # binning width
                     units = "mz")
msidata2 <- msidata[,1692:4361] #30x74
writeImzML(msidata2,
           "20190510_MetaboliteMix_HCCAp_60_900_A28_100um_25x89")
msidata2 <- readImzML("datasets/20190510_MetaboliteMix_HCCAp_60_900_A28_100um_25x89",
                      attach.only = FALSE,
                      as="MSImageSet",
                      resolution = 0.0005, # binning width
                      units = "mz")