library("MALDIquant")

setwd("~/ownCloud/MPI/R workbench")
s <- read.csv("spectrum.csv",
              header = TRUE,
              sep = ";",
              dec = ".")[,1:2]
s <- createMassSpectrum(mass = s[,1],
                        intensity = s[,2],
                        metaData = list(name = "spectrum1"))
mass(s)
intensity(s)
metaData(s)




# Data Import
data("fiedler2009subset")
length(fiedler2009subset)

# Quality Control
any(sapply(fiedler2009subset, # missing data in the spectra
           isEmpty)) 
table(sapply(fiedler2009subset, # same number of data points
             length)) 
all(sapply(fiedler2009subset, # mass difference is equal
           isRegular)) 

plot(fiedler2009subset[[1]])
plot(fiedler2009subset[[16]])

# Variance Stabilization - using square root transformation
spectra <- transformIntensity(fiedler2009subset,
                              method = "sqrt")

# Smoothing - Savitzky-Golay-Filter
spectra <- smoothIntensity(spectra,
                            method = "SavitzkyGolay",
                            halfWindowSize = 10)

# Baseline Correction
## visualisation
baseline <- estimateBaseline(spectra[[16]],
                             method = "SNIP",
                             iterations = 100)
plot(spectra[[16]])
lines(baseline,
      col = "red",
      lwd = 2)

## remove baseline if it fits
spectra <- removeBaseline(spectra,
                          method = "SNIP",
                          iterations = 100)
plot(spectra[[1]])
plot(spectra[[16]])

# Intensity Calibration / Normalization - Total-Ion-Current-Calibration
spectra <- calibrateIntensity(spectra,
                              method = "TIC")

# Warping / Alignment - peak based
spectra <- alignSpectra(spectra, # want to investigate the impact of different parameters please use 'determineWarpingFunctions' instead of the easier 'alignSpectra'
                        halfWindowSize = 20,
                        SNR = 2,
                        tolerance = 0.002,
                        warpingMethod = "lowess")
samples <- factor(sapply(spectra, # look for replicates in metadata
                         function(x)metaData(x)$sampleName))
avgSpectra <- averageMassSpectra(spectra, # create mean of replicates
                                 labels = samples,
                                 method = "mean")

# Peak Detection
noise <- estimateNoise(avgSpectra[[1]]) #estimate noise
plot(avgSpectra[[1]],
     xlim = c(4000,
              5000),
     ylim = c(0,
              0.002))
lines(noise,
      col = "red")
lines(noise[, 1],
      noise[, 2]*2,
      col = "blue")

peaks <- detectPeaks(avgSpectra,
                     method = "MAD",
                     halfWindowSize = 20,
                     SNR = 2) # choose signal-noise-ratio of 2
plot(avgSpectra[[1]],
     xlim = c(4000,
              5000),
     ylim = c(0,
              0.002))
points(peaks[[1]],
       col = "red",
       pch = 4)

# Peak Binning
peaks <- binPeaks(peaks, # makes similar peak values identical
                  tolerance = 0.002)

# Feature Matrix
peaks <- filterPeaks(peaks, # remove less frequent peaks to remove false positive peaks
                     minFrequency = 0.25)
featureMatrix <- intensityMatrix(peaks, avgSpectra) # missing values are interpolated from corresponding spectrum
head(featureMatrix[,1:3])


