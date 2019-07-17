library("MALDIquant")

# Import Data

spectra <- import(getPathSpecies(), verbose=FALSE)

# Quality Control
any(sapply(spectra, # should be FALSE
           isEmpty)) 
table(sapply(spectra,
             length))
all(sapply(spectra, # should be TRUE
           isRegular))
spectra <- trim(spectra) # ensure all spectra have same mass range
idx <- sample(length(spectra), # randomly samples from spectra
              size = 2)
plot(spectra[[idx[1]]]) 
plot(spectra[[idx[2]]])

# Transforming and Smoothing
spectra <- transformIntensity(spectra, #square-root-transformation
                              method = "sqrt")
plot(spectra[[1]],
     type = "b", # to show lines and datapoints
     xlim = c(2235.3,
              2252.0),
     ylim = c(45,
              100))
abline(h = 72, # middle of ylim --> 10 data points over the dashed blue line --> halfWindowSize = 10
       col = 4,
       lty = 2)
plot(spectra[[1]],
     type = "b",
     xlim = c(11220,
              11250),
     ylim = c(24,
              40))
abline(h = 32,
       col = 4,
       lty = 2)
spectra <- smoothIntensity(spectra, # application of a 21 point (2*halfWindowSize+1) Savitzky-Golay-Filter to smooth the spectra
                           method = "SavitzkyGolay",
                           halfWindowSize = 10)

# Baseline Correction
iterations <- seq(from = 25, # define iterations
                   to = 100,
                   by = 25)
col <- rainbow(length(iterations)) # define different colors for each iteration
plot(spectra[[1]], 
     xlim = c(2000,
              12000))
for (i in seq(along = iterations)) { # draw different baseline estimates -> 25 iterations are already very flexible, 50 is not flexible enough --> iterations = 25 for baseline removal
  baseline <- estimateBaseline(spectra[[1]],
                               method = "SNIP",
                               iterations = iterations[i])
  lines(baseline,
        col = col[i],
        lwd = 2)
}
legend("topright",
       legend = iterations,
       col = col,
       lwd = 1)
spectra <- removeBaseline(spectra,
                          method = "SNIP",
                          iterations = 25)
plot(spectra[[1]])

# Intesity Calibration (TIC)
spectra <- calibrateIntensity(spectra, # equalize intensities across spectra
                              method = "TIC") # application of RMS method???

# Alignment
spectra <- alignSpectra(spectra)
metaData(spectra[[1]])$spot
spots <- sapply(spectra,
                function(x)metaData(x)$spot)
species <- sapply(spectra,
                  function(x)metaData(x)$sampleName)
head(spots)
head(species)
avgSpectra <- averageMassSpectra(spectra,
                                 labels = paste0(species,
                                                 spots))

# Peak Detection
snrs <- seq(from = 1, #define snrs steps
            to = 2.5,
            by = 0.5)
col <- rainbow(length(snrs)) #define color for each step
noise <- estimateNoise(avgSpectra[[1]], #estimates noise
                       method = "SuperSmoother")
plot(avgSpectra[[1]],
     xlim = c(6000,
              16000),
     ylim = c(0,
              0.0016))
for (i in seq(along = snrs)) { # SNR of 2 looks like a good compromise between sensitivity and specitivity
  lines(noise[, "mass"],
        noise[, "intensity"]*snrs[i],
        col = col[i],
        lwd = 2)
}
legend("topright",
       legend = snrs,
       col = col,
       lwd = 1)
peaks <- detectPeaks(avgSpectra,
                     SNR = 2,
                     halfWindowSize = 10)
plot(avgSpectra[[1]],
     xlim = c(6000,
              16000),
     ylim = c(0,
              0.0016))
points(peaks[[1]],
       col = "red",
       pch = 4)

# PostProcessing
peaks <- binPeaks(peaks)
peaks <- filterPeaks(peaks, minFrequency = 0.25) # low signal-noise-ratio to keep as much features as possible
spots <- sapply(avgSpectra, #recollect information
                function(x)metaData(x)$spot)
species <- sapply(avgSpectra,
                  function(x)metaData(x)$sampleName)
species <- factor(species) # convert to factor for later crossval
featureMatrix <- intensityMatrix(peaks,
                                 avgSpectra)
rownames(featureMatrix) <- paste(species,
                                 spots,
                                 sep = ".")

# Clustering
library("pvclust")
pv <- pvclust(t(featureMatrix),
              method.hclust = "ward.D2",
              method.dist = "euclidean")
plot(pv,
     print.num = FALSE)

# Diagonal Discriminant Analysis
library("sda")
ddar <- sda.ranking(Xtrain = featureMatrix,
                    L = species,
                    fdr = FALSE,
                    diagonal = TRUE)
plot(ddar)

# Linear Discriminant Analysis
ldar <- sda.ranking(Xtrain = featureMatrix,
                    L = species,
                    fdr = FALSE,
                    diagonal = FALSE)
plot(ldar)

# Variable Selection using Cross-Validation
library("crossval")
predfun <- function(Xtrain,
                    Ytrain,
                    Xtest,
                    Ytest,
                    numVars,
                    diagonal = FALSE) {
  ra <- sda.ranking(Xtrain, # estimate ranking and determine the best numVars variables
                    Ytrain,
                    verbose = FALSE,
                    diagonal = diagonal,
                    fdr = FALSE)
  selVars <- ra[,
                "idx"][1:numVars]
  sda.out <- sda(Xtrain[, # fit and predict
                        selVars,
                        drop = FALSE],
                 Ytrain,
                 diagonal = diagonal,
                 verbose = FALSE)
  ynew <- predict(sda.out,
                  Xtest[,
                        selVars,
                        drop = FALSE],
                  verbose = FALSE)$class
  acc <- mean(Ytest == ynew) # compute accuracy
  return(acc)
}
K <- 5 # number of folds
B <- 20 # number of repetitions

set.seed(12345)
cv.dda10 <- crossval(predfun,
                     X = featureMatrix,
                     Y = species,
                     K = K,
                     B = B,
                     numVars = 10,
                     diagonal = FALSE,
                     verbose = FALSE)
cv.dda10$stat

npeaks <- c(1:15, # number of peaks
            ncol(featureMatrix))
set.seed(12345)
cvsim.dda <- sapply(npeaks,function(i) { # estimate accuracy for DDA
  cv <- crossval(predfun,
                 X = featureMatrix,
                 Y = species,
                 K = K,
                 B = B,
                 numVars = i,
                 diagonal = TRUE,
                 verbose = FALSE)
  return(cv$stat)
})
set.seed(12345)
cvsim.lda <- sapply(npeaks,function(i) { # estimate accuracy for LDA
                    cv <- crossval(predfun,
                                   X = featureMatrix,
                                   Y = species,
                                   K = K,
                                   B = B,
                                   numVars = i,
                                   diagonal = FALSE,
                                   verbose = FALSE)
                    return(cv$stat)
                    })
result.sim <- cbind(nPeaks = npeaks, # LDA and DDA perform very similar
                    "DDA-ACC" = cvsim.dda,
                    "LDA-ACC" = cvsim.lda)
result.sim