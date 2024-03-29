# Get input filename prefix
ImzMLfile <- "tmp/HR2MSI_mouse_urinary_bladder_S096"
binwidth <- 0.001
triplet <- "tmp/13_Bladder_DHB"
matrix <- "CHCA"

myrows <- paste(c(triplet, "_rows.list"), collapse="")
mycols <- paste(c(triplet, "_cols.list"), collapse="")
myvals <- paste(c(triplet, "_vals.list"), collapse="")
mypeaks <- paste(c(triplet, "_peaknames.list"), collapse="")
myspots <- paste(c(triplet, "_spotnames.list"), collapse="")

# load packages
library(Cardinal) 
library(mass2adduct)

# import data
SDHBp_spots <- readImzML(ImzMLfile,
                         attach.only = FALSE,
                         as="MSImagingExperiment",
                         resolution = binwidth, # dependent on instrument setup, usually 0.001
                         units = "mz") # or ppm
d <- msimat(rows = myrows,
            cols = mycols,
            vals = myvals,
            peaks = mypeaks,
            spots = myspots)
#-------------------------------------------------------------------------

# mass2adduct analysis
## filter for top 10000 peaks
d.filter <- filterPeaks(d,"topX",10000)
## calculate mass differences
d.diff <- massdiff(d.filter)
## match mass differences to adducts
d.diff.annot <- adductMatch(d.diff,add=adducts2)
## calculate correlation
d.diff.annot.cor <- corrPairsMSI(d,d.diff.annot)
## delete not significant correlations
d.diff.annot.cor <- subset(d.diff.annot.cor,Significance==1)
## delete all matches for different matrices
d.diff.annot.cor <- subset(d.diff.annot.cor,matches!=matrix) 
## order the list based on the mass of the parent ion
d.diff.annot.cor <- d.diff.annot.cor[with(d.diff.annot.cor,order(A)),]
summary(d.diff.annot.cor)

#-------------------------------------------------------------------------

# imaging of correlations
png(filename = "MetaMix_SDHBp_corr_Carnitine.png",
    width = 2000,
    height = 2000,
    pointsize = 12,
    bg = "white",
    res = 200)
par(mfrow = c(2,2))
image(SDHBp_spots,
      mz = 162.112,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+H",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
image(SDHBp_spots,
      mz = 184.094,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+Na",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
image(SDHBp_spots,
      mz = 200.068,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+K",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
image(SDHBp_spots,
      mz = 298.128,
      plusminus = 0.001,
      col.regions = viridis(n = 100),
      colorkey = TRUE)
mtext("M+DHB",
      col = "black",
      side = 3,
      line = 1,
      adj = 0,
      cex = 1.5)
dev.off()