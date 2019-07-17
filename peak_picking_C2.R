## apply peakpick function
tic("total")
tic("peakpicking")
register(SerialParam())
peaks <- process(peakPick(msi),
                 method = "adaptive",
                 SNR = 10,
                 bpparam = SerialParam())
toc()
tic("transformation")
## creat list of all peaks in all pixels
peaks_list <- pData(mzData(imageData(peaks)))

## unlist peaks
peaks_unlist <- unlist(peaks_list)
head(peaks_unlist)

## keep only unique peaks
peaks_unique <- unique(peaks_unlist)

## sort by increasing m/z
peaks_unique_sorted <- sort(peaks_unique,
                            decreasing = FALSE)
toc()
tic("writing")
## save as text file
write.table(peaks_unique_sorted,
            file = "peaklist.txt",
            row.names = FALSE,
            col.names = FALSE,
            sep = ",")
toc()
toc()