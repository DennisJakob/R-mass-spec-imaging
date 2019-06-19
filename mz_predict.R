# load packages
library(KEGGREST) # KEGG API
library(RcppClassic) # calculation of masses
library(Rdisop) # calculation of masses
library(ggplot2) # plotting
library(purrr) # functional programming toolkit

# import datasets
## from self created list
MetaMix <- read.csv("MetaMix.csv",
                    header = TRUE,
                    sep = ",",
                    dec = ".")
## from PathwayTools
PathwayTools <- read.csv("New_SmartTable_MOX.csv",
                          header = TRUE,
                          sep = ";",
                          dec = ".",
                          stringsAsFactors = TRUE)
PathwayTools <- subset(PathwayTools,
                        PathwayTools$KEGG != "")
PathwayTools <- subset(PathwayTools,
                        PathwayTools$Chemical.Formula != "")
PathwayTools <- PathwayTools[, c("KEGG",
                                 "Name")]
colnames(PathwayTools) <- c("ID",
                            "compound")
rownames(PathwayTools) = NULL
### if compounds cannot be found: HTTP 404 -> manual subsetting
PathwayTools <- rbind(PathwayTools[1:203, ],
                      PathwayTools[205:417, ],
                      PathwayTools[419:531, ],
                      PathwayTools[533:534, ],
                      PathwayTools[536:547, ])

# create datasets
## manually
ID <-  c("C03557", "C00318", "C00185", "C00158", "C00475", "C04022", "C02679", "C00504", "C00031", "C00092", "C00037", "C00123", "C01384", "C00392", "C00140", "C01061", "C00079", "C00093", "C05682", "C00022", "C00121", "C00178", "C00086")
compound <- c("Aminomethylphosphonate", "Carnitine", "Cellobiose", "Citric acid" , "Cytidine", "Dimethylsulfoniopropionic acid", "Dodecanoic acid", "Folic acid", "Glucose", "Glucose-6-phosphate", "Glycine", "Leucine", "Maleic acid", "Mannitol", "N-Acetylglucosamine", "Nonanoic acid", "Phenylalanine", "Phosphoglyceric acid", "Phosphonoacetic acid", "Pyruvic acid", "Ribose", "Thymine", "Urea")
MetaMix <- data.frame(ID, compound, row.names = NULL)
MetaMix$ID <- as.character(MetaMix$ID)
MetaMix$compound <- as.character(MetaMix$compound)
## from KEGG metabolic pathway
listDatabases()
pathway                  <- keggList("pathway")
map                      <- keggGet("path:ko00680")
keggID                   <- as.data.frame(attributes(map[[1]]$COMPOUND))
compounds                <- as.data.frame(map[[1]]$COMPOUND)
KEGG_compounds           <- as.data.frame(c(keggID[1],
                                            compounds[1]))
colnames(KEGG_compounds) <- c("ID",
                              "compound")

# query information from KEGG
query_MetaMix <- list()
for (i in 1:nrow(MetaMix)){
  query_MetaMix[i] <- keggGet(MetaMix[i,
                                      "ID"])
}

# add information to list
for (i in 1:nrow(MetaMix)) {
  MetaMix$formula[i] <- query_MetaMix[[i]]$FORMULA
  MetaMix$mass[i]    <- getMolecule(MetaMix$formula[i])$exactmass
} # --> check chemical formular; no ions or polymers allowed
colnames(MetaMix) <- c("ID",
                       "compound",
                       "formula",
                       "mass")

# calculate adduct masses
MetaMix$"M+H"    <- MetaMix$mass+getMolecule("H")$exactmass-0.00055
MetaMix$"M+Na"   <- MetaMix$mass+getMolecule("Na")$exactmass-0.00055
MetaMix$"M+K"    <- MetaMix$mass+getMolecule("K")$exactmass-0.00055
MetaMix$"M+DHB"  <- MetaMix$mass+(getMolecule("C7H6O4")$exactmass-getMolecule("H2O")$exactmass)+getMolecule("H")$exactmass-0.00055
MetaMix$"M-H"    <- MetaMix$mass-getMolecule("H")$exactmass+0.00055
MetaMix$"M+Cl"   <- MetaMix$mass+getMolecule("Cl")$exactmass+0.00055

MetaMix[,5:ncol(MetaMix)] <- round(MetaMix[,5:ncol(MetaMix)], digit = 4) # round masses

# save list
write.csv(prediction_PolyGlycanMix,
          file = "prediction_PolyGlycanMix.csv",
          row.names = FALSE)

# histogram
X <- combined_tissue[combined_tissue$mass < 1000, ]
ggplot(X,
       aes(x = mass)) +
  geom_histogram(binwidth = 2) +
  xlab('Exact Masses of Predicted Compounds') +
  ylab('Frequency') +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=15)) +
  theme(axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=15)) 
ggsave('KEGG_pt_FreqChart<1000.png')