# Matrix Peak Prediction

# load packages
library(KEGGREST) # KEGG API
library(RcppClassic) # calculation of masses
library(Rdisop) # calculation of masses
library(ggplot2) # plotting
library(purrr) # functional programming toolkit

# create list
n_matrix <- c(1,1,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8)
n_H2O    <- c(0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7,8)
compound <- c("1 DHAP","1 DHAP-H2O",
              "2 DHAP","2 DHAP-H2O","2 DHAP-2H2O",
              "3 DHAP","3 DHAP-H2O","3 DHAP-2H2O","3 DHAP-3H2O",
              "4 DHAP","4 DHAP-H2O","4 DHAP-2H2O","4 DHAP-3H2O","4 DHAP-4H2O",
              "5 DHAP","5 DHAP-H2O","5 DHAP-2H2O","5 DHAP-3H2O","5 DHAP-4H2O","5 DHAP-5H2O",
              "6 DHAP","6 DHAP-H2O","6 DHAP-2H2O","6 DHAP-3H2O","6 DHAP-4H2O","6 DHAP-5H2O","6 DHAP-6H2O",
              "7 DHAP","7 DHAP-H2O","7 DHAP-2H2O","7 DHAP-3H2O","7 DHAP-4H2O","7 DHAP-5H2O","7 DHAP-6H2O","7 DHAP-7H2O",
              "8 DHAP","8 DHAP-H2O","8 DHAP-2H2O","8 DHAP-3H2O","8 DHAP-4H2O","8 DHAP-5H2O","8 DHAP-6H2O","8 DHAP-7H2O","8 DHAP-8H2O")

matrix_peaks <- data.frame(compound, n_matrix, n_H2O, row.names = NULL)

matrix_peaks$"mass"  <- matrix_peaks$n_matrix*getMolecule("C8H8O3")$exactmass-matrix_peaks$n_H2O*getMolecule("H2O")$exactmass
matrix_peaks$"M+H"   <- matrix_peaks$n_matrix*getMolecule("C8H8O3")$exactmass-matrix_peaks$n_H2O*getMolecule("H2O")$exactmass+getMolecule("H")$exactmass-0.00055
matrix_peaks$"M+Na"  <- matrix_peaks$n_matrix*getMolecule("C8H8O3")$exactmass-matrix_peaks$n_H2O*getMolecule("H2O")$exactmass+getMolecule("Na")$exactmass-0.00055
matrix_peaks$"M+K"   <- matrix_peaks$n_matrix*getMolecule("C8H8O3")$exactmass-matrix_peaks$n_H2O*getMolecule("H2O")$exactmass+getMolecule("K")$exactmass-0.00055
matrix_peaks$"M+NH4" <- matrix_peaks$n_matrix*getMolecule("C8H8O3")$exactmass-matrix_peaks$n_H2O*getMolecule("H2O")$exactmass+getMolecule("NH4")$exactmass-0.00055

# save list
write.csv(matrix_peaks,
          file = "DHAP_peaks.csv",
          row.names = FALSE)





