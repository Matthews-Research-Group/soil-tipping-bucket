# Load the BioCro framework and the testing package
library(testthat)
library(BioCro)

# Load this package and test it
library(BioCroWaterPotential)
test_check('BioCroWaterPotential')
