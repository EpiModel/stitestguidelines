
# Update stiPrEP stack

# Install EpiModel Stack
install.packages("EpiModel", dependencies = TRUE)

# Install latest Dev Versions of Packages
install.packages("remotes")

remotes::install_github(c("statnet/network",
                          "statnet/statnet.common",
                          "statnet/ergm",
                          "statnet/tergm",
                          "statnet/EpiModel",
                          "statnet/EpiModelHPC",
                          "statnet/tergmLite",
                          "EpiModel/EpiABC"))

remotes::install_github("EpiModel/EpiModelHIV-p", ref = "syph_ept_recalib")

