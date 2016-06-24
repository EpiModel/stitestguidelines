
# Update stiPrEP stack

system("git pull")
devtools::install_github("statnet/EpiModel")
devtools::install_github("statnet/EpiModelHPC")
devtools::install_github("statnet/tergmLite", subdir = "tergmLite")
devtools::install_github("statnet/EpiModelHIV")
