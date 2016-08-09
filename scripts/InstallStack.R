
# Update stiPrEP stack

system("git pull")
devtools::install_github("statnet/EpiModel")
devtools::install_github("statnet/EpiModelHPC")
devtools::install_github("statnet/tergmLite", subdir = "tergmLite")
devtools::install_github("statnet/EpiModelHIV", ref = "port-sti")


# system("scp scripts/burnin/*.burn.[Rs]* hyak:/gscratch/csde/sjenness/sti")
# system("scp scripts/burnin/abc.parms.1pct.rda hyak:/gscratch/csde/sjenness/sti")
# system("scp source/*.* hyak:/gscratch/csde/sjenness/sti/source/")
