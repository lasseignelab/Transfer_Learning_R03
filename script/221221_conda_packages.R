#SR_TAU_CELL
.libPaths("/data/user/jfisher7/.conda/envs/SR_TAU_CELL/lib/R/library")

#make sure the environment is clean
rm(list=ls())

info<- sessionInfo()
print(info)


# Listing packages
packages<- installed.packages()[,c(1,3)]
print("Listing packages")
print(packages)
