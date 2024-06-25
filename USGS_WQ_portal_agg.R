
# Load the dataRetrieval package
library(dataRetrieval)
# Using defaults:
parameterCd <- "00076"
parameterINFO <- readNWISpCode(parameterCd)


dataPH <- readWQPdata(
  statecode = "US:55",
  characteristicName = "pH"
)