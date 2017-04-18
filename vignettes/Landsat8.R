## ---- message=FALSE------------------------------------------------------
library(water)
aoi <- createAoi(topleft = c(500000, -3644000), bottomright = c(526000, -3660000))

## ---- warning=FALSE, fig.width = 7---------------------------------------
raw_data_folder <- system.file("extdata", package="water")
image <- loadImage(path=raw_data_folder, aoi=aoi, sat="L8")
image.SR <- loadImageSR(path=raw_data_folder, aoi=aoi)
plot(image)

## ---- warning=FALSE, fig.width = 7---------------------------------------
csvfile <- system.file("extdata", "INTA.csv", package="water")
MTLfile <- system.file("extdata", "LC82320832016040LGN00_MTL.txt", package="water")
WeatherStation <- read.WSdata(WSdata = csvfile, 
                              datetime.format =  "%Y/%m/%d %H:%M", 
                              columns = c("datetime", "temp",
                              "RH", "pp", "radiation", "wind"), 
                              lat=-33.00513, long= -68.86469, elev=927, height= 2,
                              MTL=MTLfile)
Energy.Balance <- METRIC.EB(image.DN = image, image.SR = image.SR,
                            plain=TRUE, aoi=aoi, n = 5, WeatherStation = WeatherStation, 
                            ETp.coef = 1.2, sat="L8", alb.coeff = "Olmedo", LST.method = "SW", 
                            LAI.method = "metric2010", Z.om.ws = 0.03, MTL = MTLfile)
plot(Energy.Balance)

