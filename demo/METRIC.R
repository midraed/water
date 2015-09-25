library(water)

waterOptions(writeResults = F)


csvfile <- system.file("extdata", "apples.csv", package="water")
MTLfile <- system.file("extdata", "L7.MTL.txt", package="water")
WeatherStation <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y",
                  lat=-35.42222, long= -71.38639, elev=201, height= 2.2,
                  MTL = MTLfile)

print(WeatherStation)

plot(WeatherStation, alldata=FALSE)


aoi <- createAoi(topleft = c(273110, -3914450), 
                 bottomright = c( 288050, -3926650), EPSG = 32619)

image.DN <- L7_Talca[[c(1:5,7)]]

B6 <- L7_Talca[[6]]

plot(image.DN)

Energy.Balance <- METRIC.EB(image.DN = image.DN, plain=TRUE, 
                            WeatherStation = WeatherStation, 
                            MTL=MTLfile, sat="L7", thermalband=B6)


plot(Energy.Balance)

ET.24 <- ET24h(Rn = Energy.Balance$NetRadiation, G=Energy.Balance$SoilHeat,
               H = Energy.Balance$LatentHeat,
               Ts=Energy.Balance$surfaceTemperature, WeatherStation, 
               ETr.daily=6.57)
