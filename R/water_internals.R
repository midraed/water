get.sat <- function(path){
  band <- substr(list.files(path=path, pattern=paste0("^L[EC]\\d+\\w+\\d+_B2.TIF")), 0,3)
  if(length(band)==0){
    print(paste("ERROR: I expected something like landsat.band = LC82320832013319LGN00_BX.TIF in ", path))
    return()
  }
  if(band =="LC8"){return("L8")}
  if(band =="LE7"){return("L7")}
  if(band != "LE7" & band != "LC8"){
    print("Can't establish sat from band names")
    return()}
}


save.load.clean <- function(imagestack, stack.names=NULL, file, ...){
  if(missing(file)){file <- paste(deparse(substitute(raw.image)),".tif", sep="")}
  writeRaster(imagestack, filename = file, ...)
  stack <- stack(file)
  names(stack) <- stack.names
  removeTmpFiles(h=0)
  return(stack)
}