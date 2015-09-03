.onLoad <- function(libname, pkgname) {
  # Wickham it's the man!
  op <- options()
  op.water <- list(
    water.overwrite = "TRUE",
    water.inmemory = "FALSE",
    water.SRTMrepo = "NULL",
    water.autoAoi = "TRUE"
  )
  toset <- !(names(op.water) %in% names(op))
  if(any(toset)) options(op.water[toset])
  
  invisible()
}