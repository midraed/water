.onLoad <- function(libname, pkgname) {
  # Wickham it's the man!
  op <- options()
  op.water <- list(
    water.overwrite = "TRUE",
    water.autoWrite = "TRUE",
    water.destfolder = ".",
    water.SRTMrepo = "NULL",
    water.autoAoi = "TRUE"
  )
  toset <- !(names(op.water) %in% names(op))
  if(any(toset)) options(op.water[toset])
  packageStartupMessage("This package writes function's results to working directory. You can change")
  packageStartupMessage("output folder, or completely disable this feature using waterOptions()")

  invisible()
}