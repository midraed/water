#' Parameters for the selection of anchor pixels
#' @description
#' Parameters of calcAnchors() for METRIC.EB(anchors.control = ).
#' @param deltaTemp  temperature margin (K, default 5): cold candidates must be
#' colder than the 5\% quantile of Ts over well-vegetated pixels plus
#' deltaTemp, and hot candidates warmer than the maximum Ts over bare-soil
#' pixels minus deltaTemp
#' @param minDist    minimum distance (m) between two anchor pixels of the same
#' type
#' @param WSbuffer   maximum distance (m) to the weather station. If no pixel
#' is within WSbuffer, the whole image is used, with a warning
#' @param min.candidates  minimum number of hot and of cold candidates below
#' which a warning is given (the area may be too small or too homogeneous)
#' @param cold.NDVI      anchors.method = "percentile": cold candidates have
#' NDVI at or above this quantile of the NDVI of the area (default 0.95)
#' @param cold.LAI       "percentile": and LAI of at least cold.LAI (default 4)
#' @param cold.fraction  "percentile": the coldest fraction of the cold
#' candidates are the cold anchors (default 0.2)
#' @param hot.NDVI       "percentile": hot candidates have NDVI at or below
#' this quantile (default 0.1), above 0 (no water) and albedo 0.10 to 0.35
#' @param hot.fraction   "percentile": the hottest fraction of the hot
#' candidates are the hot anchors (default 0.2)
#' @param homogeneity    "percentile": maximum standard deviation of Ts (K) in
#' the 5 x 5 pixels around a candidate (default 1.5)
#' @return A list of class waterControl with the parameters.
#' @family METRIC model functions
#' @seealso \code{\link{calcAnchors}}
#' @examples
#' anchors.control(minDist = 1000)
#' @export
anchors.control <- function(deltaTemp = 5, minDist = 500, WSbuffer = 30000,
                            min.candidates = 20, cold.NDVI = 0.95, cold.LAI = 4,
                            cold.fraction = 0.2, hot.NDVI = 0.1, hot.fraction = 0.2,
                            homogeneity = 1.5){
  positive(deltaTemp, "deltaTemp"); positive(minDist, "minDist", zero = TRUE)
  positive(WSbuffer, "WSbuffer"); positive(min.candidates, "min.candidates", zero = TRUE)
  for(p in c("cold.NDVI", "cold.fraction", "hot.NDVI", "hot.fraction")){
    v <- get(p)
    if(!is.numeric(v) || length(v) != 1 || !isTRUE(v > 0 && v <= 1)){
      stop(p, " must be a number above 0 and up to 1", call. = FALSE)
    }
  }
  positive(cold.LAI, "cold.LAI", zero = TRUE); positive(homogeneity, "homogeneity")
  structure(list(deltaTemp = deltaTemp, minDist = minDist, WSbuffer = WSbuffer,
                 min.candidates = min.candidates, cold.NDVI = cold.NDVI,
                 cold.LAI = cold.LAI, cold.fraction = cold.fraction,
                 hot.NDVI = hot.NDVI, hot.fraction = hot.fraction,
                 homogeneity = homogeneity),
            class = c("waterControl", "list"))
}

#' Parameters for the iteration of H and the aerodynamic resistance
#' @description
#' Parameters of calcH() for METRIC.EB(calcH.control = ).
#' @param maxit     maximum number of iterations
#' @param u200.min  minimum wind speed at 200 m (m/s). With a lower wind
#' speed extrapolated from the station, u200.min is used, with a warning. 2
#' m/s corresponds to about 1 m/s at 2 m over grass. A choice of this package,
#' not of Allen et al. (2007)
#' @param damping   weight of the previous aerodynamic resistance in each
#' iteration, 0 (default, no damping) to below 1: r_ah = (1 - damping) * new +
#' damping * previous. Damping keeps the same solution and helps runs that
#' oscillate instead of converging. A choice of this package
#' @return A list of class waterControl with the parameters.
#' @family METRIC model functions
#' @seealso \code{\link{calcH}}
#' @examples
#' calcH.control(maxit = 30)
#' @export
calcH.control <- function(maxit = 20, u200.min = 2, damping = 0){
  if(!is.numeric(maxit) || length(maxit) != 1 || !isTRUE(maxit >= 1) || maxit %% 1 != 0){
    stop("maxit must be a whole number, 1 or more", call. = FALSE)
  }
  positive(u200.min, "u200.min")
  if(!is.numeric(damping) || length(damping) != 1 || !isTRUE(damping >= 0 && damping < 1)){
    stop("damping must be a number from 0 to below 1", call. = FALSE)
  }
  structure(list(maxit = maxit, u200.min = u200.min, damping = damping),
            class = c("waterControl", "list"))
}

#' @export
print.waterControl <- function(x, ...){
  utils::str(unclass(x), give.head = FALSE)
  invisible(x)
}

# a single positive number (or zero)
positive <- function(x, name, zero = FALSE){
  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) && (x > 0 || (zero && x == 0))
  if(!ok){
    stop(name, " must be a number above 0", if(zero) " (or 0)", call. = FALSE)
  }
  invisible(x)
}

# A control list: the defaults of the constructor called name updated with the
# elements of control (a list from it or a plain list); unknown names stop.
# The constructor is taken from the package namespace, as the arguments of
# METRIC.EB() have the same names
controlList <- function(control, name){
  fun <- get(name, envir = environment(controlList), mode = "function")
  if(is.null(control)){return(fun())}
  if(!is.list(control)){stop(name, " must be a list, e.g. from ", name, "()", call. = FALSE)}
  unknown <- setdiff(names(control), names(formals(fun)))
  if(length(unknown) || (length(control) && is.null(names(control)))){
    stop("Unknown elements in ", name, ": ", paste(unknown, collapse = ", "),
         ". Use ", paste(names(formals(fun)), collapse = ", "), call. = FALSE)
  }
  do.call(fun, unclass(control))
}
