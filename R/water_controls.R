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
#' @return A list of class waterControl with the parameters.
#' @family METRIC model functions
#' @seealso \code{\link{calcAnchors}}
#' @examples
#' anchors.control(minDist = 1000)
#' @export
anchors.control <- function(deltaTemp = 5, minDist = 500, WSbuffer = 30000){
  positive(deltaTemp, "deltaTemp"); positive(minDist, "minDist", zero = TRUE)
  positive(WSbuffer, "WSbuffer")
  structure(list(deltaTemp = deltaTemp, minDist = minDist, WSbuffer = WSbuffer),
            class = c("waterControl", "list"))
}

#' Parameters for the iteration of H and the aerodynamic resistance
#' @description
#' Parameters of calcH() for METRIC.EB(calcH.control = ).
#' @param maxit  maximum number of iterations
#' @return A list of class waterControl with the parameters.
#' @family METRIC model functions
#' @seealso \code{\link{calcH}}
#' @examples
#' calcH.control(maxit = 30)
#' @export
calcH.control <- function(maxit = 20){
  if(!is.numeric(maxit) || length(maxit) != 1 || !isTRUE(maxit >= 1) || maxit %% 1 != 0){
    stop("maxit must be a whole number, 1 or more", call. = FALSE)
  }
  structure(list(maxit = maxit), class = c("waterControl", "list"))
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
