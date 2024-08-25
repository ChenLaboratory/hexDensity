#' CPP functions
#' 

#' @export
test1 = function(y, ys, n, a) {
  res = .Call(`expSmooth`,y,ys,n,a)
}

#' @export
test2 = function() {
  res = .Call(`expSmooth2`)
}