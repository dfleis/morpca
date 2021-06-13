#' Prints 'Hello, world!'
#'
#' @export
hello <- function(id, x) {
  if (id == 1) {
    .Call("hello_world", x, PACKAGE = "morpca")
  } else if (id == 2) {
    hello_world2(x)
  } else if (id == 3) {
    hello_world(x)
  } else if (id == 4) {
    .Call("hello_world2", x, package = "morpca")
  }
  #hello_world()
}
