### ===========================================================================
### UTILITY FUNCTIONS
###
### Contains functions that might be or are useful in this project, but have
### nothing in particular to do with the theme of non-parametric density
### estimation. Will be better documented and test soon, hopefully...
### ===========================================================================

#' Recycles arguments.
#'
#' @param ... A list of arguments to be recycled.
#' @param prototype an optional argument. If given, repeats all arguments
#' up to the length of the prototype. If an element of the list has the name,
#' it is used. If not, the variable itself is used.
#' @details Recycles arguments so that all vectors are equally long. If a
#' prototype is given, each vector will have the same size as the prototype.
#' @examples
#'     a = 1:3
#'     b = letters[2:9]
#'     c = 9:20
#'
#'     # Returns a list where each element has length 5.
#'     recycle(a, b, c, prototype = 5)
#'
#'     # Each element has the same length as c.
#'     recycle(a, b, c, prototype = "c")
#'     recycle(a, b, c, prototype = c)

recycle = function(..., prototype) {

  dots = list(...)
  arg_names = utils::tail(sapply(as.list(match.call()),
                                 function(x) as.character(x)), -1)

  if(missing(prototype)) {
    names(dots) = names(arg_names)
    max_length = max(sapply(dots, length))
  } else {
    arg_names = utils::head(arg_names, -1)
    names(dots) = arg_names

    ## The rules work as follows: If it is a name, check the supplied list
    ## first. If itsn't there, use ordinary scoping to check. If it is a
    ## non-negative number, use it. If it's none of these, throw an error.

    subst_proto = deparse(substitute(prototype))
    if(subst_proto %in% arg_names) {
      max_length = length(dots[[subst_proto]])
    } else if (length(prototype) > 1) {
        max_length = length(prototype)
    } else if (is.numeric(prototype)){
        if(prototype >= 0) max_length = ceiling(prototype)
        else stop("supply a valid type. prototype is numeric and negative.")
    } else if (is.character(prototype)) {
        max_length = length(dots[[prototype]])
    } else {
        stop("supply a valid type for the prototype argument.")
    }
  }

  names(dots) = arg_names
  lapply(dots, rep, length.out = max_length, USE.NAMES = TRUE)

}

#' Puts default arguments into ellipses: ...
#'
#'
#' @param x a list of default arguments.
#' @param y a list of supplied arguments
#' @param type should x and y be merged (with y having priority),
#' or the elements x be a template filled with values from y?
#' @return a merged list where conflicts are solved in favour
#' of supplied.

listmerge = function(x, y, type = c("merge", "template")) {

  type = match.arg(type)

  if(length(y) == 0) {
    return(x)
  }

  ## Keep and not-keep are quite different.
  if(type == "merge") {
    matches = match(names(y), names(x))
    elements_to_discard = matches[!is.na(matches)]
    if(length(elements_to_discard) == 0)  {
      combined = c(y, x)
    } else{
      combined = c(y, x[-elements_to_discard])
    }
    return(combined)
  }

  if(type == "template") {
    matches = match(names(x), names(y))
    x[!is.na(matches)] = y[matches[!is.na(matches)]]
    return(x)
  }

}
