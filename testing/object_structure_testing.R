library("magrittr")
kdensity(mtcars$mpg, kernel = "gamma", start = "normal", bw = "nrd0") -> object



update(object, x = rexp(100)) %>% plot
object[["kernel"]] = "gaussian"

(object$bw = "ucv")

object %>% plot

f = function(x, y) {
  match.call()
}

current = list(x = object$x,
               bw = object$bw_str,
               adjust = object$adjust,
               kernel = object$kernel_str,
               start = object$start_str,
               support = object$support,
               na.rm = object$na.rm,
               normalized = object$normalized)

passed = list(sko = "haha",
              start = "gumbel")

listmerge(current, passed, type = "template")
