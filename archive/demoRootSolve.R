library(rootSolve)

fun <- function (x) cos(2*x)^3
curve(fun(x), 0, 8)
abline(h = 0, lty = 3)
uni <- uniroot(fun, c(0, 8))$root
points(uni, 0, pch = 16, cex = 2)

All <- uniroot.all(fun, c(0, 8))
points(All, y = rep(0, length(All)), pch = 16, cex = 2)

fun2 <- function(x) (x-3)^3
curve(fun2(x), 0, 8)
abline(h = 0, lty = 3)
uni <- uniroot(fun2, c(0, 8))$root
points(uni, 0, pch = 16, cex = 2)


# p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
# fun.1 <- function(x) x^2 + x
# fun.2 <- function(x) -1 * x + 10
# fun.3 <- function(x) 3 * x + 2
#
# p + layer(geom = "path",        # Default. Can be omitted.
#           stat = "function",
#           fun = fun.1,          # Give function
#           mapping = aes(color = "fun.1") # Give a meaningful name to color
# ) +
#   scale_x_continuous(limits = c(-5,5)) +
#   scale_color_manual(name = "Function", values = c("blue"))
