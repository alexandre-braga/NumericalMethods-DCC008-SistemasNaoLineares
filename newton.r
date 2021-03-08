#!/usr/bin/env Rscript

f <- function(x1) {
    return (x1^2 + log(x1))
}

df <- function(x1) {
    return (2*x1 + 1/x1)
}

x1 <- 0.75
erro <- 1
tol <- 0.0001
k <- 0

print(sprintf("%s    %s    %s    %s    %s", "k", "x", "f(x)", "df(x)", "erro"))

while (erro > tol) {
    x2 <- x1 - f(x1)/df(x1)
    erro <- abs(x2-x1)
    k <- k+1
    print(sprintf("%d    %.4f    %.4f    %.4f    %s", k, x2, f(x2), df(x2), formatC(erro, format="e")))
    x1 <- x2
}
