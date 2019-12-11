logit <- function(x) {
    x[x>0.999999999] = 0.999999999
    x[x<0.000000001] = 0.000000001
    log(x/(1-x))
}
logistic <- function(x) 1/(1+exp(-x))