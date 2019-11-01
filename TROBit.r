Areas = c(0.1, 0.6, 0.3, 0.2, 0.7, 0.4)

detail = 100

bootstraps = 1000

randomise_treeLoc <- function(bootN) {
    print(bootN)
    cs = seq(0, 1, length.out = detail + 1)[-1] - 1/(2*detail)
    cs = cbind(rep(cs, each = detail), rep(cs, detail), 0)
    
    xs   = runif(length(Areas), 0, 1)
    ys   = runif(length(Areas), 0, 1)
   
    radius =  sqrt(Areas/pi)
   
    testTree <- function(x, y, r) {
        test = (x - cs[,1])^2 + (y - cs[,2])^2 <= r^2
        cs[test,3] = 1
        return(cs)
    }
    
    for (i in 1:length(Areas)) cs = testTree(xs[i], ys[i], radius[i])
    mean(cs[,3])
}

boots = sapply(1:bootstraps, randomise_treeLoc)