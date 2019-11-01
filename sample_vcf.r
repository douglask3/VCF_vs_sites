qunats = c(0.0, 0.01, 0.05, 0.1, 0.25)
Ygrid_size = 250 * 250


sample_pc <- function(pc) {
    Y = rep(0,Ygrid_size)

    nc = round(pc * Ygrid_size)
    prob = (1:Ygrid_size)^100
    index = sample(1:Ygrid_size, nc, FALSE, prob = )
    Y[index] = 1
    
     
    sampleY <- function() {
         xp = sample(1:(250*250 - 50 *100), 1)
         mean(Y[xp:(xp + 50 *100)])
    }

    xs = replicate(1000,sampleY())
    return(xs)
}

make.transparent <- function(col, transparency) {
     ## Semitransparent colours
     tmp <- col2rgb(col)/255
     rgb(tmp[1,], tmp[2,], tmp[3,], alpha=1-transparency)
}   

matrix2list <- function(x) split(x, rep(1:ncol(x), each = nrow(x)))

pcs = seq(0, 1, 0.05)
xs = sapply(pcs, sample_pc)

plot(pcs, apply(xs, 2, mean), type = 'l')

xl = apply(xs, 2, quantile, qunats)
xu = apply(xs, 2, quantile, 1 - rev(qunats))

addPolygon <- function(xl, xu, alpha) 
    polygon(c(pcs, rev(pcs)), c(xl, rev(xu)), col = make.transparent('black',  alpha), border = NA)
    
mapply(addPolygon, matrix2list(t(xl)), matrix2list(t(xu)), seq(0.99 , 0.95, length.out = length(qunats)))

