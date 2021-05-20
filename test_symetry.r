x = seq(0.001, 0.999, 0.001)

logskew <- function(x, a, b) 
    log(x^a/((1-x)^b))
    
x = logskew(x, 1, 1)

plot(c(0, 1), c(-10, 10), type ='n', xlab = '', ylab ='')
lines(x, logskew(x, 1, 1))
lines(x, logskew(x,   5, 1)-logskew(0.5,   5, 1), col = 'red'  , lty = 2, lwd = 2)
lines(x, logskew(x, 1/5, 1)-logskew(0.5,  1/5, 1), col = 'red'  , lty = 3, lwd = 2)

lines(x, logskew(x, 1,   5)-logskew(0.5, 1,   5), col = 'blue' , lty = 2, lwd = 2)
lines(x, logskew(x, 1, 1/5)-logskew(0.5, 1, 1/5), col = 'blue' , lty = 3, lwd = 2)

lines(x, logskew(x, 5,   5), col = 'green', lty = 2, lwd = 2)
lines(x, logskew(x, 1/5,   1/5), col = 'green', lty = 3, lwd = 2)