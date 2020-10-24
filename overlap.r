source("libs/make.transparent.r")

ntrees = 50
z = 1-runif(ntrees, 0, 1)^5
newPlot <- function(p, name) {
    plot(0.5, 0.5, cex =1000, pch = 19, col = "brown",
         xlim = c(0, 1), ylim = c(0, 1),
         axes = FALSE, xlab = '', ylab = '')
    
    x =  sample(seq(0, 1, 0.01), replace = TRUE, ntrees, prob = (1:101)^p)
    y =  sample(seq(0, 1, 0.01), replace = TRUE, ntrees)#, prob = (1:101)^p)
    #y = runif(ntrees, 0, 1)
    
    points(x, y, cex = z*10+3, col = make.transparent("green", 0.7), pch = 19)
    mtext(side = 3, line = -2, name, lwd = 1, cex = 1.5, adj = 0.1)
}

#png("figs/overlap.png", width = 25, height = 4, res = 300, units = 'in')
#par(mfrow = c(1,4), mar = rep(0.5, 4))

#sapply(c(0, 1, 2, 3), newPlot)
#dev.off()

png("figs/overlap.png", width = 7.2, height = 7.2*1.5/2, res = 300, units = 'in')
layout(rbind(1:2, 3:4), heights = c(1, 0.5))
par(mar = rep(0.5, 4), oma = c(0, 0.5, 0, 0))
mapply(newPlot, c(0, 2), c("Unenforeced overlap", "Enforced overlap"))


plotLine <- function(y) {
    plot(xaxt = 'n', yaxt = 'n', y, type = 'l', axes = FALSE, lwd = 2, lty = 2)
    lapply(1:2, axis,at = c(-9E9, 9E9))
    u <- par("usr") 
    arrows(u[1], u[3], u[2], u[3], code = 2, xpd = TRUE) 
    arrows(u[1], u[3], u[1], u[4], code = 2, xpd = TRUE)
}

plotLine(rep(1,2))
mtext(side = 2, 'Prob.',)
plotLine(c(0,1))
dev.off()