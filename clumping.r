
x = 1:100
y = x# rep(x, 100)
#x = rep(x, each = 100)
z = sample(c(0, 1), replace = TRUE, size = length(x)^2)

z = matrix(z, ncol = length(x))

col = c("#018571", "#a6611a")

plotFun <- function(name) {
    image(x = x-0.5, y = y-0.5, z = z, col = col, xlab = '', ylab = '', axes = FALSE)

    addBox <- function(...) {
        x = sample(0:60, size = 1)
        y = sample(0:60, size = 1)
        lines(c(x, x, x+40, x+40, x), c(y, y+40, y+40, y, y), lwd = 5,...)
    }

    addBox()
    addBox(lty = 2)
    addBox(lty = 3)
    addBox(lty = 4)
    mtext(side = 3, line = -2, name, lwd = 2, cex = 3, adj = 0.1)
}
png("figs/clumping.png", width = 7.2, height = 14.2, res = 300, units = 'in')
par(mar = rep(0.5, 4), mfcol = c(2, 1))
plotFun('Unenforced clumping')

z[] = 0
z[50:100,] = 1

plotFun('Enforced clumping')
dev.off()

