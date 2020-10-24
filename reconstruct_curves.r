
titles = c('Unenforced clumping/overlap', 'Enforced clumping', 'Enforced overlap', 'Enforced clumping/overlap')
paramCols = c('#000000', '#0000BC', '#BC0000', '#BC00BC')
param_files = paste0('outputs/stan-All-', paramCols, '.csv')
paramCols[1] = '#FFFFFF'
#paramCols[2] = '#0000FF'
LCNames = c('Forest', 'Close Shrub', 'Open Shrub', 'Woody Savanna', 'Savannas', 'Grassland')
hist_files = paste0('outputs/histLC-', LCNames, '.csv')



paramMe <- function(param_file, col, addX, yside = 2, title = '') {
    params = read.csv(param_file)

    hist = read.csv(hist_files[1])
    x = hist[,2]/100

    corFull <- function(i) {
        out = i['VCF0'] + logit(x) * i['VCFD']
        out = (out - i['alpha'])/i['beta']
        out = out1 = exp(out)
        test = is.infinite(out)
        xt = seq(0.00, 1, 0.005)
        x1 = xt^i['tau1']
        x2 = (1-xt)^i['tau2']
        out = sapply(out, function(i) xt[which.min(abs(i * x2 - x1))])
     
        
        out[test] = 1 
        out = logistic((logit(out) - i['VCF0'])/i['VCFD'])
        out = out - x
        return(out)
    }

    outs = apply(params, 1, corFull)
     
    VegTypeMe <- function(hist_file) {
        hist = read.csv(hist_file)
        LUouts = apply(outs, 2, '*', hist[,3])/sum(hist[,3])


        sumQuat <- function(test = 'nan') {
            sumF <- function(i) {
                if (test == 'pos') i = i[i>0]
                if (test == 'neg') i = i[i<0]
                
                sum(i)
            }
            out = apply(LUouts, 2, sumF)
            out = quantile(out, c(0.1, 0.5, 0.9))
            return(out)
        }

        LUouts_neg = sumQuat('neg')
        LUouts_pos = sumQuat('pos')
        LUouts_net = sumQuat('')
        return(cbind(LUouts_pos, LUouts_neg, LUouts_net))
    }

    pnn = lapply(hist_files, VegTypeMe)

    median = sapply(pnn, function(i) i[2,]*100)
    lower = sapply(pnn, function(i) i[1,]*100)
    upper = sapply(pnn, function(i) i[3,]*100)

    barCenters  = barplot(median, beside = TRUE, ylim = c(-10, 30),
                          col = c('#FFFFFF', '#BBBBBB', '#777777'), yaxt = 'n')
    if (col == '#FFFFFF') col = "#000000"                      
    barplot(median, beside = TRUE, ylim = c(-10, 30), add = T,
                          col = paste0( col, '44'), yaxt = 'n')                      
    axis(yside)
    mtext(title,  side = 3, line = -1, adj = 0.1, col = col)
    
    
    arrows(barCenters, lower, barCenters,
           upper, lwd = 1.5, angle = 90,
           code = 3, length = 0.05)
    
    lapply(seq(-10, 30, 5), function(i) lines(c(-9E9, 9E9), c(i, i), lty = 2, col = 'grey'))

    if (addX) {
        text(x = barCenters, y =  -11, adj = 1, srt = 45, xpd = NA,
             unlist(lapply(LCNames, c, c('', ''))))
        col = "black"
    } else {
        text(x = barCenters, y =  31, adj = 0, srt = 45, xpd = NA,
             unlist(lapply(LCNames, c, c('', ''))))
        col = c("black", "black", "white")
    }
    
    
}

png('figs/barPlotCheck.png', height = 9, width = 7.2, res = 300, units = 'in')
    par(mfrow = c(2, 2), mar = c(1, 0.5, 1, 0.5), oma = c(4, 2.7, 4, 1.5))
    mapply(paramMe, param_files, paramCols, c(F, F, T, T), c(2, 4, 2, 4), titles)
    mtext(side = 2, outer = TRUE, 'Change in cover (%)', line =1.5)
dev.off()
