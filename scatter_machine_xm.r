library(reldist)
#library(rstanarm)
library(fields)
library(rstan)
source("libs/is_p_star.r")
source("libs/make.transparent.r")
source("libs/make_col_vector.r")
source("libs/logit_logistic.r")
source("libs/text.units.r")
source("libs/atans.r")
graphics.off()

trbFile = "data/trobit_vcf_comparison.csv"
satFile = "/prj/CROPNET/FFG/tcc30/"
vegType = c('All' = 'black', "savanna" = "#d95f02", "forest" = '#1b9e77')

clumpSize = 30
overlSize = 100

nboots = 100
tfile0 = 'temp/bootstrappingXm-test'

satPatt = c('plot_stats', 'buff120')
satLoc = c('Af', 'Oz', 'SA')

satFiles = list.files(satFile, pattern = satPatt[1], full.names=TRUE)

trbDat = read.csv(trbFile)[,c('site', 'CAI', 'size', 'forest_type')]

openCSV <- function(file) {
    dat = read.csv(file, header = FALSE)
    dat = dat[,!apply(is.na(dat), 2, all)]
    dat = t(dat)
    hd = dat[1,]
    dat = matrix(as.numeric(dat[-1, ]), ncol = ncol(dat))
    
    colnames(dat) = hd
    return(dat)
}

openLoc <- function(loc) {
    files = satFiles[grepl(loc, satFiles) & grepl(satPatt[2], satFiles)]
    sat = lapply(files, openCSV)
    do.call('+', sat)
}

sat = lapply(satLoc, openLoc)
nr = max(sapply(sat, nrow))

sizeSat <- function(dat) 
    dat = rbind(dat, matrix(NA, ncol = ncol(dat), nrow = nr-nrow(dat)))

sat = lapply(sat, sizeSat)
sat = do.call(cbind, sat)
grid = matrix(0, ncol = clumpSize, nrow = clumpSize)
xgrid = matrix((1:clumpSize)-clumpSize/2-0.5, ncol = clumpSize, nrow = clumpSize)
ClumpArea = clumpSize*clumpSize
clumping <- function(id, p = 0) {
    tfile = paste0(tfile0, '-slumping2-', id, '-', '-', p, clumpSize, '-', nboots, '.Rd')
    
    if (file.exists(tfile)) { 
        load(tfile)
        return(outs)
    }
    sat = sat[,which(colnames(sat) == id)]/80
    sat[sat >1] = 1
    trSize = trbDat[which(trbDat[,1] == id), 3]
    x = seq(0, 2, length.out = length(sat))
    y = 0
    xgrid = xgrid /(trSize * 100); ygrid = t(xgrid)
    bootstrap <- function(...) {
        theta =  atans(x, y, units = "radians") + runif(1, 0, 2*pi) 
        d = sqrt(x^2 + y^2)
        forSat <- function(d, th, sat) { 
            if (is.na(sat)) return(c(NaN, NaN))     
            if (p == 0) prob = rep(1, ClumpArea) else prob = (1:ClumpArea)^p     
            index = sample(1:ClumpArea, round(ClumpArea*sat), FALSE, prob = prob)
            dist = grid
            grid[index] = 1
            xgrid = xgrid + d*sin(th); ygrid = ygrid + d*cos(th)
            test = (abs(xgrid) + abs(ygrid))<(trSize)/2
            
            return(c(mean(grid[test]), sum(test)))
        }
       
        out = mapply(forSat, d, theta, sat)
        sum(out[1,] * out[2,], na.rm = TRUE)/sum(out[2,], na.rm = TRUE)
    }
    outs = sapply(1:nboots, bootstrap)
    save(outs, file = tfile)
    return(outs)
}

sclumped1 = lapply(colnames(sat), clumping, p = 100)
sclumped0 = lapply(colnames(sat), clumping)

overlSize = 10
nc = overlSize^2
overlapping <- function(id, p = 0) {
    tfile = paste0(tfile0, '-overlapping2-', id, '-', '-', p, nc, '-', nboots, '.Rd')
    
    if (file.exists(tfile)) { 
        load(tfile)
        return(outs)
    }
    bootstrap <- function(...) {
        index = which(trbDat[,1] == id)
        cai = trbDat[index,2]
        prob = (1:nc)^p
        index = rep(0, nc)
    
        for (i in 1:round(nc * cai)) {
            i = sample(1:nc, 1,replace = TRUE, prob = prob)
            index[i] = 1
        }
        mean(index)
    }
    outs = sapply(1:nboots, bootstrap)
    return(outs)
}

soverlap0 = lapply(colnames(sat), overlapping)
soverlap1 = lapply(colnames(sat), overlapping, 1)

addBars <- function(x, y, ...) {
    addBar <- function(i) {
        x = x[,i]; y = y[,i]
        lines(x, y, ...)
        addHand <- function(FUN = min) {
            if (x[1] == x[2]) lines(x[1] + c(-1, 1), rep(FUN(y), 2))
            else              lines(rep(FUN(x), 2), y[1] + c(-1, 1))
        }
        lapply(c(min, max), addHand)
    }
    lapply(1:ncol(x), addBar)
}

plotAndFit <- function(trb, sat, id) {
    mnSD <- function(ls) 
        sapply(ls, quantile, c(0.159, 0.5, 0.841))
    
    trb = mnSD(trb)*100; sat = mnSD(sat)*100
    plot(c(0, 100), c(0, 100), type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
    
    forVeg <- function(veg, col) {
        if (veg != 'All') browser()
        xf = logit(trb/100)
        yf = logit(sat/100)

        w =  sqrt((xf[3,]-xf[1,])^2 + (yf[3,]-yf[1,])^2)

    #browser()
    
        temp_file = paste0(tfile0, veg, '-rstan2-', id, '.Rd')
        if (file.exists(temp_file)) { load(temp_file) } else {             
            fit <- stan(file = "transform.stan",
                        data = list(n = ncol(sat), y = logit(sat[2,]/100),
                                    x = trb[2,]/100),#data.frame(y, x/100),
                        chains = 10,
                        warmup = 1000,
                        iter = 10000,
                        cores = 1,
                        refresh = 250,
                        init = rep(list(list(tau1 = 1, tau2 = 1, alpha = 0, 
                                             beta = 1, sigma = 1)), 10),
                        control = list(max_treedepth = 10,
                                       adapt_delta = 0.95))
            save(fit, file = temp_file)
        }
        params = data.frame( rstan::extract(fit))
        params = params[sample(1:nrow(params),1000, replace = FALSE),]
    
        x = seq(0, 1, 0.01)
        xf = logit(x)
        reconFit <- function(ps)
            logistic(ps['alpha'] +log((x^ps['tau1'])/(((1-x)^ps['tau2']))) * ps['beta'])
        
    
        y = apply(params, 1 , reconFit)
        y = apply(y, 1, quantile, c(0.05, 0.5, 0.95))*100
        x = x*100
        polygon(c(x, rev(x)), c(y[1,], rev(y[3,])), 
                col = make.transparent(col, 0.8), border = col)
        lines(x, y[2,], lwd = 2)
    }
    mapply(forVeg, names(vegType ), vegType)
    grid()
    points(trb[2,], sat[2,], pch =19, cex = 0.1)
    addBars(trb[c(2,2), ], sat[c(1, 3), ])
    addBars(trb[c(1,3), ], sat[c(2, 2), ])
    lines(c(0, 100), c(0,100), lwd = 2, lty = 2)
}

png("figs/scatter_machine.png", height = 7, width = 7.2, units = 'in', res = 300)
#par(mar = rep(0.25, 4), oma  = c(2.5, 3.5, 3.5, 2.5), mfrow = c(2,2))
layout(rbind(1:2, 3:4, 5), heights = c(1, 1, 0.3))
par(mar = c(1, 1, 1, 0.5), oma = c(0, 4.5, 2.5, 1.5))

plotAndFit(soverlap0, sclumped0, 'overlap0-sclumped0')
axis(2)
axis(3)
mtext(side = 3, line = 2, 'Unenforced clumping')
mtext(side = 2, line = 2, 'Unenforced overlap')

plotAndFit(soverlap0, sclumped1, 'overlap0-sclumped1')
axis(3)
axis(4)
mtext(side = 3, line = 2, 'Enforced clumping')

plotAndFit(soverlap1, sclumped0, 'overlap1-sclumped0')
axis(1)
axis(2)

mtext(side = 2, line = 2, 'Enforced overlap')
plotAndFit(soverlap1, sclumped1, 'overlap1-sclumped1')
axis(1)
axis(4)
#sat = do.call(rbind, lapply(satFile, read.csv))

mtext("Sexton cover (%)", side = 2, line = 2.6, outer = TRUE)
plot.new()
mtext("Trobit cover (%)", side = 1, line = -3.5)

dev.off()



