source('libs/plotStandardMap.r')
source("../LimFIRE/src/libs/biomeInfo.r")
library(rasterExtras)
library(gitBasedProjects)
graphics.off()

#teow/ = '../LimFIRE/data/official_teow/official/'
cols = c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
limits = c(0, 1, 2, 5, 10, 20, 50)

dcols = c('#40004b','#762a83','#9970ab','#c2a5cf','#f7f7f7','#a6dba0','#5aae61','#1b7837','#00441b')
dlimits = c(-20, -16, -12, -8, -4, -2, -1, 1, 2, 4, 8, 12, 16, 20)
vegTypeNames = c("EG NL Forest", "EG BL Forest", "Dec NL Forest",
                    "Dec BL Forest", "Mixed Forest",
                    "Close Shrub", "Open Shrub", "Woody Savanna", "Savannas", "Grassland")
forMask <- function(id) {
    load(paste0("outputs/gridded_VCF_correction-stan-test2-mask", id, ".Rd"))
    ### histogram
    
    conHist = exps[[1]][-(1:13)][1:10]
    expHist = lapply(exps, function(i) i[-(1:13)][11:20])
    expHist = lapply(1:length(expHist[[1]]), function(i) lapply(expHist, function(j) j[[i]]))
     
    lty_area = unlist(exps[[1]][4:13])
    plotLUtype <- function(con, exp, name) {
        plot(c(0, 100), c(0, 1), yaxt = 'n', xlab = '', cex = 1000,
             pch = 19, col = "white", xaxs = 'i')
        mtext(side = 3, line = -1.5, adj = 0.9, name)
        
        x = seq(0.5, 99.5, 1)
        addPoly <- function(y, col = "black") {
            test = y > 0
            y = y[test]/max(y)
            x = x[test]
            polygon(c(x, rev(x)), c(y, rep(0, length(y))), border = col,
                   col = make.transparent(col, 0.67))
        }
        addPoly(con)
        #mapply(addPoly, exp, c("red", "green", "blue", "pink"))
        #browser()
    }
    if (id == 1 && F) {
        png("figs/vegTypes-comb.png", width = 7.2, height = 7.2*5/3, res = 300, units = 'in')
            par(mfrow = c(5, 2), mar = rep(0.5, 4), oma = c(3,3,0, 0))
            browser()
            vegCom = list(1:5, 6, 7, 8, 9, 10)
            vegTypeNames = c('Forest', vegTypeNames[6:10])
            combine <- function(vt, Hist) Reduce('+', Hist[vt])
            conHisti = lapply(vegCom, combine, conHist)
            #expHisti = lapply(vegCom, combine, expHist)
            
            mapply(plotLUtype, conHisti, expHist[1:5],
                  vegTypeNames)
        
        dev.off()
    }
    out = pout = nout = c()
    mask = !is.na(do.call('sum', c(control, lapply(exps, function(i) i[[1]]))))
    output.csv <- function(r, name) {
        ri = r[[1]]       
        
        if (is.null(out)) {
            out = cbind(xyFromCell(ri, 1:length(ri)), ri[])
            out = out[mask[],]
            colnames(out) = c("lon", "lat", "control")
        } else {
            x = ri[mask] + out[,3]
            colnames(x) = paste0(name, '-',c("10%", "50%", "90%"))
            out = cbind(out, x)
        }
        out <<- out
        #if (ncol(x) == 3) colnames = c(colnames, "value")
        #    else colnames = c(colnames, "5%", "50%", "95%")
        #colnames(x) = colnames
        #fname = paste0("outputs/", name, "mask", id, ".csv")
        #

        if (is.list(r) && length(r) >= 3) {
            outNegPos <- function(i, x0) {
                #fname = paste0("outputs/", "landTypeChange-", name, "-",
                #                pname, "mask", id, ".csv")
                x = r[[i + 1]][[1]]
                
                if (id == 2) x = x*250*250/(1000*1000*1000000)
                else x = sweep(x, 2, unlist(r[4:13]), '/')
                colnames(x) = vegTypeNames
                rownames(x) = paste0(name, '-',rownames(x))
                
                return(rbind(x0, x))
                #write.csv(x, file = fname)
            } 
            nout <<- outNegPos(1, nout)
            pout <<- outNegPos(2, pout)
        }
    
    }
    titles = c('clumping0_overlap0', 'clumpingMax_overlap0',
               'clumping0_overlapMax', 'clumpingMax_overlapMax')

    output.csv(control, 'control')
    mapply(output.csv, exps, titles)
    
    write.csv(out, file = paste0('outputs/output4histogram-mask_', id, '-git_rev_', gitVersionNumber(), '.csv'), row.names = FALSE)
    
    write.csv(rbind(nout, lty_area), file = paste0('outputs/output_negetativeChange-mask_', id, '-git_rev_', gitVersionNumber(), '.csv'), row.names = TRUE)
    write.csv(rbind(pout, lty_area), file = paste0('outputs/output_positiveChange-mask_'  , id, '-git_rev_', gitVersionNumber(), '.csv'), row.names = TRUE)
    
    files = list.files("outputs/", pattern = "stan")
    titles = c('no clumping/overlap', 'Max. clumping', 'Max. overlap', 'Max. clumping/overlap')
    txt.cols = sapply(files[grepl('All', files)],
                      function(i) strsplit(strsplit(i, 'All-')[[1]][2], '.csv')[[1]])
    
    #control[is.na(biomeAssigned)] = NaN
    plotMap <- function(r, ...) {
        r0 = r
    
        plotStandardMap(r, ...)
        #contour(savanna, levels = 0.5, drawlabels = FALSE, add = TRUE)
        
        return(r)
    }

    png(paste0("figs/VCF_maps-mask",id,".png"),
        height = 3.7, width = 7.2, units = 'in', res = 300)
        layout(rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), 9), heights = c(1, 0.45, 1, 1, 0.45))
        par(mar = c(0, 0, 1.1, 0))
        
        plotMap(100*control, cols = cols, limits = limits, title3 = 'VCF')

        exps = lapply(exps, function(i) 100*i[[1]])# (i - control) *100)
        sig = layer.apply(exps, function(i)
                    (sum(i[[c(1,3)]] < 0)==2) + (sum(i[[c(1, 3)]]>0) == 2))
        
        sig_cols = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb',
                               '#41b6c4','#1d91c0','#225ea8','#253494','#081d58')
        plotMap(sum(sig), cols = sig_cols, 
                limits = 0.5:3.5, title3 = 'VCF')
        par(mar = c(0, 0, 0.33, 0))        
        StandardLegend(cols, limits, 100*control, extend_max = FALSE, maxLab = 100, units = '%')
        
        StandardLegend(sig_cols, 0:3, sig, labelss = c(0, paste(1:4, '                       ')),
                        extend_max = FALSE)
        par(mar = c(0, 0, 1.1, 0)) 
        
        ehist = mapply(plotMap, exps, title3 = titles, txt.col = txt.cols,
            MoreArgs = list(cols = dcols, limits = dlimits))
        
       
        control = control*100       
        ehist = lapply(ehist, function(i) i + control)
        par(mar = c(0, 0, 0.33, 0))
        StandardLegend(dcols, dlimits, exps[[2]], extend_min = TRUE, units = '%')
    dev.off()
    breaks = seq(0, 130, 1)

    getHist <- function(r)   { 
        if (nlayers(r) == 3) r = r[[2]]
        # r = raster::crop(r, c(-20, 50, -30, 30)) 
        out = hist(r[][r[] > 9], plot = FALSE, breaks = breaks)$counts
        out = out#/max(out)
    }

    chist = getHist(control)
    ehist = lapply(ehist, getHist)

    my = max(chist, unlist(ehist))

    png("figs/hist.png", height = 5, width = 5, res = 300, units = 'in')
        plot(c(10, 100),c(0, my), type = 'n', yaxt = 'n', xaxs = 'i', xlab = '', ylab = '')
        mtext('VCF cover (%)', side = 1, line = 2)

        x = breaks[-1] + diff(breaks)/2
        polyFun <- function(y, line = TRUE, fill = FALSE, ...) {
            y = smooth.spline(x, y, spar=0.35)[[2]]#predict(loess(y~x))
            if (line) lines(x, y, ...)
            if (fill) polygon(c(x, rev(x)), c(y, rep(0, length(x))), border = NA, ...)
        }
        polyFun(chist, col = "grey", fill = TRUE, line = FALSE)
        mapply(polyFun, ehist, col = make.transparent(txt.cols, 0.0), fill = FALSE, lwd = 2)
    dev.off()

    png("figs/hist0.png", height = 5, width = 5, res = 300, units = 'in')
        plot(c(10, 100),c(0, my), type = 'n', yaxt = 'n', xaxs = 'i', xlab = '', ylab = '')
        mtext('VCF cover (%)', side = 1, line = 2)
    polyFun(chist, col = "black", fill = TRUE, line = FALSE)
    dev.off()
}

forMask(1)
forMask(2)
