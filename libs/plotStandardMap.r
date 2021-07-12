source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs('../rasterextrafuns/rasterPlotFunctions/R/')
source("libs/return_multiple_from_functions.r")
source("libs/sd.raster.r")

library(raster)

library(plotrix)
library(mapdata)
library(mapplots)

library(rgdal)

#SA_ste <- readOGR(dsn = "data/South_America", layer = "South_America")
rivers <- readOGR(dsn = "data/majorrivers_0_0", layer = "MajorRivers")

StandardLegend <- function(cols, limits, dat, rightx = 0.95, extend_max = TRUE, oneSideLabels = NA, add = FALSE, ylabposScling = 1,...) {
    if (add)        
        plot_loc = c(0.41, rightx, 0.1, 0.13)
    else 
        plot_loc = c(0.05, rightx, 0.3, 0.56)
    add_raster_legend2(cols, limits, dat = dat, add = add,
                       transpose = FALSE, srt = 0, oneSideLabels= oneSideLabels,
                       plot_loc = plot_loc,
                       ylabposScling = ylabposScling, extend_max = extend_max, ...)
}

lineBox <- function(x, y, ...) 
    lines(c(x[1], x[2], x[2], x[1], x[1]), c(y[1], y[1], y[2], y[2], y[1]),...)

plotStandardMap <- function(r, cols, limits, e = NULL, add_legend = FALSE,
                            limits_error = c(0.5, 0.500000001),
                            title2 = '', title3 = '', txt.col = "black",
                            xlim = c(-120, 160), ylim = c(-30, 30),quick = TRUE, ...) {
     
    if (nlayers(r) > 1 && is.null(e)) {
        if (nlayers(r) == 3) {
            if (any(r[] <0, na.rm = TRUE) && any(r[] > 0, na.rm = TRUE)) {
                
                e = 1-((sum(r<0) == 3) + (sum(r>0) == 3))
                
                #e = r[[2]]
                #e[!is.na(r[[2]])] = 0
                #e = abs(r[[3]] - r[[1]])/max(abs(r[[c(1,3)]]))/2
                #e[r[[3]]>0 & r[[1]] <0] = 1                
            } else  e = 1-r[[1]]/r[[3]]
            
            r = r[[2]]
        } else {
            e = sd.raster(r)
            r = mean(r)
        }
    } 
    
    r[r>9E9] = NaN
    if (!is.null(e)) e[is.na(r)] = NaN
    
    plot(xlim, ylim, xlab = '', ylab = '', axes = FALSE, type ='n', xaxs = 'i', yaxs = 'i')
    grid()
    plot_raster_from_raster(r, e = e,interior = FALSE,#coast.lwd = NULL,
                            cols = cols, limits = limits, add_legend = FALSE,
                            quick = quick, ePatternRes = 40, ePatternThick = 0.5,
                            limits_error = limits_error, add = TRUE, ...)
    
    #plot(rivers, col = c(rep("#FFFFFF00", 9), "black", rep("#FFFFFF00", 88)), add = TRUE, lwd = 2.5)
    #plot(rivers, col = c(rep("#FFFFFF00", 9), "white", rep("#FFFFFF00", 88)), add = TRUE, lwd = 0.5)
    #lines(SA_ste)
    #lines(rivers, col = "white", lwd = 0.67)
    #lines(SA_ste, lty = 2)
    #lineBox(-c(71.25, 63.75), -c(11.25,  6.25))     
    #lineBox(-c(61.25, 53.75), -c(11.25,  6.25))   
    #lineBox(-c(48.25, 43.25), -c( 8.75,  1.25))
    #lineBox(-c(66.25, 58.75), -c(18.75, 13.75))
    addCoastline(quick = quick, r, ...)
    polygon(c(-62.5, -35, -35, -62.5), c(-56, -56, -50, -50), border = NA, col = "white")
    mtext(title3, adj = 0.01, line = 0, col = txt.col)
    mtext(title2, side = 2, line = -0.75, col = txt.col)
    if (add_legend) {
        add_raster_legend2(cols, limits, dat = r,
                           transpose = FALSE, srt = 0, oneSideLabels= FALSE,
                           plot_loc = c(0.35, 0.99, 0.09, 0.12),  ylabposScling=0.8, ...)
    }
}

addCoastline <- function(quick = TRUE, mask = NULL, ...) {
    
    if (is.null(mask)) mask = raster('data/seamask.nc') else mask = (is.na(mask) + raster('data/seamask.nc'))>1
   
    if (!quick) mask = raster::disaggregate(mask, fact = 5, method = "bilinear")
    mask = mask>0.5
    
    plot_raster_from_raster(mask+1, add = TRUE, 
                             cols = c("white", "transparent"),readyCut = TRUE,
                             limits =  NULL, quick = TRUE, interior = FALSE, 
                             coast.lwd = NULL, add_legend = FALSE, ...)
    #
    #contour(mask, add = TRUE, drawlabels = FALSE, lwd = 0.5)  
    
    ployBox <- function(x, y)
        polygon(c(x[1], x[2], x[2], x[1]), c(y[1], y[1], y[2], y[2]), col = "white", border = "white")
        
    ployBox(c(-180, -90), c(-60, 0))
    ployBox(c(-180, -120), c(-60, 25))
    ployBox(c(-50, -19), c(10, 25))
    ployBox(c(-50, -13.5), c(27.5, 34))
    ployBox(c(115, 125), c(-8, -7))
    ployBox(c(104, 111), c(2.5, 8))
    ployBox(c(122, 128), c(2.5, 5)) 
}
