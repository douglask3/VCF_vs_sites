source('libs/plotStandardMap.r')
source("../LimFIRE/src/libs/biomeInfo.r")
library(rasterExtras)
graphics.off()

teow = '../LimFIRE/data/official_teow/official/'
cols = c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
limits = c(0, 1, 2, 5, 10, 20, 50)

dcols = c('#40004b','#762a83','#9970ab','#c2a5cf','#f7f7f7','#a6dba0','#5aae61','#1b7837','#00441b')
dlimits = c(-35, -30, -25, -20, -15, -10, -5, -2, -1, 1, 2, 5, 10, 15, 20, 25, 30, 35)

load("outputs/gridded_VCF_correction.Rd")

output.csv <- function(r, name) {
    ri = r[[1]]
    x = cbind(xyFromCell(ri, 1:length(ri)), r[])
    x = x[!mask[],]
    colnames = c("lon", "lat")
    if (ncol(x) == 3) colnames = c(colnames, "value")
        else colnames = c(colnames, "5%", "50%", "95%")
    colnames(x) = colnames
    fname = paste0("outputs/", name, ".csv")
    write.csv(x, file = fname, row.names = FALSE)
    
}
titles = c('clumping0_overlap0', 'clumpingMax_overlap0',
           'clumping0_overlapMax', 'clumpingMax_overlapMax')

mask = is.na(sum(do.call(addLayer, exps)) + control)
output.csv(control, 'control')
mapply(output.csv, exps, titles)
browser()

if (T) {
teow = readOGR(dsn = teow, layer = "wwf_terr_ecos")
ext  = extent (-180, 180, -90, 90)
xy   = abs(apply(as.matrix(bbox(ext)), 1, diff))
n    = 2
r    = raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)

## Rasterize the shapefile
biome = rasterize(teow, r, 'BIOME')
biome = raster::resample(biome, control)
biomeAssigned = biome
biomeAssigned[] = NaN
for (i in 2:length(biomes))
    biomeAssigned[any(layer.apply(biomes[[i]], function(j) biome == j))] = i

savanna = (biomeAssigned > 2) & (biomeAssigned < 6)
}

files = list.files("outputs/")
titles = c('no clumping/overlap', 'Max. clumping', 'Max. overlap', 'Max. clumping/overlap')
txt.cols = sapply(files[grepl('Savanna', files)],
                      function(i) strsplit(strsplit(i, 'Savanna-')[[1]][2], '.csv')[[1]])
control[is.na(biomeAssigned)] = NaN
plotMap <- function(r, ...) {
    r0 = r
    r[!savanna] = NaN
    plotStandardMap(r, ...)
    contour(savanna, levels = 0.5, drawlabels = FALSE, add = TRUE)
    r[!savanna] = 0
    return(r)
}
png("figs/VCF_maps.png", height = 3.7, width = 7.2, units = 'in', res = 300)
    layout(rbind(c(1, 0), c(2, 0), c(3, 4), c(5, 6), 7), heights = c(1, 0.45, 1, 1, 0.45))
    par(mar = c(0, 0, 1.1, 0))
    plotMap(100*control, cols = cols, limits = limits, title3 = 'VCF')
    par(mar = c(0, 0, 0.33, 0))
        StandardLegend(cols, limits, 100*control, extend_max = FALSE, maxLab = 100, units = '%')
    par(mar = c(0, 0, 1.1, 0)) 
    exps = lapply(exps, function(i) (i - control) *100)
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

