source("libs/make_col_vector.r")

##################
## define stuff ##
##################
height_cover_file = 'data/histogram/plotcover.csv'
vcf_cover_file    = 'data/histogram/VCF_data_annual_means.csv'

# cols used for colour codeing sites. Can be less than number of sites, 
#   it will figure it out later. Just provide at least 2 colours.
cols = c('#1b9e77','#d95f02','#7570b3')

## number of bins (how detailed) in your histergram
nbins = 20

##################
## open   stuff ##
##################
height_cover = read.csv(height_cover_file)
vcf_covers   = read.csv(vcf_cover_file)

heights = height_cover[,1]
covers  = height_cover[,-1]

###########################################################
## functions required for plotting. Dont look to closley ##
###########################################################

cols = make_col_vector(cols, ncols = ncol(covers))

addCoverLines <- function(cover, vcf, col) {
    
    lines(heights, cover, col = col, lwd = 2)
    
    diff = abs(cover - vcf)
    vcf_height = which(diff == min(diff))
    vcf_height = heights[vcf_height]
    vcf_height = range(vcf_height)
    
    lines(c(0, max(vcf_height)), c(vcf, vcf), col = col, lty = 2)
    
    alpha =  1 - exp(-diff(vcf_height))
    
    polygon(c(vcf_height, rev(vcf_height)), c(-1, -1, vcf, vcf),
            col = make.transparent(col, alpha), border = NA)
            
    return(vcf_height)
}

hist.spread <- function(x, ..., bin_width = 1) {
        
    diff = 1/(apply(x, 1, diff))
    diff = diff/min(diff)
    diff =  diff/bin_width
    diff = round(diff)
    
    xs = apply(x, 1, function(i) seq(i[1], i[2], by = bin_width))
    xs =  mapply(rep, xs, diff)
    
    hist(unlist(xs), ...)
    return(xs)
}

##################
## set up plot  ##
##################
png('hist_machine.png', height = 8, width = 6, res = 300, units = 'in')
    par(mfrow = c(2, 1), mar = c(1, 3.5, 0, 0), oma = c(2, 0, 0, 1))
    
    plot(c(0, max(heights)), c(0, max(covers)), type = 'n',
         xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', xaxt = 'n')
    grid()
    mtext(side = 2, line = 2, 'cover (%)')
    
    legend('topright',legend = vcf_covers[,'site'], col = cols, lwd = 2)

    vcf_height = mapply(addCoverLines, covers, vcf_covers[,'value'], cols)
    
    xs = hist.spread(t(vcf_height), nbins, bin_width = 0.1, xlab = '', main = '',ylab = '',
                yaxt = 'n', , xaxs = 'i', xlim = c(0, max(heights)))
    grid()          
    mtext(side = 1, line = 2, 'height (m)')
dev.off()