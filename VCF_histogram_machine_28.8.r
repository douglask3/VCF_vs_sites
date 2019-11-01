library(reldist)
##################
## define stuff ##
##################
#setwd("C:/Users/Rahayu Adzhar/Desktop/30.7")
graphics.off()
height_cover_file = 'data/plotcover-23.08.csv'
vcf_cover_file    = 'data/VCF_data_annual_means_21.8.csv'

# cols used for colour codeing sites. Can be less than number of sites, 
# it will figure it out later. Just provide at least 2 colours.
cols = list(forest  = c("#008837", "#7b3294"),
            savanna = c("#018571", "#a6611a"))

## Information required to find site specific VCF hight threshold
bin_width = 0.1 # Incrememnt of hight thresholds plotted and tested
height_max = 30 # max height to do analysis over
p_threshold = 0.05 # p-value signifcant threshold
n_bootstraps = 100 # how many times VHT is bootstrapped


## Information for converting VCF size pixel to height probablity distribution at NORDESTE sized pizes
detail = 1
VCF_grid_size = 250*250 * detail^2 ## size of a vcf pixel
NOR_grid_size = 100*50  * detail^2 ## size of a nordeste pixel
n_bootstraps_grid_size = 1000  # number of times we'll test the uncertanty
pc_test_width = 0.01 # fractional cover increments tested
vcf_clumping = 0 # How "clumped" are the trees. 0 = no clumping

vcf2Nor_quants = seq(0, 1, length.out = 10)
vcf2Nor_quants = head(q[-1], -1) # this is the probablity quantiles we sample to vcf -> trobit pd with

##################
## open   stuff ##
##################
height_cover = read.csv(height_cover_file, stringsAsFactors = FALSE)
vcf_covers   = read.csv(vcf_cover_file)


veg_type = height_cover[1,-1]
height_cover = height_cover[-1,]

heights = seq(0, height_max, by = bin_width)

test_cover = as.numeric(height_cover[,1])
test_cover = sapply(heights, function(i) which.min(abs(test_cover - i)))
covers  = height_cover[test_cover,-1]


vcf_covers[,'value'] = round(vcf_covers[,'value'], 1)
vcf_covers[, 'veg_type'] = t(veg_type)

######################
## process NORDESDE ##
######################
## finds the probablity density of % cover of a VCF size grid based on a trobit measurement

## Just to keep track of things, we going to assign each site a colour and stick that in 
## the vcf_cover data.frame, were we're storing a lot of the met info.

make_col_vector <- function(r,g,b,ncols,
							limits = NULL, whiteAt0 = TRUE) {    
	if (!is.null(limits)) {
		ncols = length(limits+1)
		if (whiteAt0 && class(r)=="character") {
			whiteIndex = which(r == "white" | r == "#FFFFFF")
			if (length(whiteIndex) == 0)  whiteIndex = which(r == "black" | r == "#000000")
			if (length(whiteIndex) != 0) {
				zeroIndex  = which(limits[-1]>0 & head(limits,-1)<0)

				if ( length(zeroIndex)==0) {
					if (limits[1] > 0)r = r[whiteIndex:length(r)]
						else r = r[1: whiteIndex]

					return(make_col_vector_gubbins(r,ncols = length(limits)))
				}

				negCols = make_col_vector_gubbins(r[1:whiteIndex],
												  ncols = zeroIndex+1)
				posCols = make_col_vector_gubbins(r[whiteIndex:length(r)],
												  ncols = length(limits) - zeroIndex + 1)
				return(c(negCols,posCols[-1]))
				}
			}
	}
	return(make_col_vector_gubbins(r, g, b, ncols) )
}

make_col_vector_gubbins <- function(r, g, b, ncols) {
    library(colorspace)
	if (class(r)=="character") {
		col=col2rgb(r)
		r=col[1,]/255
		g=col[2,]/255
		b=col[3,]/255
	}

	col_vec <- function(a) approx(seq(1,ncols,length.out=length(a)),a,1:ncols)$y

	cols=colorspace::RGB(col_vec(r),col_vec(g),col_vec(b))

	return(hex(cols))
}

vcf_covers[, 'cols'] = 'black'
for (type in names(cols)) {
    test = veg_type == type
    ncols = sum(test)
    col = cols[[type]]
    if (ncols == 1) col = col[1]
    else if (ncols > 2) col = make_col_vector(col, ncols = ncols)
    vcf_covers[test, 'cols'] = col
}

## function for converting VCF to NORD grid
covert_from_VCF_grid <- function(vc_pc, grid_size) {
    print(vc_pc)
    Y = rep(0,VCF_grid_size)
    nc = round(vc_pc * VCF_grid_size)
    prob = (1:VCF_grid_size)^vcf_clumping
    testY <- function(pc) {        
        index = sample(1:VCF_grid_size, nc, FALSE, prob = prob)
        Y[index] = 1

        testXpos <- function(i) {
            mn = mean(Y[i:(i+ grid_size)])
            abs(mn - pc) < 0.1
        }
        mean(sapply(1:(VCF_grid_size - grid_size), testXpos))
    }

    pcs = rep(seq(0, 1, pc_test_width), n_bootstraps_grid_size)
    
    Ys = sapply(pcs, testY)
    
    wtd.quantile(pcs, q = vcf2Nor_quants, weight = Ys)
}

## This is a slightly hacky way of doing this. What we do is basically create a copied site
# for each probablity interval, and assign a probality range for each vcf -> NOR prob pair (i.e c(0.1, 0.9), (0.2, 0.8) etc)
vcf_covers0 = vcf_covers
NOR_equivilent = sapply(vcf_covers[,'value']/100, covert_from_VCF_grid, NOR_grid_size)*100

covers_new = c()
NOR_equivilent_new = c()
vcf_covers_new = c()

nsite_split = nrow(NOR_equivilent)/2
for (site in 1:ncol(covers)) {
    for (splt in 1:nsite_split) {
        splts = c(splt, (nsite_split*2)-splt+1)
        NOR_equivilent_new = rbind(NOR_equivilent_new,  NOR_equivilent[splts, site])
        covers_new = cbind(covers_new, as.numeric(covers[, site]) * 0.8)
        vcf_covers_new = rbind(vcf_covers_new, vcf_covers[site, ])
    }
}
colnames(vcf_covers_new) = colnames(vcf_covers)
vcf_covers = vcf_covers_new
NOR_equivilent = NOR_equivilent_new
covers = covers_new

#########################
## histergram machine  ##
#########################
## useful for what some next
make.transparent <- function(col, transparency) {
     ## Semitransparent colours
     tmp <- col2rgb(col)/255
     rgb(tmp[1,], tmp[2,], tmp[3,], alpha=1-transparency)
}

## Here, for a given vcf range, we find the hieght range where vcf and NORD match
addCoverLines <- function(cover, vcf, col) {
    lines(heights, cover, col = col, lwd = 2)
    diff = abs(sapply(vcf, '-', cover))
    maxHeight = heights[which(cover == 0)[1] - 1]
    test_bound <- function(vcfi, diffi) {
        if (vcfi > cover[1]) {
            vcf_height = c(0, bin_width)
        } else {
            vcf_height = which(diffi == min(diffi))
            vcf_height = heights[vcf_height]
            vcf_height = range(vcf_height)
        }
    }
    vcf_height = range(mapply(test_bound, vcf, data.frame(diff), SIMPLIFY = FALSE))
    if (vcf_height[2] > maxHeight) vcf_height[2] = maxHeight
    alpha =  1 - exp(-diff(vcf)*0.25)
    
    polygon(c(0, rev(range(vcf_height)), 0), rep(range(vcf), each = 2), 
            col = make.transparent(col, alpha), border = NA)
    
    #lines(c(0, max(vcf_height)), c(vcf, vcf), col = col, lty = 2)
    
    
    alpha =  1 - exp(-diff(vcf_height))
    
    polygon(c(vcf_height, rev(vcf_height)), c(-1, -1, range(vcf)),
            col = make.transparent(col, alpha), border = NA)
    
    return(vcf_height)
}



png('hist_machine.png', height = 7, width = 5, res = 300, units = 'in')
    par(mfrow = c(2,1), mar = c(1, 3.5, 0, 0), oma = c(2, 0, 0.5, 1))
    vcf_height = matrix(NaN, nrow = 2, ncol = nrow(vcf_covers))
    for (type in names(cols)) {
        test1 = vcf_covers [,'veg_type'] == type
        test2 = vcf_covers0[,'veg_type'] == type
        covers_type = covers[, test1]
        
        plot(c(-1, max(heights)), c(0, max(covers_type)), type = 'n',
             xlab = '', ylab = '', xaxs = 'i', yaxs = 'i', xaxt = 'n', xlim = c(0, 26))
             
        grid()
        
        legend('topright',legend = vcf_covers0[test2,'site'], col = vcf_covers0[test2, 'cols'], lwd = 2, title = type)
        vcf_height[,test1] = mapply(addCoverLines, data.frame(covers_type),
                                   data.frame(t(NOR_equivilent))[,test1], vcf_covers[test1, 'cols'])
    }
    axis(1)
    mtext(side = 2, line = 0, 'cover (%)', outer = TRUE)
    mtext(side = 1, line = 0, 'height (m)', outer = TRUE)
dev.off()


hist.spread <- function(x, ..., bin_width = 1) {    
    diff = (apply(x, 1, diff)) == 0
    x[diff,2] = x[diff,2] + bin_width
    
    diff = 1/(apply(x, 1, diff))
    diff = diff/min(diff)
    diff =  diff/bin_width
    diff = round(diff)
    
    xs = apply(x, 1, function(i) seq(i[1], i[2], by = bin_width))
    xs =  mapply(rep, xs, diff)
    breaks = seq(0, max(unlist(xs)), by = bin_width)
    hist(unlist(xs), breaks = breaks, ...)
    return(xs)
}


png('histogram.png', height = 8, width = 8, res = 300, units = 'in')
    hist.spread(t(vcf_height), bin_width = bin_width, xlab = '', main = '',ylab = '',
                yaxt = 'n', xaxs = 'i', xlim = c(-1, max(heights)))
    grid()          
    mtext(side = 1, line = 2, 'height (m)')
dev.off()


#######################################
## height threhsold distrution test  ##
#######################################
rmse <- function(actual, predicted) {
    actual = actual# * 0.8
    error = actual - predicted
	sqrt(mean(as.numeric(error^2)))
}

wiggle <- function(cover, height0, height1, height2) {
    
    #if (height0 <= height2) height0 = height2 - 0.1
    #height0 = max(height0, height1, height2)
    if (height0 < height1 || height0 < height2) browser()
    heights_t = heights/height0
    height1 = mean(height1/height0)
    height2 = height2/height0
    
    c = (height1/(1-height1)) * ((1/height2) - 1)
    heights_t = height0*heights_t/(heights_t + c - heights_t * c)
    newCover <- function(h) cover[which.min(abs(h - heights_t))]
    sapply(heights, newCover)
}


## for all heights, calculate p-value to see if that height is the mean height threshold, and what the NMSE would be if it was
test_height <- function(height, covers, ptest = FALSE) {
    if (ptest) {
        print(height)
        vcf_height = data.frame(vcf_height)
        nsample = ncol(covers)/nsite_split
        p_value_boot <- function(i) {
            smpl = apply(vcf_height, 2, function(i) runif(1, i[1], i[2]))
            
            ts = (mean(smpl) - height)/(sd(smpl)/sqrt(nsample))
            dt(ts, df = nsample-1)
        }
        prob = mean(sapply(1:33, p_value_boot))
    } else prob = NaN
    rmse = rmse(vcf_covers[,'value'], unlist(covers[heights == height,]))
    return(c(prob, rmse))
}

tests = sapply(heights, test_height, covers, TRUE)
p_values          = tests[1,]
p_values = p_values#/max(p_values)
NMSE_actualHeight = tests[2,]


## now find the NMSE values for all heights if TROBUt sites all had the same hieght distru
heights0 = heights[apply(covers, 2, function(i) min(which(i == 0)))] - 0.1
pred_height = heights[which.max(p_values)]
heights0[heights0 < pred_height] = pred_height + 0.2

runWiggle <- function(target_height)
    mapply(wiggle, data.frame(covers), heights0, data.frame(vcf_height), target_height)


covers_NoGlobalHeight  = lapply(1:n_bootstraps, function(i) runWiggle(runif(ncol(covers), 0, heights0)))
covers_YesGlobalHeight = runWiggle(pred_height)



#NMSE_YesGlobalHeight = sapply(heights - 5.8, test_height_relative2min)
NMSE_YesGlobalHeight = sapply(heights, test_height, covers_YesGlobalHeight)[2,]
#red line
#5.8 is the mean height value for the threshold of the VCF readings (??? needs to be changed as we add more data sis)

  
  #this is for the blue line, where height is random, so we're getting random height values here
#plots of noise because we're using random values across v small data set, so 'bootstrap' it to get repeats
  #in this case ten times (1:10) and that's what yes is
  
  
NMSE_NoGlobalHeight = sapply(covers_NoGlobalHeight, function(i) sapply(heights, test_height, i)[2,])
NMSE_NoGlobalHeight = apply(NMSE_NoGlobalHeight, 1, mean)

#we have 10 values for every height increment, so get the average representative value 
  
#########################
## plot hypthosis test ##
#########################
png('height_threshold_distrubtion_test.png', height = 6, width = 7, res = 300, units = 'in')
    layout(rbind(1, 2), heights = c(0.3, 1))
    par(mar = c(0, 3.1, 0, 3.1), oma = c(3, 0, 0, 0))
    
    ## test which heights are significantly different from mean of thresholds
    plot(heights, p_values, ylim = c(min(p_values), 1), type = 'l', log = 'y', 
         xaxs = 'i', axes = FALSE, ylab = '', bg = "white")
    grid()
    axis(4)
    mtext(side = 4, line = 2, 'p-value', cex = par("cex"))
    lines(c(-9E9, 9E9), c(p_threshold, p_threshold), lty = 2)
    sig_height = range(heights[p_values > p_threshold])
    sig_height = paste0("significant heights: ", paste(sig_height, collapse = '-'), 'm')
    mtext(sig_height, side = 3, line = -2.6, adj = 0.9)
    
    index = which(p_values > 0.05)
    polygon(c(heights[index], rev(range(heights[index]))), c(p_values[index], rep(p_threshold, 2)), col = "grey")

    ## plot the curves that go and test the spread of height thresholds
    
    ## setup plot
    MaxY = max(NMSE_YesGlobalHeight, NMSE_actualHeight)
    
    plot (range(heights), c(0, MaxY), type = 'n', xlab = '', ylab = '', axes = FALSE, xaxs = 'i', yaxs = 'i')
    grid()
    axis(1)
    axis(2)
    mtext(side = 1, 'height (m)', line = 2)
    mtext(side = 2, 'NMSE', line = 2)
    
    ## I'm adding the differences between curved first so they dont get in the way of the curces themselves
    pheights = c(heights, rev(heights))
    polygon(pheights, c(NMSE_YesGlobalHeight, rev(NMSE_actualHeight)), col = make.transparent('red', 0.67), border = NA)
    polygon(pheights, c(NMSE_NoGlobalHeight , rev(NMSE_actualHeight)), col = make.transparent('blue', 0.67), border = NA)

    ## add the curves
    lines(heights, NMSE_YesGlobalHeight, lwd = 2, col = 'red'  )
    lines(heights, NMSE_NoGlobalHeight , lwd = 2, col = 'blue' )
    lines(heights, NMSE_actualHeight   , lwd = 2, col = 'black')

    ## calculate the different of actual height distrbution and run metric, and add to plot.
    diff = c(sum(abs(NMSE_YesGlobalHeight - NMSE_actualHeight)),
             sum(abs(NMSE_NoGlobalHeight  - NMSE_actualHeight)) ) * bin_width/100
    diff[3] = diff[1] / diff[2]
    
    pw = par("usr")
    lx = pw[1] + diff(pw[1:2])*0.5
    ly = pw[3] + diff(pw[3:4])*c(0.925, 0.95, 0.7, 0.65, 0.6)
    
    legend(lx, ly[1], legend = rep(' ', 3), col = make.transparent(c('red', 'blue', 'white'), 0.67), lwd = 20, bty = 'n')
    legend(lx, ly[2], legend = c( 'globally consistent', 'actual thresholds','no global threshold'), col = c('red', 'black', 'blue'), lwd = 2, bty = 'n')
    text(lx, ly[3:5], paste(c('|actual - consistent|', '|actual - no threshold|', 'ratio'), ': ', round(diff, 2), c('m', 'm', '')), pos = 4)

   
dev.off()
