############
## setup  ##
############
library(reldist)
library(fields)
source("libs/is_p_star.r")
source("libs/make.transparent.r")
source("libs/logit_logistic.r")

col_choices  = c("savanna" = "#d95f02", "forest" = '#1b9e77')
detail = 0.1
VCF_grid_size = round(c(250,250) * detail) ## size of a vcf pixel
TRB_grid_size = round(c(100,100) * detail) ## size of a trobit pixel
n_bootstraps_grid_size  = 3  # number of times we'll test the uncertanty
pc_test_width = 0.01
vcf_clumpings = c(0, 2, 4, 8) # How "clumped" are the trees. 0 = no clumping
CAI_shade_ps = c(0, 0.1, 0.2, 0.5, 1)
pd_sample = seq(0.1, 0.9, 0.1)

var = "trobit_pct"
var = "CAI"

grabe_cache = TRUE

#############
## process ##
#############
## open data
dat = read.csv( 'data/trobit_vcf_comparison.csv')

if (var == "trobit_pct")
    dat[, 'mvcf_pct'] = dat[,'mvcf_pct'] / 0.8

## finds the probablity density of % cover of a VCF size grid based on a trobit measurement
VCF_grid = matrix(0,VCF_grid_size[1], VCF_grid_size[1])
VCF_no_cells = VCF_grid_size[1] * VCF_grid_size[2]
pcs = rep(seq(0, 1, pc_test_width), n_bootstraps_grid_size)

find_TRB_area <- function(CAI_shade_p) {
    
    if (var == "CAI") {
        CAI_to_cover <- function(cai) {
            #grid = matrix(0, ncol = 10, nrow = 10)
            nc = 100#(nrow(grid) * ncol(grid))
            prob = prob = (1:nc)^CAI_shade_p
            index = rep(0, nc)
            for (i in 1:round(nc * cai)) {
                i = sample(1:nc, 1,replace = TRUE, prob = prob)
                #prob[i] = prob[i]/CAI_shade_p
                index[i] = 1
            }
            cover = sum(index)/nc
        }
        
        temp_file = paste('temp/TRB_cover_eq', 'CAIpower', CAI_shade_p, '.Rd', sep = '-')
        if (file.exists(temp_file) && grabe_cache) load(temp_file)
        else {
            cover_eq = sapply(1:100, function(i) sapply(dat[,var], CAI_to_cover))
            cover_eq = cover_eq * 100
            save(cover_eq, file = temp_file)
        }
    }
    return(cover_eq)
}


test_clumping <- function(vcf_clumping, CAI_shade_p = 0, cont = NULL) {
    if (var == "CAI") cover_eq = find_TRB_area(CAI_shade_p)
    prob = (1:VCF_no_cells)^vcf_clumping
    if (!is.null(cont)) {
        test = dat$continent == cont
        dat = dat[test, ]
        if (var == "CAI") cover_eq = cover_eq[test,]
    }
    covert_from_VCF_grid <- function(vc_pc, grid_size, plot_clumping = FALSE) {
        nc = round(vc_pc * VCF_no_cells)
        print(vc_pc)
        testY <- function(pc) {
            index = sample(1:VCF_no_cells, nc, FALSE, prob = prob)
            VCF_grid[index] = 1
            if (plot_clumping) {
                grid = list(x = 1:VCF_grid_size[1], y = 1:VCF_grid_size[2], z = VCF_grid)
                int_grid = list(seq(1, VCF_grid_size[1], 0.2), seq(1, VCF_grid_size[2], 0.2))
                loc = make.surface.grid(int_grid)
                look = interp.surface(grid, loc)
                
                image(as.surface( loc, look), col = c('#a6611a','#018571'), axes = FALSE)
                
                x = sample(1:(VCF_grid_size[1] - grid_size[1] - 1), 4, replace = FALSE)
                y = sample(1:(VCF_grid_size[2] - grid_size[2] - 1), 4, replace = FALSE)
                addbox <- function(i) {                    
                    x = x[i] + c(0, 0, grid_size[1], grid_size[1], 0)
                    y = y[i] + c(0, grid_size[2], grid_size[2], 0, 0)
                    lines(x, y, col = 'black', lty = i, lwd = 1.5)
                }
                lapply(1:4, addbox)
                return()
            }
            testXpos <- function(i, j) {
                mn = mean(VCF_grid[i:(i+ grid_size[1]),j:(j+ grid_size[1])])
                abs(mn - pc) < 0.1
            }
            
            testYpos <- function(i) sapply(1:(VCF_grid_size[2] - grid_size[2]), function(j) testXpos(i,j))
            mean(sapply(1:(VCF_grid_size[1] - grid_size[1]), testYpos))
        }
        
        testY(0.5)
        if (plot_clumping) return()
        temp_file = paste('temp/', 'climping', vcf_clumping, 'vcf_pc', vc_pc, 'pc_test_width', pc_test_width, 'nboots', 
                          n_bootstraps_grid_size, 'VCF_grid', VCF_grid_size[1], VCF_grid_size[2], 
                          'TRB_grid', TRB_grid_size[1], TRB_grid_size[2], 'cont', cont, '.csv', sep = '-')
                          
        
        if(file.exists(temp_file) && grabe_cache) {
            Ys = read.csv(temp_file)[,1]
        } else {            
            Ys = sapply(pcs, testY)
            write.csv(Ys, temp_file, row.names = FALSE)
        }
        return(Ys)
    }
    if (length(CAI_shade_ps) <= 1)covert_from_VCF_grid(0.5, TRB_grid_size, TRUE)
    TRB_equivilent = t(sapply(dat[, 'mvcf_pct']/100, covert_from_VCF_grid, TRB_grid_size)*100)
    

    ##########
    ## plot ##
    ##########
    cols = col_choices[(dat[, "forest_type"] == "forest")+1]

    ## set up 
    plot(c(0, 100), c(0, 100), xlab = "", ylab = "", type = 'n', xaxt = 'n', yaxt = 'n')
    grid()
    lines(c(0, 100), c(0, 100), lty = 2, lwd = 2)
    lines(c(0, 100), c(0, 80), lty = 2, lwd = 2)
    #text(0, c(95, 90), c('k:', 'x0:'), adj = 0)

    ## add points
    find_VCF_pd_point <- function(x, pd_sample) reldist::wtd.quantile(pcs, q = pd_sample, weight = x) * 100
    VCF_med = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.5)
    VCF_low  = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.32)
    VCF_high = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.68)
    
    if (var == "CAI") {
        TRB_mid  = apply(cover_eq, 1, quantile, 0.5)
        TRB_low  = apply(cover_eq, 1, quantile, 0.32)
        TRB_high = apply(cover_eq, 1, quantile, 0.68)
        arrows(TRB_low, VCF_med, TRB_high, VCF_med,
           length=0.05, angle=90, code=3, col = cols)
    } else {
        TRB_mid = dat[,"trobit_pct"]
    }
    
    points(TRB_mid, VCF_med, pch = 19, col = cols)
    points(TRB_mid, VCF_med)
    ## add error bars    
    arrows(TRB_mid, VCF_low, TRB_mid, VCF_high,
           length=0.05, angle=90, code=3, col = cols)
           
    ## add trend lines
    add_trend_line <- function(type = NULL, l = 0) {
        
        x = TRB_mid / 100 # x is our trobit medians
        y = pcs # y is our vcf value
        
        # if we're looking at just savanna or forest, we'll just grab those bits of data
        if (!is.null(type)) { 
            test = dat[,"forest_type"] == type
            x = x[test];# y = y[test, ];
            TRB_equivilent = TRB_equivilent[test,]
            col = col_choices[type]
        } else col = 'black'
        
        # cos our domain is bounded [0,1], we'll transform our data (and hence our domain) first, and we'll transform it back when we're done
        xf = logit(x); yf = logit(y);
        
        # do the regression
        ## When we start to consider tree clumping, we may need to start weighting points
        ## by the error bars. However, for now, all VCF/TROBUT uncertanty should be
        ## equal for all sites under transformation, so does not need considering
        
        #xf = rep(xf, ncol(TRB_equivilent))
        #yf = rep(yf, each = nrow(TRB_equivilent))
        
        yf[yf < (-10)] = -10
        yf[yf >   10 ] =  10
        #w = as.vector(unlist(TRB_equivilent))
        y  = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.5)/100
        VCF_low  = logit(apply(TRB_equivilent, 1, find_VCF_pd_point, 0.32)/100)
        VCF_high = logit(apply(TRB_equivilent, 1, find_VCF_pd_point, 0.68)/100)
        
        w = 1/(VCF_high - VCF_low)
        yf = logit(y);
        fit = lm(yf ~ xf, weights = w)
        
        ## get the best fit and confidence out, and transform it back for plotting
        x = logit(seq(min(x), max(x), length.out = 100))
        y = predict(fit, newdata = data.frame(xf = x), interval = "confidence")
        
        x = logistic(x)*100
        y = logistic(y)*100
        
        ## add best fit
        lines(x, y[,1], col = col, lwd = 2)
        polygon(c(x, rev(x)), c(y[,2],rev(y[,3])), col = make.transparent(col, 0.67), border = NA)
        text(x = 10, y = 100 - l,  paste("R2:", round(summary(fit)$r.squared, 3)), col = col)
        
        p =  summary(fit)[[5]][,4][2]
        p = paste(round(p, 3), is_p_star(p))
        text(x = 25, y = 100 - l,  p, col = col)
        return(fit)
    }

    ## run for all and "forest", "savanna"
    fitAll = add_trend_line(l = 0)
    fitFor = try(add_trend_line("forest", l = 10))
    fitSav = try(add_trend_line("savanna", l = 20))
    
    
    
    if (length(CAI_shade_ps)>1) {
        if (CAI_shade_p == tail(CAI_shade_ps,1)) axis(4)
        
        side = 2
    } else side = 3
    if (CAI_shade_p == CAI_shade_ps[1]) {
        axis(2)
        mtext(side = side, paste("Clumping", vcf_clumping), line = 2.2)
    }
    if (length(vcf_clumpings)>1 && vcf_clumping == vcf_clumpings[1]) {        
        mtext(side = 3, paste("Leaf overlap", CAI_shade_p), line = 2)
        axis(3)
    }
    if (vcf_clumping == tail(vcf_clumpings,1)) axis(1)
    return(list(fitAll, fitFor, fitSav))
}

## add the legend
graphics.off()
run4Continent <- function(cont = NULL) {
    fname = paste("figs/tribit_vs_VCF", '-', cont, ".png")
    if (is.null(CAI_shade_ps)) CAI_shade_ps = 0
    if (var == "CAI" && length(CAI_shade_ps) > 1)
        lmat = matrix(1:(length(vcf_clumpings) * length(CAI_shade_ps)), ncol = length(CAI_shade_ps))
    else
        lmat = t(matrix(1:(length(vcf_clumpings)*2), nrow = 2))[,c(2,1)]
    
    
    heights = c(rep(1, nrow(lmat)), 0.3)
    lmat = rbind(lmat, max(lmat) + 1)
    
    png(fname, height = sum(heights)*2.33, width = 3.75 * ncol(lmat), res = 300, units = 'in')
        layout(lmat, heights = heights)
        
        par(mar = c(1, 1, 1, 0.5), oma = c(0, 4.5, 2.5, 1.5))
        fits = lapply(CAI_shade_ps, function(CAI_shade_p) lapply(vcf_clumpings, test_clumping, CAI_shade_p, cont))
        
        plot.new()
        par(mar = rep(0, 4))    
        legend('bottom', legend = names(col_choices), col = col_choices, pch = 19, horiz = TRUE, bty = 'n')
        mtext("Trobit cover (%)", side = 1, line = -3.5)
        mtext("VCF cover (%)", side = 2, line = 2.6, outer = TRUE)
    dev.off()
}

#lapply( unique(dat$continent), run4Continent)
fits = run4Continent()