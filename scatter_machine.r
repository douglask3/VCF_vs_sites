############
## setup  ##
############
library(reldist)
library(rstanarm)
library(fields)
source("libs/is_p_star.r")
source("libs/make.transparent.r")
source("libs/make_col_vector.r")
source("libs/logit_logistic.r")
source("libs/text.units.r")


col_choices  = c("savanna" = "#d95f02", "forest" = '#1b9e77')
detail = 0.2 # 0.2
VCF_grid_size = round(c(250,250) * detail) ## size of a vcf pixel
TRB_grid_size = round(c(100,100) * detail) ## size of a trobit pixel
n_bootstraps_grid_size  = 50  # 5 number of times we'll test the uncertanty
pc_test_width = 0.01 # 0.01
vcf_clumpings = c(0, 1000  ) # How "clumped" are the trees. 0 = no clumping
CAI_shade_ps = c(0,  1)
pd_sample = seq(0.1, 0.9, 0.1)

var = "trobit_pct"
var = "CAI"

grabe_cache = TRUE

#############
## process ##
#############
## open data
dat = read.csv( 'data/trobit_vcf_comparison.csv')

#if (var == "trobit_pct")
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

plot.window <- function() {
    plot(c(0, 100), c(0, 100), xlab = "", ylab = "", type = 'n', xaxt = 'n', yaxt = 'n')
    grid()
    lines(c(0, 100), c(0, 100), lty = 2, lwd = 2)
#    lines(c(0, 100), c(0, 80), lty = 2, lwd = 2)
}

bestFit <- function(x, y, col, alpha = 0.95, lty = 1, summary = FALSE, ...) {
    if (lty == 1) lines(x, y[,2], col = col, lwd = 2)
    lines(x, y[,1], col = make.transparent(col, alpha/2), lty = lty, lwd = 1)
    lines(x, y[,3], col = make.transparent(col, alpha/2), lty = lty, lwd = 1)
    polygon(c(x, rev(x)), c(y[,1],rev(y[,3])), col = make.transparent(col, alpha), border = NA)
    if (summary) {
        x = logit(x/100)
        y = logit(y/100)
    }
}

test_clumping <- function(vcf_clumping, CAI_shade_p = 0, cont = NULL) {
    if (var == "CAI") cover_eq = find_TRB_area(CAI_shade_p)
    # prob = NULL#round((1:VCF_no_cells)/VCF_no_cells)
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
            if (vcf_clumping > 100) index = 1:nc
            else index = sample(1:VCF_no_cells, nc, FALSE, prob = prob)
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
    if (length(CAI_shade_ps) <= 1) covert_from_VCF_grid(0.5, TRB_grid_size, TRUE)
    TRB_equivilent = t(sapply(dat[, 'mvcf_pct']/100, covert_from_VCF_grid, TRB_grid_size)*100)
    x = logit(dat[, 'mvcf_pct']/100)
    VCFfit <- function(...) {
        
        y = logit(apply(TRB_equivilent,1 , function(i) sample(pcs[i>0], 1, prob = i[i>0])))
        fit = lm(y~x)
        
        # points(1/(1+exp(-x)), 1/(1+exp(-predict(fit))), pch =19, col = make.transparent('red', 0.95), cex = 1)
        #fit = glm(y ~ x, data = data.frame(x = dat[, 'mvcf_pct']/100, y = random_equiv), family = quasibinomial())
        return(c(coef(fit),  sd(predict(fit) - y)))
    }
    #dev.new()
    #plot(c(0, 1), c(0, 1), type = 'n')
    VCFcs = sapply(1:100, VCFfit)
    
    ##########
    ## plot ##
    ##########
    cols = col_choices[(dat[, "forest_type"] == "forest")+1]

    ## set up 
    plot.window()
    ## add points
    find_VCF_pd_point <- function(x, pd_sample) reldist::wtd.quantile(pcs, q = pd_sample, weight = x) * 100
    VCF_med  = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.5 )
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
         
    #dev.off()
    #browser()
    ## add trend lines
    add_trend_line <- function(type = NULL, l = 0) {
        sampleDist <- function(vs) {
            x = sample(1:ncol(vs), nrow(vs), replace = TRUE)
            x = sapply(1:length(x), function(i) vs[i, x[i]])
            return(x)
        }
        x = sampleDist(cover_eq)
        y = sampleDist(TRB_equivilent)
#TRB_mid / 100 # x is our trobit medians
        #y = pcs # y is our vcf value
        
        # if we're looking at just savanna or forest, we'll just grab those bits of data
        if (!is.null(type)) { 
            test = dat[,"forest_type"] == type
            x = x[test]; y = y[test];
            TRB_equivilent = TRB_equivilent[test,]
            col = col_choices[type]
        } else col = 'black'
        
        # cos our domain is bounded [0,1], we'll transform our data (and hence our domain) first, and we'll transform it back when we're done
        xf = logit(x/100); yf = logit(y/100);
        
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
        temp_file = paste("temp/Stan-skew22-tau12_full", vcf_clumping, CAI_shade_p, cont, type, l, '.Rd', sep = '-')
        
        #fit = nls(yf ~ a * xf + b, lower = list(a = 0, b = -9E9))
        #fit = lm(yf ~ xf, weights = w)
        
        
        if (file.exists(temp_file)) { load(temp_file) } else {            
            #fit0 = stan_glm(yf ~ xf, data = data.frame(yf, xf), prior_intercept = normal(), weights = w,
            #               family = gaussian(), prior = normal(location = 1), chain = 10, iter = 10000)
            
            fit <- stan(file = "transform.stan",
                 data = list(n = length(x), y = logit(y), x = x/100),#data.frame(y, x/100),
                 chains = 10,
                 warmup = 1000,
                 iter = 10000,
                 cores = 1,
                 refresh = 250,
                 init = rep(list(list(tau1 = 1, tau2 = 1, alpha = 0, beta = 1, sigma = 1)), 10),
                 control = list(max_treedepth = 10,
                                adapt_delta = 0.95))
            save(fit, file = temp_file)
        }
        
        ## get the best fit and confidence out, and transform it back for plotting        
        params = data.frame(extract(fit))#data.frame(fit)
        params = params[round(seq(1, nrow(params), length.out = 1000)),]
        #params = params[(nrow(params)-999):nrow(params),]
        cname = c(colnames(params), c('VCF0', 'VCFD', 'VCFE'))
        params = cbind(params, t(VCFcs))
        colnames(params) = cname
        
        addPred <- function(x, ...) {
            #x = logit(seq(0.005, 1-0.005, length.out = 100))
            
            #y = predict(fit, newdata = data.frame(xf = x),
            #            interval = "confidence")           
            
            conf = apply(params, 1, function(i)
                    i['alpha'] + i['beta']*log(x^i['tau1']/((1-x)^i['tau2'])))
            conf = t(apply(conf, 1, quantile, c(0.1, 0.5, 0.9)))            
            
            pred = c()
            for (i in 1:100) 
                pred = cbind(pred, apply(params, 1, function(i) 
                        i['alpha'] + i['beta']*log(x^i['tau1']/((1-x)^i['tau2'])) + rnorm(1, 0, i['sigma'])))
            #pred = posterior_predict(fit, newdata = data.frame(xf = x))
            pred = t(apply(pred, 1, quantile, c(0.1, 0.5, 0.9)))
            
            x = x*100
            #x    = logistic(x   )*100
            conf = logistic(conf)*100
            pred = logistic(pred)*100
        
            ## add best fit               
            bestFit(x, conf, col, ...) 
            bestFit(x, pred, col, lty = 2, ...)
            return(list(x, conf, pred))
        }
        xi = seq(0.5, 99.5, length.out = 100)/100
        addPred(xi, alpha = 0.99)
        xi = seq(min(x), max(x), length.out = 100)/100
        out = addPred(xi, alpha = 0.87)
        x = out[[1]]; conf = out[[2]]; pred = out[[3]]
        #text(x = 10, y = 100 - l,  paste("R2:", round(cor(predict(fit), yf)^2, 3)), col = col)
        
        #grad = params[,2]        
        #hx = seq(min(1, min(grad)) -0.03, max(1, max(grad)) + 0.03, 0.01)        
        #p = hist(grad, breaks = hx, plot = FALSE)$count
        #p = p[which.min(abs(hx-1))-1]/sum(p)
        #p = summary(fit)[[5]][,4][2]
        #p = paste(round(p, 3), is_p_star(p))
        #text(x = 30, y = 100 - l,  p, col = col)
        
        return(list(cbind(x, conf, pred), params))#,  cbind(posterior_interval(fit)[1:2,], coef(fit))[,c(1, 3, 2)], fit))
    }

    ## run for all and "forest", "savanna"
    fitAll = add_trend_line(l = 0)
    fitFor = add_trend_line("forest", l = 10)
    fitSav = add_trend_line("savanna", l = 20)
    
    if (length(CAI_shade_ps)>1) {
        if (CAI_shade_p == tail(CAI_shade_ps,1)) axis(4)
        
        side = 2
    } else side = 3
    if (CAI_shade_p == CAI_shade_ps[1]) {
        axis(2)
        if (vcf_clumping > 100) txt = bquote(paste( .("Clumping   "), .(bquote(infinity))))
            else  txt = paste("Clumping", vcf_clumping)
        mtext(side = side, txt, line = 2.2)
    }
    if (length(vcf_clumpings)>1 && vcf_clumping == vcf_clumpings[1]) {        
        mtext(side = 3, paste("Canopy overlap", CAI_shade_p), line = 2)
        axis(3)
    }
    if (vcf_clumping == tail(vcf_clumpings,1)) axis(1)
    return(list(fitAll, fitFor, fitSav))
}

## add the legend
graphics.off()

run4Continent <- function(cont = NULL) {
    fname = paste0("figs/tribit_vs_VCF", '-', cont, ".png")
    if (is.null(CAI_shade_ps)) CAI_shade_ps = 0
    if (var == "CAI" && length(CAI_shade_ps) > 1)
        lmat = matrix(1:(length(vcf_clumpings) * length(CAI_shade_ps)), ncol = length(CAI_shade_ps))
    else
        lmat = t(matrix(1:(length(vcf_clumpings)*2), nrow = 2))[,c(2,1)]
    
    
    heights = c(rep(1, nrow(lmat)), 0.3)
    lmat = rbind(lmat, max(lmat) + 1)
    
    png(fname, height = 6*7.2/8, width = 7.2, res = 300, units = 'in')
        layout(lmat, heights = heights)
        
        par(mar = c(1, 1, 1, 0.5), oma = c(0, 4.5, 2.5, 1.5))
        fits = lapply(CAI_shade_ps, function(CAI_shade_p) lapply(vcf_clumpings, test_clumping, CAI_shade_p, cont))
        
        plot.new()
        par(mar = rep(0, 4))    
        legend('bottom', legend = names(col_choices), col = col_choices, pch = 19, horiz = TRUE, bty = 'n')
        mtext("Trobit cover (%)", side = 1, line = -3.5)
        mtext("VCF cover (%)", side = 2, line = 2.6, outer = TRUE)
    dev.off()
    return(fits)
}

#lapply( unique(dat$continent), run4Continent)
fits = run4Continent()

fits = fits[c(1, length(fits))]
for (i in 1:length(fits)) fits[[i]] = fits[[i]][c(1, length(fits[[i]]))]

cols_cai = c("black", "red")
cols_clu = c("black", "blue")

cols_cai =  make_col_vector(cols_cai, ncols = length(fits))
cols_clu =  make_col_vector(cols_clu, ncols = length(fits[[1]]))

cols = lapply(cols_cai, function(col1) lapply(cols_clu, function(col2) make_col_vector(c(col1, col2), ncol =3)[2]))
bestFit_fun <- function(fit, y0, type, col = "black", addParams = FALSE, name = '', ...) { 
    
    bestFit(fit[[type]][[1]][,1], fit[[type]][[1]][, 2:4], 
            alpha = 0.8, col = col, summary = TRUE, ...)
    bestFit(fit[[type]][[1]][,1], fit[[type]][[1]][, 5:7],
            alpha = 0.97, col = col, lty = 2, summary = TRUE, ...)
    
    if (addParams) {
        ys = 90 - y0 - c(0, 5)
        
        text(x = rep(c(3, 16.5, 30), each = 2), y = ys,
             round(fit[[type]][[2]], 4)[c(2, 1), ], adj = c(0,0), col = col, cex = 0.67)
        text(x = -2, y = ys, c('a', 'b'), col = col, adj = c(0, 0), cex = 0.67)
    }
    fname = paste0('outputs/stan-', name, '-', col, '.csv')
    write.csv(data.frame(fit[[type]][[2]]), file = fname)
}

plotType <- function(type, name, addParams = FALSE) {
    plot.window()
    mtext(name, side = 3, line = -1)
    y0s = list(list(5, 17), list(31, 44))
    if (addParams) text(x = c(3, 18, 33), y = 90, c('5%', '50%', '95%'), adj = c(0, 0), cex = 0.67)
    mapply(function(i, j, k) mapply(bestFit_fun, i, j, k, type = type,
                                    MoreArgs = list(addParams = addParams, name = name)),
                                    fits, y0s, cols)
}

png("figs/Clumping_canopy_overlap_extremes.png", height = 6*7.2/8, width =7.2, res = 300, units = 'in')
    par(mfrow = c(2,2), mar = c(1, 1, 1, 0.5), oma = c(3, 4.5, 2.5, 1.5))
    plotType(1, "All")
    axis(2)
    plotType(2, "Forest")
    axis(1)
    plotType(3, "Savanna", addParams = FALSE)
    axis(2)
    axis(1)

    mtext("Trobit cover (%)", side = 1, line = 2.5)
    mtext("VCF cover (%)", side = 2, line = 1.5, outer = TRUE)

    plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE)
    vcf_clumpings[vcf_clumpings>100] = "~infinity~"
    bestFit(c(0.1, 0.4), cbind(c(0.8, 0.8), c(0.85, 0.85), c(0.9, 0.9)), cols[[1]][[1]])
    text.units(adj = 0, x = 0.5, y = 0.85, paste("Clumping:", vcf_clumpings[1], " \nCanopy overlap:", CAI_shade_ps[1]))

    bestFit(c(0.1, 0.4), cbind(c(0.6, 0.6), c(0.65, 0.65), c(0.7, 0.7)), cols[[2]][[1]])
    text.units(adj = 0, x = 0.5, y = 0.65, paste("Clumping:", vcf_clumpings[1], " \nCanopy overlap:", tail(CAI_shade_ps,1)))


    bestFit(c(0.1, 0.4), cbind(c(0.4, 0.4), c(0.45, 0.45), c(0.5, 0.5)), cols[[1]][[2]])
    text.units(adj = 0, x = 0.5, y = 0.45, paste("Clumping:", tail(vcf_clumpings,1), " \nCanopy overlap:", CAI_shade_ps[1]))

    bestFit(c(0.1, 0.4), cbind(c(0.2, 0.2), c(0.25, 0.25), c(0.3, 0.3)), cols[[2]][[2]])
    text.units(adj = 0, x = 0.5, y = 0.25, paste("Clumping:", tail(vcf_clumpings,1), " \nCanopy overlap:", tail(CAI_shade_ps,1)))
dev.off()