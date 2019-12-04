############
## setup  ##
############
library(reldist)
library(fields)
source("libs/is_p_star.r")
source("libs/make.transparent.r")

col_choices  = c("savanna" = "#d95f02", "forest" = '#1b9e77')
detail = 0.25

VCF_grid_size = round(c(250,250) * detail) ## size of a vcf pixel
TRB_grid_size = round(c(100,100) * detail) ## size of a trobit pixel
n_bootstraps_grid_size  = 10  # number of times we'll test the uncertanty
pc_test_width = 0.01
vcf_clumpings = c(0, 2, 4, 8) # How "clumped" are the trees. 0 = no clumping
pd_sample = seq(0.1, 0.9, 0.1)

grabe_cache = TRUE

## some required functions
logit <- function(x) log(x/(1-x))
logistic <- function(x) 1/(1+exp(-x))

#############
## process ##
#############
## open data
dat = read.csv( 'data/trobit_vcf_comparison.csv')
dat[, 'mvcf_pct'] = dat[,'mvcf_pct'] / 0.8


## finds the probablity density of % cover of a VCF size grid based on a trobit measurement
VCF_grid = matrix(0,VCF_grid_size[1], VCF_grid_size[1])
VCF_no_cells = VCF_grid_size[1] * VCF_grid_size[2]
pcs = rep(seq(0, 1, pc_test_width), n_bootstraps_grid_size)

test_clumping <- function(vcf_clumping) {
    prob = (1:VCF_no_cells)^vcf_clumping
    
    covert_from_VCF_grid <- function(vc_pc, grid_size, plot_clumping = FALSE) {
        print(vc_pc)
        nc = round(vc_pc * VCF_no_cells)
        
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
                          'TRB_grid', TRB_grid_size[1], TRB_grid_size[2], '.csv', sep = '-')
                          
        
        if(file.exists(temp_file) && grabe_cache) {
            Ys = read.csv(temp_file)[,1]
        } else {
            
            Ys = sapply(pcs, testY)
            write.csv(Ys, temp_file, row.names = FALSE)
        }
        return(Ys)
        #wtd.quantile(pcs, q = pd_sample, weight = Ys)
    }
    covert_from_VCF_grid(0.5, TRB_grid_size, TRUE)
    TRB_equivilent = t(sapply(dat[, 'mvcf_pct']/100, covert_from_VCF_grid, TRB_grid_size)*100)
    
    
    ## for each vcf %, we get back the 0.1, 0.5, and 0.9 of the prob density function.
    ## 0.5 will be the point, while 0.1 to 0.9 will be our error bars
    
    #colnames(TRB_equivilent) = paste0('trob_vcf', pcs)

    ## tag on vcf % and forest type. We then have all the info we need for plotting
    #TRB_equivilent = cbind(dat[, c("trobit_pct", 'forest_type')], TRB_equivilent)


    ##########
    ## plot ##
    ##########
    cols = col_choices[(dat[, "forest_type"] == "forest")+1]

    ## set up 
    plot(c(0, 100), c(0, 100), xlab = "", ylab = "", type = 'n', xaxt = 'n')
    grid()
    lines(c(0, 100), c(0, 100), lty = 2, lwd = 2)
    #text(0, c(95, 90), c('k:', 'x0:'), adj = 0)

    ## add points
    find_VCF_pd_point <- function(x, pd_sample) reldist::wtd.quantile(pcs, q = pd_sample, weight = x) * 100
    VCF_med = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.5)
    
    points(dat[,"trobit_pct"], VCF_med, pch = 19, col = cols)
    points(dat[,"trobit_pct"], VCF_med)

    ## add error bars    
    VCF_low  = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.32)
    VCF_high = apply(TRB_equivilent, 1, find_VCF_pd_point, 0.68)
    arrows(dat[,"trobit_pct"], VCF_low, dat[,"trobit_pct"], VCF_high,
           length=0.05, angle=90, code=3, col = cols)
           
    ## add trend lines
    add_trend_line <- function(type = NULL) {
        
        x = dat[,"trobit_pct"] / 100 # x is our trobit medians
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
    }

    ## run for all and "forest", "savanna"
    add_trend_line()
    add_trend_line("forest")
    add_trend_line("savanna")
    mtext(side = 3, paste("Clumping", vcf_clumping))
}
## add teh legend
graphics.off()
png("figs/tribit_vs_VCF.png", height = 10, width = 7.5, res = 300, units = 'in')
    layout(rbind(2:1, 4:3, 6:5, 8:7))
    par(mar = c(1, 1, 1, 0.5), oma = c(3, 3 , 1, 0))
    lapply(vcf_clumpings, test_clumping)
    axis(1)
    legend('topleft', legend = names(col_choices), col = col_choices, pch = 19)
    mtext("Trobit cover (%)", side = 1, line = 2.5)
    mtext("VCF cover (%)", side = 2, line = 1.5, outer = TRUE)
dev.off()