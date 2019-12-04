############
## setup  ##
############
library(reldist)

col_choices  = c("savanna" = "#d95f02", "forest" = '#1b9e77')
detail = 0.25

VCF_grid_size = round(c(250,250) * detail) ## size of a vcf pixel
TRB_grid_size = round(c(100,100) * detail) ## size of a trobit pixel
n_bootstraps_grid_size  = 2  # number of times we'll test the uncertanty
pc_test_width = 0.04
vcf_clumping = 0 # How "clumped" are the trees. 0 = no clumping

grabe_cache = TRUE

## some required functions
logit <- function(x) log(x/(1-x))
logistic <- function(x) 1/(1+exp(-x))

is_p_star <- function(P) {
    if (P > 0.1) out = ' '
    else if (P > 0.05) out = '.'
    else if (P > 0.01) out = '*'
    else if (P > 0.001) out = '**'
    else out = '***'
    out
}

make.transparent <- function(col, transparency) {
     ## Semitransparent colours
     tmp <- col2rgb(col)/255
     rgb(tmp[1,], tmp[2,], tmp[3,], alpha=1-transparency)
}


#############
## process ##
#############
## open data
dat = read.csv( 'data/trobit_vcf_comparison.csv')
dat[, 'mvcf_pct'] = dat[,'mvcf_pct'] / 0.8


## finds the probablity density of % cover of a VCF size grid based on a trobit measurement
VCF_grid = matrix(0,VCF_grid_size[1], VCF_grid_size[1])
VCF_no_cells = VCF_grid_size[1] * VCF_grid_size[2]
prob = (1:VCF_no_cells)^vcf_clumping

covert_from_VCF_grid <- function(vc_pc, grid_size) {
    print(vc_pc)
       
    nc = round(vc_pc * VCF_no_cells)
    
    testY <- function(pc) {
        index = sample(1:VCF_no_cells, nc, FALSE, prob = prob)
        VCF_grid[index] = 1

        testXpos <- function(i, j) {
            mn = mean(VCF_grid[i:(i+ grid_size[1]),j:(j+ grid_size[1])])
            abs(mn - pc) < 0.1
        }
        
        testYpos <- function(i) sapply(1:(VCF_grid_size[2] - grid_size[2]), function(j) testXpos(i,j))
        mean(sapply(1:(VCF_grid_size[1] - grid_size[1]), testYpos))
    }

    temp_file = paste('temp/', 'climping', vcf_clumping, 'vcf_pc', vc_pc, 'pc_test_width', pc_test_width, 'nboots', 
                      n_bootstraps_grid_size, 'VCF_grid', VCF_grid_size[1], VCF_grid_size[2], 
                      'TRB_grid', TRB_grid_size[1], TRB_grid_size[2], '.csv', sep = '-')
                      
    pcs = rep(seq(0, 1, pc_test_width), n_bootstraps_grid_size)
    
    if(file.exists(temp_file) && grabe_cache) {
        Ys = read.csv(temp_file)[,1]
    } else {
        
        Ys = sapply(pcs, testY)
        write.csv(Ys, temp_file, row.names = FALSE)
    }
    #c(wtd.mean(pcs, Ys) + c(-1, 0, 1) * sqrt(wtd.var(pcs, Ys)))
    wtd.quantile(pcs, q = c(0.32, 0.5, 0.68), weight = Ys)
}
TRB_equivilent = t(sapply(dat[, 'mvcf_pct']/100, covert_from_VCF_grid, TRB_grid_size)*100)

## for each vcf %, we get back the 0.1, 0.5, and 0.9 of the prob density function.
## 0.5 will be the point, while 0.1 to 0.9 will be our error bars
colnames(TRB_equivilent) = c('trob_vcf_10', 'trob_vcf_50', 'trob_vcf_90')

## tag on vcf % and forest type. We then have all the info we need for plotting
TRB_equivilent = cbind(TRB_equivilent, dat[, c("trobit_pct", 'forest_type')])


##########
## plot ##
##########
cols = col_choices[(TRB_equivilent[, "forest_type"] == "forest")+1]

## set up 
plot(c(0, 100), c(0, 100), xlab = "Trobit cover (%)", ylab = "VCF cover (%)", type = 'n')
grid()
lines(c(0, 100), c(0, 100), lty = 2, lwd = 2)
#text(0, c(95, 90), c('k:', 'x0:'), adj = 0)

## add points
points(TRB_equivilent[,4], TRB_equivilent[,2], pch = 19, col = cols)
points(TRB_equivilent[,4], TRB_equivilent[,2])

## add error bars
arrows(TRB_equivilent[,4], TRB_equivilent[,1], TRB_equivilent[,4], TRB_equivilent[,3],
       length=0.05, angle=90, code=3, col = cols)
       
## add trend lines
add_trend_line <- function(type = NULL) {
    
    x = TRB_equivilent[,4] / 100 # x is our trobit medians
    y = TRB_equivilent[,2] / 100 # y is our vcf value
    y.low =  TRB_equivilent[,1]/100; y.high =  TRB_equivilent[,3]/100

    # if we're looking at just savanna or forest, we'll just grab those bits of data
    if (!is.null(type)) { 
        test = TRB_equivilent[,5] == type
        x = x[test]; y = y[test];  y.low = y.low[test]; y.high = y.high[test]
        col = col_choices[type]
    } else col = 'black'
    
    # cos our domain is bounded [0,1], we'll transform our data (and hence our domain) first, and we'll transform it back when we're done
    xf = logit(x); yf = logit(y); ye = logit(y.high) - logit(y.low)
    
    # do the regression
    ## When we start to consider tree clumping, we may need to start weighting points
    ## by the error bars. However, for now, all VCF/TROBUT uncertanty should be
    ## equal for all sites under transformation, so does not need considering
    fit = lm(yf ~ xf)
    
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

## add teh legend
legend('topleft', legend = names(col_choices), col = col_choices, pch = 19)