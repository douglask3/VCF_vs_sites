######################################
## set up the stuff we want to grab ##
######################################
library("MODISTools")
latlon = read.csv('data/nordeste_locations.csv')

start_date = "2014-01-01"
end_date   = "2019-12-31"
km_ab      = 0.0
km_lr      = 0.0
band       = "Percent_Tree_Cover"

#################################################
## the function we need to download stuff with ##
#################################################
download_site = function(latlon) {
    out = mt_subset(product   = "MOD44B",
              lon       = as.numeric(latlon[3]),
              lat       = as.numeric(latlon[2]),
              band      = band,
              start     = start_date,
              end       = end_date,
              km_lr     = km_lr,
              km_ab     = km_ab,
              site_name = latlon[1],
              internal  = TRUE,
              progress  = FALSE)
    browser()
}
              
box =  apply(latlon[,c("lat", "long")], 2, range) 
lats  = do.call(seq, as.list(c(box[,"lat" ], 0.001)))
longs = do.call(seq, as.list(c(box[,"long"], 0.001)))

xy = cbind(rep(lats, each = length(longs)), longs)


#########################
## download and output ##
######################### 
## dowload 
out = apply(latlon, 1, download_site)   
out = do.call(rbind, out)

## output everything
write.csv(out, 'VCF_data_full.csv')  

## output the fields we're mainly interested in
summ_out = out[c("site", "latitude", "longitude", "calendar_date", "value")]
write.csv(summ_out, 'VCF_data_summary.csv') 


## calculate and output annual averge
siteMean <- function(site) {
    value = mean(site[,"value"])
    c( site[1, 1:3], value = value)
}

site_out = split(summ_out, summ_out[,1], drop = TRUE)
site_out = t(sapply(site_out, siteMean))
write.csv(site_out, 'VCF_data_annual_means.csv') 