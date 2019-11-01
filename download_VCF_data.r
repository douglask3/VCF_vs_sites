######################################
## set up the stuff we want to grab ##
######################################
library("MODISTools")
latlon = read.csv('data/trobit_cai.csv')
fname_out = 'VCF_data_trobit_sites'

start_date = "2014-01-01"
end_date   = "2019-12-31"
km_ab      = 0.0
km_lr      = 0.0
band       = "Percent_Tree_Cover"

#################################################
## the function we need to download stuff with ##
#################################################
download_site = function(latlon)
    mt_subset(product   = "MOD44B",
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
              
  
#########################
## download and output ##
######################### 
## dowload 
out = apply(latlon, 1, download_site)   
out = do.call(rbind, out)

## output everything
write.csv(out, paste0(fname_out, '_full.csv'))

## output the fields we're mainly interested in
summ_out = out[c("site", "latitude", "longitude", "calendar_date", "value")]
write.csv(summ_out, paste0(fname_out, '_summary.csv') )


## calculate and output annual averge
siteMean <- function(site) {
    value = mean(site[,"value"])
    c( site[1, 1:3], value = value)
}

site_out = split(summ_out, summ_out[,1], drop = TRUE)
site_out = t(sapply(site_out, siteMean))
write.csv(site_out, paste0(fname_out, '_means.csv'))