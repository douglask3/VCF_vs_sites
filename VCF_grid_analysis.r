library(gdalUtils)
library(raster)

VCF_dir = 'data/VCF/'
files = list.files(VCF_dir, pattern = "MOD44B")
temp_file = "temp/dat.tif"
params_files = list.files("outputs/", pattern = "Savanna", full.names = TRUE)

newproj <- "+proj=longlat +datum=WGS84"
dat_global = NULL
mask_global = NULL
grab_cache = TRUE
gridFile <- function(file) {
    temp_file_file = paste0("temp/VCF_stitching_control-", file, ".Rd")
    if (file.exists(temp_file_file) && grab_cache) {
        load(temp_file_file)
    } else {      
    
        sds = get_subdatasets(paste0(VCF_dir, file))
        gdal_translate(sds[1], dst_dataset = temp_file)

        dat = raster(temp_file)/200
        ext = extent(dat)
        ext = ext + c(-1, 1, -1, 1) * c(diff(ext[1:2]), diff(ext[3:4])) * 0.25
        dat = dat0 = extend(dat, ext)
        
        mask = dat
        test = is.na(dat)
        dat [test] = 0
        mask[test] = 0
        mask[!test] = 1
        dat = projectRaster(dat , crs=newproj, res = 0.01)
        
        if (is.null(dat_global)) {
            dat_global = dat
            mask_global = mask
        } else {
            ext_global = ext_global0 = extent(dat_global)
            ext_dat = extent(dat)
            if (!all(ext_global == ext_dat)) {
                if (ext_global[1] > ext_dat[1]) ext_global[1] = ext_dat[1]
                if (ext_global[2] < ext_dat[2]) ext_global[2] = ext_dat[2]
                if (ext_global[3] > ext_dat[3]) ext_global[3] = ext_dat[3]
                if (ext_global[4] < ext_dat[4]) ext_global[4] = ext_dat[4]
                if (!all(ext_global == ext_global0)) {
                    dat_global = extend(dat_global, ext_global)
                    test = is.na(dat_global)
                    dat_global[test] = 0
                    mask_global = extend(mask_global, ext_global)
                    mask_global[test] = 0
                }
                dat = extend(dat, ext_global)
            }
       
            mask = dat
            test = is.na(dat)
            mask[!test] = 1
            mask[test] = 0
            dat[test] = 0
            dat_global = dat_global + dat
            mask_global = mask_global + dat
        }
        save(dat_global, mask_global, file = temp_file_file)
    }
    dat_global <<- dat_global
    mask_global <<- mask_global  
    print(file)
    browser()
}
lapply(files, gridFile)


trans <- function(VCF, a, b) {
    tr = (1/VCF) - 1
    tr = (tr * exp(b))^(1/a)
    tr = 1/(tr + 1)
    return(tr)
}

mask = dat
mask[!is.na(mask)] = 1
mask[ is.na(mask)] = 0
browser()
#dat = trans(dat, -0.03578928, b = -2.56274018)

#

low_mask = projectRaster(mask, crs=newproj, res = 0.5)
dat      = projectRaster(dat , crs=newproj, res = 0.5)/mask
