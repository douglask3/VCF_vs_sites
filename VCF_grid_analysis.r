library(gdalUtils)
library(raster)

VCF_dir = 'data/VCF2/'
files = list.files(VCF_dir, pattern = "MOD44B")
temp_file = "temp/dat.tif"
params_files = list.files("outputs/", pattern = "Savanna", full.names = TRUE)

newproj <- "+proj=longlat +datum=WGS84"
dat_global = NULL
mask_global = NULL
grab_cache = TRUE
gridFile <- function(file) {
    print(file)
    temp_file_file = paste0("temp/VCF_stitching_control-", file, ".nc")
    if (file.exists(temp_file_file) && grab_cache) {
        out = brick(temp_file_file)
        return(out)
    } 
    
    sds = get_subdatasets(paste0(VCF_dir, file))
    gdal_translate(sds[1], dst_dataset = temp_file)

    dat = raster(temp_file)/200
    ext = extent(dat)
    ext = ext + c(-1, 1, -1, 1) * c(diff(ext[1:2]), diff(ext[3:4])) * 0.25
    
    gridAndMask <- function(dat, tmask = FALSE) {
        dat = dat0 = extend(dat, ext)        
        if (tmask) mask = dat
        test = is.na(dat)
        dat [test] = 0
        
        aggPrj <- function(r) 
            raster::aggregate(r, fact = 100)

        dat = aggPrj(dat)

        if (tmask){
            mask[test] = 0
            mask[!test] = 1
            mask = aggPrj(mask)
            dat = addLayer(dat, mask)
        }
        
        return(dat)
    }
    out = gridAndMask(dat, TRUE)
    
    writeRaster(out, file = temp_file_file, overwrite = TRUE)
    return(out)
}

out = lapply(files[1:3], gridFile)
browser()

r = raster(crs = newproj, res = 0.2)
r[] = 0
for (ot in out) {
    ri = projectRaster(ot[[1]], crs = newproj)
    ri = resample(ri, r)
    ri[is.na(ri)] = 0
    r = r +ri
}
    

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
