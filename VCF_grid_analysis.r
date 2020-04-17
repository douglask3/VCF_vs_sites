library(gdalUtils)
library(raster)

VCF_dir = 'data/VCF2/'
files = list.files(VCF_dir, pattern = "MOD44B")
temp_file = "temp/dat.tif"
params_files = list.files("outputs/", pattern = "Savanna", full.names = TRUE)

newproj <- "+proj=longlat +datum=WGS84"
dat_global = NULL
mask_global = NULL
grab_cache = FALSE

ens_nos = seq(1, 5000, by = 2000)

trans <- function(VCF, a, b) {
    tr = (1/VCF) - 1
    tr = (tr * exp(b))^(1/a)
    tr = 1/(tr + 1)
    return(tr)
}

gridFile <- function(file) {
    print(file) 
    temp_file_all = paste0("temp/VCF_stitching2-ALL-", file,'-',
                           paste0(ens_nos, collapse = '_'), ".Rd")
    
    if (file.exists(temp_file_all) && grab_cache) {
        load(temp_file_all)
        return(out)
    }
    sds = get_subdatasets(paste0(VCF_dir, file))
    gdal_translate(sds[1], dst_dataset = temp_file)

    dat = raster(temp_file)
    dat[dat>100] = NaN
    dat = dat/100
    ext = extent(dat)
    ext = ext + c(-1, 1, -1, 1) * c(diff(ext[1:2]), diff(ext[3:4])) * 0.25
    
    gridAndMask <- function(dat, tmask = FALSE, params = NULL,
                            pfile = NaN, ens = 0) {
        temp_file_file = paste0("temp/VCF_stitching2-",
                                file,'-', pfile, '-', ens, ".nc")
        
        if (file.exists(temp_file_file) && grab_cache) {
            out = brick(temp_file_file)
            return(out)
        }
        if (!is.null(params))
            dat = trans(dat, params[ens, 3], params[ens, 2])
        dat = dat0 = extend(dat, ext)        
        if (tmask) mask = dat
        test = is.na(dat)
        dat [test] = 0
        
        aggPrj <- function(r) {
            r = raster::aggregate(r, fact = 40)
            projectRaster(r, crs = newproj)
        }
        print("yay")
        dat = aggPrj(dat)

        if (tmask){
            mask[test] = 0
            mask[!test] = 1
            mask = aggPrj(mask)
            dat = addLayer(dat, mask)
        }     
          
        dat = writeRaster(dat, file = temp_file_file,
                    overwrite = TRUE)    
        return(dat)
    }
    control = gridAndMask(dat, TRUE)   
    
    forParams <- function(pfile) {
        params = read.csv(pfile)
        name = strsplit(pfile, '-#')[[1]][2]
        name = strsplit(name, '.csv.')[[1]][1]
        out = lapply(ens_nos, function(i)
                        gridAndMask(dat, tmask = FALSE, params,
                                    pfile = name, ens = i))
    }
    out = lapply(params_files, forParams)
    out = c(control, out)
    save(out, file = temp_file_all)
    return(out)
}

out = lapply(files, gridFile)
browser()
exts = sapply(out, function(i) extent(out[[1]][[1]])[])

r = raster(crs = newproj, res = 0.2)

r[] = 0
for (ot in out) {
    ri = projectRaster(ot[[1]], crs = newproj)
    ri = resample(ri, r)
    ri[is.na(ri)] = 0
    r = r +ri
}
    



mask = dat
mask[!is.na(mask)] = 1
mask[ is.na(mask)] = 0
browser()
#dat = trans(dat, -0.03578928, b = -2.56274018)

#

low_mask = projectRaster(mask, crs=newproj, res = 0.5)
dat      = projectRaster(dat , crs=newproj, res = 0.5)/mask
