library(gdalUtils)
library(raster)
library(rasterExtras)
library(snow)

VCF_dir = 'data/VCF2/'
files = list.files(VCF_dir, pattern = "MOD44B")
temp_file = "temp/Tdat"
params_files = list.files("outputs/", pattern = "Savanna", full.names = TRUE)

newproj <- "+proj=longlat +datum=WGS84"

temp_file_base = "temp/VCF_masked4-stitching80-2-"

grab_cache = TRUE

muliCore = F

ens_nos = round(seq(1, 5000, length.out = 1001))
ens_nos = ens_nos[seq(1, length(ens_nos), length.out = 51)]

LC_dir = 'data/Cover/'
LC_files = list.files(LC_dir)

rFinal =  raster(crs = newproj, res = 0.5)
rFinal[] = 0
rFinal = crop(rFinal, c(-180, 180, -30, 30))

trans <- function(VCF, a, b) {
    tr = (1/VCF) - 1
    tr = (tr * exp(b))^(1/a)
    tr = 1/(tr + 1)
    return(tr)
}

gridFile <- function(file) {
    
    temp_file_all = paste0(temp_file_base, '-All-', file,'-',
                           paste0(range(ens_nos), collapse = '_to_'), '--',
                           paste0(unique(diff(ens_nos)), collapse = '_'), ".Rd")
    
    if (file.exists(temp_file_all) && grab_cache) {
        load(temp_file_all)
        return(out)
    }
    #return(NULL)
    #print(file) 
    library(gdalUtils)
    library(raster)

    temp_file_vcf = paste0(temp_file, file, '.tif')
    if (!file.exists(temp_file_vcf)) {        
        sds = get_subdatasets(paste0(vcf_dir, file))
        gdal_translate(sds[1], dst_dataset = temp_file_vcf)
    }
    dat = raster(temp_file_vcf)

    MODISgrid = strsplit(file, '.', fixed = TRUE)[[1]][c(3, 5)]
    LC_file = LC_files[grepl(MODISgrid[1], LC_files)]
     if (length(LC_file) == 0) {
        out = NULL
        save(out, file = temp_file_all)
        return(out)
    }
    if (length(LC_file) > 1) {
        test  = grepl(substr(MODISgrid[2], 1, 6), LC_file)
        if (sum(test) == 0) LC_file = LC_file[1]
        else {
            LC_file = LC_file[test]
            if (length(LC_file) > 0) browser()
        }
    }
    temp_file_lc = paste0(temp_file, LC_file, '.tif')
    if (!file.exists(temp_file_lc)) {        
        sds = get_subdatasets(paste0(LC_dir, LC_file))
        gdal_translate(sds[1], dst_dataset = temp_file_lc)
    }
    lu = raster(temp_file_lc)
    
    if (!any(sapply(unique(lu), function(i) any(i == (6:10))))) {
        out = NULL
        save(out, file = temp_file_all)
        return(out)
    }
    lu = raster::resample(lu, dat, method = "ngb")
    lu_mask = lu >= 6 & lu <= 10
    
    #u_mask = raster::resample(lu_mask, dat, method = "ngb")  
    
    dat[dat>100] = NaN
    dat[is.na(lu_mask)] = NaN
    dat = dat/80    
    
    dat_all = dat
    dat[lu_mask < 0.5] = NaN

    test_na = !is.na(dat)
    lu_v = lu[test_na]

    ext = extent(dat)
    ext = ext + c(-1, 1, -1, 1) * c(diff(ext[1:2]), diff(ext[3:4])) * 0.25
    dat_all = raster::extend(dat_all, ext)
    gridAndMask <- function(dat, tmask = FALSE, params = NULL,
                            pfile = NaN, ens = 0) {
        temp_file_file = paste0(temp_file_base,
                                file,'-', pfile, '-', ens, ".nc")
        temp_file_save =  paste0(temp_file_base,
                                file,'-', pfile, '-', ens, "Rd")
        
        if (file.exists(temp_file_save) && grab_cache) {
            load(temp_file_save)
            return(dat)
        }
        #print("start")
        #start_time = Sys.time()
        if (!is.null(params)) {
            vals = trans(dat[test_na], params[ens, 3], params[ens, 2]) - dat[test_na]
            #dat_i = addLayer(dat_i, dat_i, dat_i)
            #dat_i[[1]][dat_i[[1]] > 0] = 0
            #dat_i[[2]][dat_i[[2]] < 0] = 0
            #print(Sys.time() - start_time)
            posnegLU <- function(ty) {            
                test = lu_v == ty
                vals = vals[test]
                neg = sum(vals[vals<0], na.rm = TRUE)
                pos = sum(vals[vals>0], na.rm = TRUE)
                
                return(c(neg, pos, sum(test)))
            }
            negPos = sapply(6:10, posnegLU)            
            
            #print(Sys.time() - start_time)
            dat[test_na] = vals
           
        }
        #print(Sys.time() - start_time)
        dat0 = dat
        dat = raster::extend(dat, ext)  
          
        if (tmask) mask = addLayer(dat, dat)      
        
        aggPrj <- function(rin) {            
            #print(Sys.time() - start_time)
            #r0 = rin
            rin = raster::aggregate(rin, fact = 200)
            rin = projectRaster(rin, crs = newproj)
            test = class(try(intersect(rin, rFinal),
                             silent = TRUE))
            
            if (test == "try-error") return(NULL)
            rin =  raster::resample(rin, rFinal) 
            gc()
            return(rin)          
        }
        #print(Sys.time() - start_time)
        nl = nlayers(dat)
         
        if (nl>1) {
            test = is.na(dat[[1]])
            for (i in nl) dat[[i]][test] = 0
            dat0 = dat
            dat = layer.apply(dat, aggPrj)
        } else {
            test = is.na(dat)
            dat0 = dat
            dat[test] = 0            
            dat = aggPrj(dat)            
        }

        #print(Sys.time() - start_time)

        if (!is.null(dat)) {
            if (tmask){               
                mask[[1]][ test ] = 0
                mask[[1]][!test ] = 1
                
                test2 = is.na(dat_all)
                mask[[2]][ test2] = 0
                mask[[2]][!test2] = 1 
                
                mask = layer.apply(mask, aggPrj)
                dat = addLayer(dat, mask)
            }     
          
            #print(Sys.time() - start_time)
            dat = writeRaster(dat, file = temp_file_file,
                        overwrite = TRUE)  
            if (exists("negPos")) dat = list(dat, negPos)
            save(dat, file = temp_file_save)
            #print(Sys.time() - start_time)  
        }       
        return(dat)
    }
    control = gridAndMask(dat, TRUE)   
    if (is.null(control)) {
        out = NULL
        save(out, file = temp_file_all)
        return(out)
    }
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

if (muliCore) {
    print("yay")

    cl = makeSOCKcluster(rep("localhost", 3))
        clusterExport(cl = cl, list("ens_nos", "grab_cache",
                                    "temp_file","VCF_dir", "LC_files", "LC_dir",
                                    "params_files", "newproj","temp_file_base",
                                    "layer.apply", "trans", "rFinal"))
        
        outs = parLapply(cl, files[1:100], gridFile)
    stopCluster(cl)
} else {
    outs = lapply(files[101:200], gridFile)
}

outs = outs[!sapply(outs, is.null)]
#browser()
rDiffs = addLayer(rFinal, rFinal, rFinal)
compare4Mask <- function(mi = 1) {
    mask = sum(layer.apply(outs, function(i) i[[1]][[mi+1]]), na.rm = TRUE)
    control = sum(layer.apply(outs, function(i) i[[1]][[1]]), na.rm = TRUE) / mask
    
    
    experiment <- function(id) {
        library(raster)
        library(rasterExtras)
        out = lapply(outs, function(i) i[[id+1]])
     
        ensMember <- function(i, index = 1) {
            temp_file_ens = paste0(temp_file_base,"-6mask_", mi, 
                                   '-', id, "-", i, ".nc")            
            print(temp_file_ens)
            
            if (file.exists(temp_file_ens)) return(brick(temp_file_ens))            
            ensR = sum(layer.apply(out, function(ri) ri[[i]][[index]]), na.rm = TRUE) / mask            
            ensR = writeRaster(ensR, file = temp_file_ens, overwrite = TRUE)            
            return(ensR)
        }            
        enss = layer.apply(1:length(out[[1]]), ensMember)
        
        qenss = apply(enss[], 1, quantile, c(0.05, 0.5, 0.95),
                      na.rm = TRUE)
        rDiffs[] = t(qenss)

        lc_assess <- function(out) #{       
            out[[2]][1:2,]
        #    ltc = ltc0 = out[[2]]            
            #if (any(ltc[1:2,] > 1)) browser()
            #ltc = sweep(ltc[1:2,],2, ltc[3,], '*')
        #}    
        ltcs = lapply(out, lapply, lc_assess)
        ens_change = lapply(1:length(ltcs[[1]]), function(eno)
                            Reduce('+', lapply(ltcs, function(i) i[[eno]])))
        if (mi == 1) {
            area = lapply(out, function(i) i[[2]][[2]][3,])
            area = Reduce('+', area)
        }
        ens_change = array(unlist(ens_change), dim = c(2,5, 6))
        qchange = apply(ens_change, 1,
                        function(i) list(apply(i, 1, quantile, 
                        c(0.05, 0.5, 0.95))))

        
        return(c(rDiffs, qchange,area))
    }
    
    if (F) {
        cl = makeSOCKcluster(rep("localhost", 3))
            clusterExport(cl = cl, list("outs", "rDiffs", "temp_file_base"))
        
            exps = parLapply(cl, 1:4, experiment)
        stopCluster(cl)
    } else {
        exps = lapply(1:4, experiment)
    }    
    save(control, exps, mask,
         file = paste0("outputs/gridded_VCF_correction-6mask", mi, ".Rd"))
}

compare4Mask(2)
compare4Mask(1) 
