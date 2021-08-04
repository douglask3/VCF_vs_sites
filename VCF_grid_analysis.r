library(gdalUtils)
library(raster)
library(rasterExtras)
library(snow)
source("libs/logit_logistic.r")

VCF_dir = 'data/VCF2/'
files = list.files(VCF_dir, pattern = "MOD44B")
temp_file = "temp/Tdat"
params_files = list.files("outputs/", pattern = "trsnlateCurve1", full.names = TRUE)

newproj <- "+proj=longlat +datum=WGS84"

temp_file_base = "temp/VCF_masked-stan-full-stitching80-newsiteInfo-maxMinMed3-"

grab_cache = TRUE
muliCore = F
ncores = 3

#ens_nos = round(seq(1, 5000, length.out = 1001))
#ens_nos = seq(1, 1000, length.out = 51)
#ens_nos = ens_nos[1:5]

LC_dir = 'data/Cover/'
LC_files = list.files(LC_dir)

rFinal =  raster(crs = newproj, res = 0.5)
rFinal[] = 0
rFinal = crop(rFinal, c(-180, 180, -30, 30))

LU_keys = 1:10

trans <- function(VCF, params) {

    findCloseReplace <- function(y) 
        params[which.min(abs(y - params[,1])), 2]
    out = sapply(VCF, findCloseReplace)
    return(out)
    VCF0 = VCF
    
    #VCF = VCF0
    VCF = params[['VCF0']] + logit(VCF)*params[['VCFD']]
    #VCF01 = logistic(VCF)
    VCF = VCF1 = (VCF-params[['alpha']])/params[['beta']]
    VCF = exp(VCF)

    x = seq(0.00, 1, 0.005)
    f1 = x^params[['tau1']]
    x2 = (1-x)^params[['tau2']]

    transBack <- function(y) {
         if( is.infinite(y)) return(1)
        x[which.min(abs(y*x2-f1))]
    }
    VCF = sapply(VCF, transBack)
    VCF = logistic((logit(VCF) - params[['VCF0']])/params[['VCFD']])
    
    return(VCF)
}

gridFile <- function(file) {
    
    temp_file_all = paste0(temp_file_base, '-All-', file,'-', ".Rd")
                           #paste0(range(ens_nos), collapse = '_to_'), '--',
                           #paste0(unique(diff(ens_nos)), collapse = '_'), ".Rd")
    print(file)
    if (file.exists(temp_file_all) && grab_cache) {
        load(temp_file_all)
        if (is.null(out) || (
            length(out) >1 && is.raster(out[[2]][[1]][[1]][[1]][[1]]))) return(out)
    }
    
    #return(NULL)
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
    
    if (!any(sapply(unique(lu), function(i) any(i == (LU_keys))))) {
        out = NULL
        save(out, file = temp_file_all)
        return(out)
    }
    
    lu = raster::resample(lu, dat, method = "ngb")
    lu_mask = lu >= min(LU_keys) & lu <= max(LU_keys)
    
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

    LUMasks = lapply(LU_keys,  function(ty) lu_v == ty)
    
    vals = dat[test_na]
    vals[vals >1] = 1
    vals[vals <0] = 0
    vals0 = vals
    gridAndMask <- function(dat, tmask = FALSE, params = NULL,
                            pfile = NaN, ens = 0) {
        
        temp_file_file = paste0(temp_file_base,
                                file,'-', pfile, '-', ens, ".nc")
        temp_file_save =  paste0(temp_file_base,
                                file,'-', pfile, '-', ens, ".Rd")
        
        if (file.exists(temp_file_save) && grab_cache) {
            load(temp_file_save)            
            if (is.raster(dat[[1]]) || is.raster(dat[[1]][[1]]))return(dat)
        }
        print(temp_file_save)
        
        if (!is.null(params)) {            
            vals  = trans(vals0, params[,1+c(1, ens)]) 
            valsi = vals - vals0
            posnegLU <- function(ty) {
                valsi = valsi[ty]
                neg = sum(valsi[valsi<0], na.rm = TRUE)
                pos = sum(valsi[valsi>0], na.rm = TRUE)
                
                return(c(neg, pos, sum(ty)))
            }
            negPos = sapply(LUMasks, posnegLU)        
            dat[test_na] = valsi
        }
        histLU <- function(ty) 
            hist(vals[ty], plot = FALSE, breaks = seq(0, 1, 0.01))[c("counts","mids")]
        
        hists = lapply(LUMasks, histLU)
        
        dat0 = dat
        dat = raster::extend(dat, ext)  
          
        if (tmask) mask = addLayer(dat, dat)      
        
        aggPrj <- function(rin) {          
            rin = raster::aggregate(rin, fact = 200)
            rin = projectRaster(rin, crs = newproj)
            test = class(try(intersect(rin, rFinal),
                             silent = TRUE))
            
            if (test == "try-error") {
                print(paste(c("no overlap. Extent:" , extent(rin)[1:4]), collapse = ' '))
                return(NULL)
            }
            rin =  raster::resample(rin, rFinal) 
            gc()
            return(rin)          
        }
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
          
            dat = writeRaster(dat, file = temp_file_file,
                        overwrite = TRUE)  
            if (exists("negPos")) dat = list(dat, negPos)
            dat = list(dat, hists)
            
            save(dat, file = temp_file_save)
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
        print(pfile)
        params = read.csv(pfile)/100
        
        name = strsplit(pfile, '#')[[1]][2]
        name = strsplit(name, '.csv')[[1]][1]
        
        out = lapply(2:4, function(i)
                        gridAndMask(dat, tmask = FALSE, params,
                                    pfile = name, ens = i))
        
    }
    
     
    out = lapply(params_files, forParams)
    
    out = list(control, out)
    save(out, file = temp_file_all)
    return(out)
}

if (muliCore) {
    print("yay")

    cl = makeSOCKcluster(rep("localhost", ncores))
        clusterExport(cl = cl, list("ens_nos", "grab_cache",
                                    "temp_file","VCF_dir", "LC_files", "LC_dir",
                                    "params_files", "newproj","temp_file_base","LU_keys",
                                    "layer.apply", "trans", "rFinal", "logit", "logistic"))
        
        outs = parLapply(cl, files[1:30], gridFile)
    stopCluster(cl)
} else {
    outs = lapply(files, gridFile) #2,6, [seq(1, length(files), by = 1)]
}

outs = outs[!sapply(outs, is.null)]

rDiffs = addLayer(rFinal, rFinal, rFinal)

compare4Mask <- function(mi = 1) {
    
    mask = sum(layer.apply(outs, function(i) i[[1]][[1]][[mi+1]]), na.rm = TRUE)
    control = sum(layer.apply(outs, function(i) i[[1]][[1]][[1]]), na.rm = TRUE) / mask   
    
    cnt = lapply(outs, function(i) i[[1]][[2]])
    cnt = lapply(1:length(cnt[[1]]), function(id)
                    Reduce('+', lapply(cnt, function(h) h[[id]][[1]])))
    
    experiment <- function(id) {
        library(raster)
        library(rasterExtras)
        
        hists = lapply(outs, function(i) list(i[[2]][[id]][[2]][[2]])[[1]][[1]])
        hists = lapply(1:length(hists[[1]]), function(ty)
                        Reduce('+', lapply(hists, function(i) i[[ty]][[1]]) ))
        out = lapply(outs, function(i) i[[2]][[id]])
        
        ensMember <- function(i) {
            temp_file_ens = paste0(temp_file_base,"-7mask_", mi, 
                                   '-', id, "-", i, ".nc")            
            print(temp_file_ens)            
            if (file.exists(temp_file_ens) & F) return(brick(temp_file_ens))                     
            
            ensR = sum(layer.apply(out, function(ri) ri[[i]][[1]][[1]]),
                        na.rm = TRUE) / mask      
            
            #ensR = writeRaster(ensR, file = temp_file_ens, overwrite = TRUE)            
            return(ensR)
        }  
       
        enss = layer.apply(1:length(out[[1]]), ensMember)
        #browser()
        #qenss = apply(enss[], 1, quantile, c(0.1, 0.5, 0.9),
        #              na.rm = TRUE)
        rDiffs[] = enss[]#t(qenss)

        lc_assess <- function(out)     
            out[[1]][[2]][1:2,]       
        
        ltcs = lapply(out, lapply, lc_assess)
        
        ens_change = lapply(1:length(ltcs[[1]]), function(eno)
                            Reduce('+', lapply(ltcs, function(i) i[[eno]])))
        if (mi == 1) {             
            area = lapply(out, function(i) i[[1]][[1]][[2]][3,])
            area = unlist(Reduce('+', area))
        } else area = 0
        ens_change = array(unlist(ens_change), dim = c(2,10, 3))
        
        qchange = apply(ens_change, 1,
                        function(i) list(apply(i, 1, quantile, 
                        c(0.0, 0.5, 1.0)))) 
        
        return(c(rDiffs, qchange,c(area), c(cnt), c(hists)))
    }
    
    if (F) {
        cl = makeSOCKcluster(rep("localhost", 4))
            clusterExport(cl = cl, list("outs", "rDiffs", "temp_file_base"))
            
            exps = parLapply(cl, 1:4, experiment)
        stopCluster(cl)
    } else {
        exps = lapply(1:4, experiment)
    }   
     
    save(control, exps, mask,
         file = paste0("outputs/gridded_VCF_correction-stan-test2-mask-minMedMax-", mi, ".Rd"))
}

compare4Mask(2)
compare4Mask(1) 
