sd.raster <- function(ldata,pmean=TRUE) {
    llayers=nlayers(ldata)
    ones=rep(1,llayers)
    lmean=stackApply(ldata,ones,'mean',na.rm=TRUE)
    ldelt=ldata-lmean

    ldelt=ldelt*ldelt

    lvarn=stackApply(ldelt,ones,'sum',na.rm=TRUE)/llayers

    lvarn=sqrt(lvarn)
    if (pmean) lvarn=lvarn/abs(lmean)
    return(lvarn)

}
