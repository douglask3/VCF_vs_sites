text.units <- function(x, y, txt, adj = NULL, ...) {
    txt = as.character(txt)
    for (i in unitExpressionSearches) txt = findAndReplace(txt, i)
    
    if(length(txt)==1) text(x,y, txt, adj = adj,...)
    else {
        mtextNsplit(txt, adj = adj,...,  FUN = text, x = x, y = y)
    }
}


mtext.units <- function(txt,...) {
    txt = as.character(txt)
    for (i in unitExpressionSearches) txt = findAndReplace(txt, i)
    if(length(txt)==1) mtext(txt,...)
    else {
        mtextNsplit(txt,...)
    }
}

findAndReplace <- function(txt, exp) {
    fexp = paste('~',exp[1],'~',sep="")
    rexp = exp[2]

    findAndReplacei <- function(t) {
        if (class(t)=="call") return(t)
        if (!grepl(fexp,t)) return(t)
        t = gsub(fexp,'~$$~~##~~$$~',t)
        t = unlist(strsplit(t,'~$$~',TRUE))

        t[t=='~##~'] = rexp
        return(t)
    }
    #if (length(txt)>1) browser()
    return(unlist(sapply(txt,findAndReplacei)))
}

mtextNsplit <- function(txt, line = 0,..., adj = NULL) {
    # sourceAllLibs(); plot(0); mtext.units("Resprouter\nFixed ~CO2~")
    
    if (any(grepl('\n',txt))) {
        txti = txtj = list()
        for (i in txt) {
            if (class(i)=="character" && grepl("\n",i)) {
                j = strsplit(i,"\n")[[1]]
                txti = c(txti, j[[1]])
                txtj = c(txtj, list(txti))
                txti = list(j[[2]])
                if (length(j)>2) browser()
            } else txti= c(txti, i)
        }
        txtj  = c(txtj, list(txti))
        line = line+1.33*(-length(txtj)/2)+1.33*(1:length(txtj))-0.5
        if (is.null(adj)) adj = 0.5
        adj = lapply(1:length(txtj) - length(txtj)/2 - 0.5, function(i) c(0, i) + adj)
        mapply(combineMtext, txtj, adj = adj, line = rev(line), MoreArgs=list(...))
    } else return(combineMtext(txt, line = line, ...))

}

combineMtext <- function(txt, ..., FUN = mtext) {
    
    combineBquaote <- function(a, b) bquote(paste( .(a), .("\n"), .(b)))
    txti=txt[[1]]
    for (i in txt[-1]) txti=combineBquaote(txti,i)
    print(txti)
    FUN(txti, ...)
}

unitExpressionSearches <- list(
    list("km2"    , bquote(km^2)         ),
    list("m2"     , bquote(m^2)          ),
    list("m-2"    , bquote(m^{-2})       ),
    list("yr-1"   , bquote(yr^{-1})      ),
    list("month-1", bquote(month^{-1})   ),
    list("DEG"    , bquote(degree)       ),
    list("alpha"  , bquote(alpha)        ),
    list("omega"  , bquote(omega)        ),
    list("sigma"  , bquote(sigma)        ),
    list("DELTA"  , bquote(Delta)        ),
    list("infinity"  , bquote(infinity)        ),
    list("CO2"    , bquote(CO[2])        ),
    list("AET/PET", bquote(over(AET,PET))),
    list("PET"    , bquote({}[PET])      ),
    list("AET"    , bquote({}^AET)       ),
    list("_max"   , bquote({}[max])      ),
    list("_cmax"   , bquote({}[cmax])      ),
    list("_0"     , bquote({}[0])        ),
    list("_W"     , bquote({}[W])        ),
    list("_omega" , bquote({}[omega])    ),
    list("_ig"    , bquote({}[ig])       ),
    list("_S"     , bquote({}[S])        ),
    list("_ECM"   , bquote({}[ECM])      ),
    list("_pas"   , bquote({}[pas])      ),
    list("_Tree"  , bquote({}[Tree])     ),
    list("_popDenS", bquote({}[popDenS]) ),
    list("_popDenIg", bquote({}[popDenIg])))

