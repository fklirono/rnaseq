###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################

msigdb2link<-function(g, trim=NULL){
    #  convert vector of MSigDB gene set names to HTML links pointing to the corresponding description
    #  if trim is set then gene set name is trimmed to the given length adding ... at the end
    if(!is.null(trim)){
        require(data.table)
        orig<-g
        n<-which(nchar(g)>trim)
        g[n]<-paste0(substring(g[n], 1, trim-3), '...')
        d<-data.table(original=orig, trimmed=g)
        d[, link:=paste0('[', trimmed, '](http://software.broadinstitute.org/gsea/msigdb/cards/', original, ')')]
        return(d$link)
    } else {
        return( sapply(g, function(x){ paste0('[', x, '](http://software.broadinstitute.org/gsea/msigdb/cards/', x, ')') }))
    }
}

