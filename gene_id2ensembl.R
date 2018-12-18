###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################

gene_id2ensembl<-function(g){
    #  convert vector of gene_ids to HTML links at Ensembl browser
    return( sapply(g, function(x){ paste0('[', x, '](http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=', x, ')') }))
}

