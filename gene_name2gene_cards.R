###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################

gene_name2gene_cards<-function(g){
    #  convert vector of gene_names to HTML links at gene_cards
    return( sapply(g, function(x){ paste0('[', x, '](http://www.genecards.org/cgi-bin/carddisp.pl?gene=', x, ')') }))
}

