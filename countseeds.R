###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




countseeds<-function(SEEDS, SEQS){
    #  this function counts the occurrence of miRNA seeds (SEEDS) in sequences (SEQS) even when the seeds of a 
    #  given miRNA match overlapping positions.
    #
    #  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #  !                                                                              !
    #  !  Generally one would want to count overlapping seed positions as one match,  !
    #  !  using the following python script:                                          !
    #  !                                                                              !
    #  !      countseeds.py --targets targets.fa --seeds seeds.fa counts.tsv          !
    #  !                                                                              !
    #  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #
    #  SEEDS = DNAStringSet of identical length miRNA seeds (reverse-complemented) named after the corresponding miRNA 
    #          different seeds with same miRNA name are allowed
    #  SEQS = DNAStringSet of named sequences where seed occurrences will be counted on
    require(data.table)
    require(Biostrings)


    #  convert to strings since rownames will be defined by them and matchPattern() functions do not allow for duplicate names
    M<-data.table(data.frame(seed=unname(sapply(SEEDS, toString)), mir=names(SEEDS)))


    #  count seed occurrences
    M.pd<-PDict(M$seed)
    v<-t(as.data.frame(vcountPDict(M.pd, SEQS, max.mismatch=0, with.indels=F, fixed=T)))
    rownames(v)<-names(SEQS)
    colnames(v)<-M$seed


    #  collect transcript names and seed counts per seed 
    v<-apply(v, 2, function(x){ y<-x[x!=0];  data.table(tx_name=list(names(y)), ncount=list(unname(y))) } )
    stopifnot( all.equal( names(v), M$seed ) )
    v<-do.call(rbind, v)
    v$seed<-M$seed
    v$mir<-M$mir


    #  uncollapse transcript names and counts to have unique entries for each transcript:seed interaction
    #  N.B. this eliminates empty interactions without a transcript or count
    v<-v[, .(tx_name=unlist(tx_name), ncount=unlist(ncount)), by=.(seed, mir)]


    #  sum over seed counts for given miRNA:transcript interaction
    v<-v[, .(seed=list(seed), ncount=sum(ncount)), by=.(mir, tx_name)]


    #  reorder the columns
    v<-v[, c('mir', 'tx_name', 'seed', 'ncount'), with=F]


    return(v)
}

