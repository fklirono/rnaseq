###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




circ2genome<-function(CIR, EXN, SEQ){
    #  function that converts to global genomic coordinates the CIR coordinates along the circRNA sequence that 
    #  is made up by pasting the exons in EXN, in the end the genomic sequences obtained should be identical to 
    #  the seed sequences in SEQ
    #
    #  KNOWN BUG: spliced seeds are not done correctly, that would mean GRangesLists and lists of lists just for 
    #             a few cases...

    WID<-cumsum(width(EXN))
    ST<-start(CIR)
    EN<-end(CIR)
    STRAND<-as.character(strand(EXN[1]))
    CHR<-as.character(seqnames(EXN[1]))

    for(s in seq_along(CIR)){
        if(STRAND=='+'){
            r<-which( (WID-ST[s])>=0 )[1]    #  identify exon rank where the seed start is found
            stopifnot( r>=1 )
            ST[s]<-start(EXN[r]) + ST[s] - 1  
            if (r>1){  #  we need to subtract the sequence part the belongs to previous exons
                ST[s]<-ST[s] - WID[ r - 1 ]
            }
        } else {
            r<-which( (WID-EN[s])>=0 )[1]    #  identify exon rank where the seed start is found
            stopifnot( r>=1 )
            ST[s]<-end(EXN[r]) - EN[s] + 1  
            if (r>1){  #  we need to subtract the sequence part the belongs to previous exons
                ST[s]<-ST[s] + WID[ r - 1 ]
            }
        }
        EN[s]<-ST[s] + width(CIR[s]) - 1
    }

    return(GRanges(seqnames=rep(CHR, length(ST)), strand=rep(STRAND, length(ST)), ranges=IRanges(start=ST, end=EN), DataFrame(seqs=SEQ)))
}
