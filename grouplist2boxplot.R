###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




grouplist2boxplot<-function(L, L.COL, YLAB='', XLAB.COL='black', YTEXT.LINE=3, ...){
    #  function that takes a named L of different groups that each contain the same number of rows for specific members, i.e. genes
    #  but not necessarily the same number of column entires across the different groups, i.e. number of samples per group can vary, 
    #  and creates a boxplot where member distributions across groups are next to each other. 
    #
    #  For example, suppose L contains the five groups ST4S, LR, IMR, HR_nMNA, MNA and suppose each group contains data.frames with two 
    #  rows for geneA and geneB but the column numbers are different, let's say, ST4S has 3 samples, LR, 2, samples, IMR 4 samples, HR_nMNA 5 samples
    #  and MNA 5 samples. Then this function will produce a boxplot where the distributions will be ordered as follows:
    #
    #      ST4S, LR, IMR, HR_nMNA, MNA   ST4S, LR, IMR, HR_nMNA, MNA               
    #               gene A                        geneB
    #
       #          L : named list of different groups with data.frame elements with rows the members and columns the samples
    #  GROUP.COLORS : colors to use for the boxplots for the different groups
    #      XLAB.COL : color vector for the x-axis labels
    #    YTEXT.LINE : mtext(..., side=2, line=YTEXT.LINE, ..) 
    #           ... : parameters to pass down to par()
    require(data.table)


    #  remove empty groups
    L<-L[ lengths(L)>0 ]


    #  stop if not all list elements have the same number of rows
    stopifnot( all(sapply(L, nrow)==nrow(L[[1]])) )

    
    #  number of columns for each list element
    NC<-sapply(L, ncol)


    #  x-axis names to use
    NAMES<-rownames(L[[1]])


    #  column-bind all groups and name the columns by the group name
    L<-do.call(cbind, L)
    colnames(L)<-rep(names(NC), NC)


    #  for each member, i.e. gene, split the row by the group respecting the original group order
    L<-apply(L, 1, function(y){ split(y, factor(colnames(L), levels=names(NC))) })


    #  unlist only the upper level so that we end up with a list of size (number of genes) * (number of groups)
    L<-unlist(L, recursive=F, use.names=T)


    #  do the boxplot now
    options(scipen=0)
    par(...)
    YTICK<-pretty(c(0, max(floor(sapply(L, max, na.rm=T)))), 4)
    plot(0:1, 0:1, xlim=c(0, length(L))+c(0, 1), type='n', ylim=range(YTICK), axes=F, ann=F, xaxs='i')
    bp<-boxplot(L, col=L.COL, ylab='', xlab='', show.names=F, frame.plot=F, medcol='lightgrey', boxwex=0.8, xpd=F, outline=F, boxcol=L.COL, range=0, add=T)
    mtext(YLAB, side=2, line=YTEXT.LINE, padj=-0.5, las=0, cex=2.0)
    mtext(text=NAMES, side=1, line=0, at=seq(median(seq_along(NC)), length(L), length(NC)), las=2, adj=1.09, cex=2.0, col=XLAB.COL)
    Y0=par('usr')[3]*0.5
    segments(x0=seq(1, length(L), length(NC))-0.4, x1=seq(length(NC), length(L), length(NC))+0.4, y0=Y0, lwd=1, col='black')
    a<-seq(length(NC)+1, length(L), 2*length(NC))-0.4
    b<-seq(2*length(NC), length(L), 2*length(NC))+0.4
    M<-matrix(c(a,b,b,a), nrow=4, byrow=T)
    for(n in seq_len(ncol(M))){
        polygon(x=M[, n], y=c(Y0, Y0, max(YTICK), max(YTICK)), lty=0, col=adjustcolor('black', alpha.f=0.1))
    }
    legend('topright', legend=names(NC), col=L.COL, bty='n', lty=1, lwd=10, pch=NA, cex=1.8, y.intersp=0.8, xpd=T)


    invisible()
}

