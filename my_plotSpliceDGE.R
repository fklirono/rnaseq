###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################




my_plotSpliceDGE<-function(DS, GENEID, GENECOL, FDR=0.05, CND){
    #       DS : DGELRT object
    #   GENEID : identifier for the gene to pick up and plot (gene_id, gene_name, etc.)
    #  GENECOL : column name that contains the GENEID identifiers
    #      FDR : cutoff for FDR-adjusted p-values
    #      CND : differential expression condition string to add to filenames, e.g. 'MNA vs HR_nMNA' which will be converted to 'MNA_HR_nMNA'
    require(data.table)
    require(edgeR)


    #  identify the exons we are going to plot
    ex<-as.data.table(DS$genes[ DS$genes[, GENECOL] %in% GENEID, ])


    #  test exons separately for differential splicing
    #  add FDR to the exons of interest found significantly differentially spliced
    ex.ds<-as.data.table(topSpliceDGE(DS, test='exon', n=Inf, FDR=FDR))
    ex.ds<-ex.ds[ unlist(ex.ds[, GENECOL, with=F]) %in% unlist(ex[, GENECOL, with=F]) ][, c('exon_id', 'FDR'), with=F]


    #  identify if the exons of interest are found individually differentially spliced as well
    #  get the log2-fold-changes of each exon against the average across conditions
    ex<-ex.ds[ ex, on='exon_id' ]  #  this will work even if ex.ds is empty
    ex$logFC<-DS$coefficients[ DS$genes[, GENECOL] %in% GENEID ]
    ex<-ex[ order(start) ]


    #  invert exon order if gene is in minus strand
    if (unique(ex$strand)=='-'){
        ex<-ex[ order(-start) ]
    }


    #  do the plot of log2FC (each exon relative to the average) for each exon and save
    x11(width=18, height=14, title='', bg='white', type='cairo', pointsize=20, antialias='subpixel', family='Arial')
    par(mar=c(12.5,7.0,1.0,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.0, cex.axis=2.0)  #  do the plot with an initial margin
    YLIM<-range(pretty(range(ex$logFC, na.rm=T), 5))
    plot(ex$logFC, type='b', lwd=4, xlab='', ylab='', main='', ylim=YLIM, cex=1.4, xaxt='n')
    BOTTOM<-max(ceiling(log10(ex$start))) + 2
    par(mar=c(BOTTOM,7.0,1.0,0.5), mgp=c(3,1,0), oma=c(0,0,0,0), las=1, xpd=F, bty='n', cex.lab=2.0, cex.axis=2.0)
    options(scipen=0)
    plot(ex$logFC, type='b', lwd=4, xlab='', ylab='', main='', ylim=YLIM, cex=1.4, xaxt='n')
    abline(h=0, lty=3, lwd=4, col='pink')
    mtext(GENEID, side=3, line=0, padj=+0.5, cex=1.8)
    axis(1, at=seq_along(ex$logFC), labels=ex$start, las=2, cex.axis=1.8)
    mtext(expression(log[2]~'fold change'), side=2, line=4, padj=-0.6, cex=2.0, las=3)
    mtext('Exon start', side=1, line=BOTTOM-1, padj=-0.1, cex=2.0, las=1)
    if( ex[, any(!is.na(FDR))]){
        n<-ex[, which( !is.na(FDR) ) ]
        points(n, ex$logFC[n], col='red3', pch=19, cex=ex[n, (-log10(FDR)-min(-log10(FDR)))/diff(range(-log10(FDR)))+1.8])
    }
    dev.print(device=pdf, file=paste0('/data/sequencing/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotDS_', GENEID, '_', sub(' vs ', '_', CND), '.pdf'), width=ceiling(nrow(ex)*0.65), height=14, bg='white', colormodel='cmyk', pointsize=20, useDingbats=F, family='Arial')
    dev.print(device=svg, file=paste0('/data/sequencing/2017-11-08_Fuchs_totalRNAseq/raw/unified/figures/plotDS_', GENEID, '_', sub(' vs ', '_', CND), '.svg'), width=ceiling(nrow(ex)*0.65), height=14, bg='white', antialias='subpixel', pointsize=20, family='Arial')


   invisible()
}

