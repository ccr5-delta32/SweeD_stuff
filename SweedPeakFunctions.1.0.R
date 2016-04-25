## Bjorn Pieper. MPIPZ Cologne. April 2016. Version 1
## These are several functions used to analyse SweeD output including:
## loadSweeD - returns a single data structure from independent SweeD
##             outputs per chromosome. The path to the SweeD_Reports needs to
##             beprovided with the argument pathand defaults to the current  
##             path. As it stands it expects immediately after the standard 
##             'SweeD_report.' the word 'chrn' where n indicates the chromosome
##             number. The total number of chromo- somes can be provided by the 
##             argument 'nchr', which defaults to 8. A common part of the file-
##             name can be provided with  'common'. A unique population part of 
##             the filename can be provided with 'pop'. And finally, an additio-
##             nal extension can be provided with 'ext'. So the final input 
##             filename is a concatenation of: 
##                  'SweeD_report.' + 'chrn' + 'common' + 'pop' + 'ext' 
##             Depending on SweeD output format the number of lines to skip 
##             before expecting the data can be provided with the argument 
##             'skip'. 
## splits - takes the output of loadSweeD as input and returns the cummu-
##          lative positions that delineate the different chromosomes.
## sweedPeaks - takes the output of loadSweeD as input and detects peaks in 
##              the CLR profile w.r.t. the local context. It returns a vector 
##              containing their indices in the input data. Peak detection is 
##              on loess smoothing, which can be provided with some parameters
##              through the arguments. 'w' best left at 1 actually. Peak 
##              detection is done in batches of CLR estimates and the batch
##              size can be provided with the 'win' argument. The span of the 
##              loess function is automatically derived to be borderline small 
##              but this may not always work well depending on SweeD grid size 
##              etc. As it stands it works like a charm for a batch size of 
##              5000 when the SweeD -grid was chosen to result in 500bp i
##              windows. The smoothing may result in some slight displacement
##              of the detected peak w.r.t. the true highest CLR in the 
##              vicinity so CLR values around the detected peak are evaluated
##              and the highest one is kept. How many CLR values above and 
##              below the detected peak are evaluated can be set with the 
##              argument 'look'. Please note that the loess function may spit 
##              out many warnings but these can be ignored.
## sigPkGenes - takes the SweeD output, the indices of the detected peaks, and
##              a gff table and return ranked list with most significant on top
##              of i) significant peaks, and ii) the genes associated with the 
##              peaks. For now the provided cutoff  simply determines the 
##              quantile above which a peak is considered significant. The 
##              argument 'thr_meth' can be set to 'peak' to determine the 
##              quantiles of the detected peaks only, or to 'all' to determine 
##              the quantiles of all estimated CLR values. The maximum distance 
##              in bp any part of a gene can be from the locationof a peak to 
##              be considered associated with a peak can be provided with the 
##              argument maxPkDist. 
## peakPlot - plots single peaks in the CLR profile while automatically adjust-
##            ing the scale of the plot to contain the entire peak. To choose a
##            peak provide it is best to provide the 'sig.peak' part of the
##            output from the function sigPkGene to the argument 'sigpks'. the
##            first column, 'peak.rank' of the 'peak.gene' part of the output 
##            from sigPkGenes then provides the indices of the peaks to provide
##            to 'pk'. If the output from sigPkGene is not provided it reverts
##            to the full vector of peak indices that needs to be provided to 
##            the argument 'pks'. the original dataframe with CLR results as 
##            returned by 'loadSweeD' needs to be provided to 'clr' and the 
##            gff table to 'gff'. The threshold quantile to be shown on the 
##            plot can be set with 'thr' and whether to take this from peaks
##            only or all CLR values can be set with 'thr_meth' by setting it 
##            to 'peak' or 'all' respectively'. The size of the plotted region 
##            can be changed by setting 'addspan'. The genes in the plotted
##            can be added to the plot by setting 'arr.gene' to TRUE (default).
##            The spacing and sizes of gene names and arrows adjusts itself 
##            to some degree but for large regions with many genesi and small
##            regions with little space overlap is unavoidable.

loadSweeD <- function(path='./', common='', pop, ext='', 
                      nchr=8, skip=1) {
  sdata <- read.table(paste(path, '/SweeD_Report.chr1', common, pop, ext, 
                            sep=''), skip=skip, header=TRUE)
  chromosome <- rep('chr1', length(sdata[,1]))
  cummpos <- sdata$Position
  sdata <- data.frame(chromosome, cummpos, sdata)
  for ( i in c(2:nchr) ) {
    tmp <- read.table(paste(path, '/SweeD_Report.chr', i, common, pop, ext, 
                            sep=''), skip=skip, header=TRUE)  
    chromosome <- rep(paste('chr', i, sep=''), length(tmp[,1]))
    cummpos <- tmp$Position + max(sdata$cummpos) 
    sdata <- rbind(sdata, data.frame(chromosome, cummpos, tmp))
  }
  return(sdata)
}
  
splits <- function(x) {
  for ( i in 2:length(x[,1])) {
    if (x$Position[i] < x$Position[i-1]) {
      if (!exists('csplit')) {
        csplit <- x$cummpos[i] - ((x$cummpos[i] - x$cummpos[i-1])/2)
      } else {
        csplit <- rbind(csplit, x$cummpos[i] - 
                        ((x$cummpos[i] - x$cummpos[i-1])/2))
      }
    }
  }
  return(csplit)
}

###############################################################################
## Peak detection                                                           ###   
###############################################################################
sweedPeaks <- function(clr, w=1, win=5000, loess.span=NULL, look=5) {
  # returns the indices of peaks in the CLR profile
  # clr - dataframe with SweeD report
  # w - smoothing
  # win - window size
  # look - positions left and right of smooth peak to evaluate for having 
  #        higher Likelihood  
  peaks <- function(x, y, w=w, span=loess.span, ...) {
    #adapted from reply by 'whuber' to:
    #http://stats.stackexchange.com/questions/36309/
    #how-do-i-find-peaks-in-a-dataset
    require(zoo)
    n <- length(y)
    y.smooth <- loess(y ~ x, span=loess.span, ...)$fitted
    y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
    delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
    i.max <- which(delta <= 0) + w
    list(x=x[i.max], i=i.max, y.hat=y.smooth)
  }

  maxxpeak <- function(look, pk, datay) {
    # The smoothing from the 'peaks' function may result in the reported
    # peak not being exactly at the position with the highest Likelihood.
    # The function maxxpeaks looks for the actual peak in the vicinity.
    # parameters:
    # look - positions to check left and right of peak (pk)
    # pk - vector of peak indices from funtion 'peaks'
    # datay - data for which peaks were detected using function 'peaks'
    if (length(pk) == 0) { return(0)}
    for (i in pk) {
      n <- length(datay)
      pl <- ifelse(i>=look, i-look, 1) 
      pr <- ifelse(i<=n-look, i+look, n) 
      if (max(datay[pl:pr], na.rm=TRUE) == 0) { next }
      if (!exists('f_PEAK')) {
        f_PEAK <- seq(pl,pr)[which(datay[pl:pr] == max(datay[pl:pr], 
                                                       na.rm=TRUE))]
      } else {
        f_PEAK <- c(f_PEAK, seq(pl,pr)[which(datay[pl:pr] == 
                                             max(datay[pl:pr], 
                                                 na.rm=TRUE))])
      }
    }
    if (f_PEAK) { 
      return (unique(f_PEAK))
    } else {
      return(0)
    }
  }
 
  n <- length(clr[,1])
  indx <- 1:n

  for (i in 1:(1+trunc(n/win))) {
    from <- i * win + 1 - win
    to <- ifelse(i < trunc(n/win)+1, from + (win-1) , n)
    if (to <= from | max(clr$Likelihood[from:to]) == 0) { next }
    if (is.null(loess.span)) {
      loess.span <- (1/(to-from))*8 
    }
    aux  <- peaks(clr$cummpos[from:to], clr$Likelihood[from:to], w, loess.span) 
    aux.realpeak <- maxxpeak(look, aux$i, clr$Likelihood[from:to])
    aux.peak.indx <- indx[from:to][aux.realpeak]
    if (i < trunc(n/win)+1) {
      from <- from + 50; to <- to + 50
      aux1 <- peaks(clr$cummpos[from:to], clr$Likelihood[from:to], w, 
                    loess.span) 
      aux1.realpeak <- maxxpeak(look, aux1$i, clr$Likelihood[from:to])
      aux1.peak.indx <- indx[from:to][aux1.realpeak]
    }
    if (!exists('f_RES')) {
      f_RES <- c(aux.peak.indx, aux1.peak.indx)
    } else {
      f_RES <- c(f_RES, aux.peak.indx, aux1.peak.indx)
    }
  }
  return(unique(f_RES[which(!is.na(f_RES))]))
}

###############################################################################
## Get the significant peaks and the underlying genes                       ###
###############################################################################
sigPkGenes <- function(clr, pks, gff, cutoff=0.995, thr_meth='peak', 
                       maxPkDist=3000) {
  # returns a list with the significant CLR peaks and genes located close to it
  # clr - complete SweeD report dataframe 
  # pks - peaks in SweeD output found using funtion sweedPeaks
  # gff - gff table
  # cutoff - significance cutoff for CLR Likelihood (quantile)
  # thr_meth - threshold method: 'peak' for quantile of the detected peaks, 'all' for
  #            quantile of all estimated likelihoods
  # maxPkDist - maximum distance between SweeD peak and gene
  pk <- clr[pks,]
  if ( thr_meth == 'peak' ) { Q <- pk$Likelihood } else { Q <- clr$Likelihood } 
  thr_sub <- quantile(Q, cutoff)
  sigpeak <- pk[which(pk$Likelihood >= thr_sub),]
  for ( chr in 1:8 ) {
    tpeak <- sigpeak[which(sigpeak$chromosome == paste('chr', chr, sep='')),]
    if ( length(tpeak[,1]) > 0 ) {
      tgff <- gff[which(gff$seqname == paste('chr', chr, sep='') & 
                        gff$feature == 'gene'),]
      for ( i in 1:length(tpeak[,1]) ) {
        select <- tgff[which(abs((tgff$start-maxPkDist) - tpeak$Position[i]) <= 
                             maxPkDist | 
                             abs((tgff$end+maxPkDist) - tpeak$Position[i]) <= 
                             maxPkDist |
                             (tpeak$Position[i]>tgff$start & tpeak$Position[i]<
                              tgff$end )),]
        if (length(select[,1]) > 0) {
          select <- cbind(peak.position = tpeak$Position[i], select,
                          row.names=NULL)
          select <- cbind(peak.likelihood = tpeak$Likelihood[i], select, 
                          row.names=NULL)
          select <- cbind(peak.rank = 0, select, row.names=NULL)
        } else {
          select <- cbind(peak.position=tpeak$Position[i], peak.likelihood=
                          tpeak$Likelihood[i], peak.rank=0, seqname=
                          tpeak$chromosome[i], source=NA, feature=NA, start=NA,
                          end=NA, score=NA, strand=NA, frame=NA, gene=NA, 
                          annotation=NA, row.names=NULL) 
        }
        if (!exists('peakgene')) {
          peakgene <- select
        } else {
          peakgene <- rbind(peakgene, select)
        }
      }
    }
  }
  peakgene <- peakgene[with(peakgene, order(-as.numeric(peak.likelihood))),]
  for ( i in 1:length(peakgene[,1]) ) {
    if ( i == 1 ) { prank = 1 
    } else if (i>1 & peakgene$peak.position[i] != peakgene$peak.position[i-1]){
      prank = prank + 1
    }
    peakgene$peak.rank[i] = prank
  }
  sigpeak <- sigpeak[with(sigpeak, order(-Likelihood)),]
  return(list(sig.peak = sigpeak, peak.genes = peakgene))
}

# peak plots
peakPlot <- function(pk, pks, sigpks=pks, clr, gff, thr, thr_meth='peak', 
                     addspan, arr.gene=TRUE) {
  # Plots the CLR profiles zoomed in to peaks
  # pk - index of the peak in sigpks to plot.peak
  # pks - all detected peaks from function sweedPeaks for thr
  # sigpks - list of peaks of interest such as from sigPkGenes$sig.peak
  #          defaults to the all detected peaks
  # clr - data.frame with the original cLR scan results 
  # gff - gff table
  # thr - genome-wide threshold for likelihood determining the
  #       range shown in the plot
  # thr_meth - determine thr quantile from either the detected peaks ('peak')
  #            or from all CLR likelihoods ('all')
  # addspan - added range to plotted positions
  # arr.gene - whether to plot arrows for genes and print the gene id below
  
  threshold <- quantile(clr$Likelihood[pks], thr)

  findSpan <- function(pk, sigpks, clr, threshold) {
    pkchr <- clr[which(clr$chromosome == sigpks$chromosome[pk]),]
    begin <- which(pkchr$Position == sigpks$Position[pk]) 
    checkNext <- function(begin, clr, stride=10, n=1, direction='up', 
                          threshold) {
      count = n
      if (direction == 'up') { add = -stride }
      else if (direction == 'down') { add = stride }
      else { stop("direction specified should be either 'up' or 'down'") }
      if (begin > 2+stride & clr$Likelihood[begin + add] >= threshold) {
        begin = begin + add
        if (n > 4) { stride = 5 } else if (n>9) { stride = 1 }
        checkNext(begin, clr, stride=stride, n+n, direction=direction, 
                  threshold)
      } else {
        return(begin)
      }
    }
    chkup  <- pkchr$Position[checkNext(begin, pkchr, 10, 1, 'up', threshold)]
    chkdwn <- pkchr$Position[checkNext(begin, pkchr, 10, 1, 'down', threshold)]
    return(round(max(c(abs(pkchr$Position[begin]-chkup), 
                       abs(chkdwn-pkchr$Position[begin]))), 0) + 3500)
  }
  
  pspan <- findSpan(pk, sigpks, clr, threshold) + addspan

  clr.sub <- clr[which(clr$chromosome == sigpks$chromosome[pk] &
                       clr$Position >= (sigpks$Position[pk] - pspan) &
                       clr$Position <= (sigpks$Position[pk] + pspan)),] 
  gff.sub <- gff[which(gff$seqname == sigpks$chromosome[pk] &
                       (gff$start >= (sigpks$Position[pk] - pspan) |
                        gff$end >= (sigpks$Position[pk] - pspan)) &
                       (gff$end <= (sigpks$Position[pk] + pspan) |
                        gff$start <= (sigpks$Position[pk] + pspan))),]
  
  ylim <- c(-(max(clr.sub$Likelihood, na.rm=TRUE)/20), 
            1.05*max(clr.sub$Likelihood, na.rm=TRUE))
  xlim <- c(sigpks$Position[pk] - pspan, sigpks$Position[pk] + pspan)
  symb <- 16
  symb.cex <- 0.5
  col1 <- rgb(255, 50, 50, 150, maxColorValue=255) # red
  col2 <- rgb(50, 255, 0, 180, maxColorValue=255) # green 
  col3 <- rgb(50, 50, 255, 150, maxColorValue=255) # blue
  col4 <- rgb(150, 150, 150, 180, maxColorValue=255) # grey
  cols <- c(col1, col2, col3, col4)
  xlab <- paste('Physical position on', sigpks$chromosome[pk], '(bp)')
  
  plot(1, type='n', ylim=ylim, xlim=xlim, xlab=xlab, 
       ylab='Composite Likelihood Ratio', 
       las=ifelse(max(clr.sub$Likelihood) < 1000, 1, 0))   
  abline(h=c(0, threshold), col=c('black', cols[2]), lty=c(1,2))
  points(clr.sub$Position, clr.sub$Likelihood, pch=symb, cex=symb.cex)
  lines(clr.sub$Position, clr.sub$Likelihood, col=cols[1])
  lines(rep(sigpks$Position[pk],2), c(sigpks$Likelihood[pk], 0), lty=4, 
        col='magenta')
  
  if (arr.gene == TRUE) { 
    gff.sub.gene <- gff.sub[which(gff.sub$feature == 'gene'),]
    texty <- c(1.5*ylim[1]/3, 2.5*ylim[1]/3, 3.5*ylim[1]/3)
    ctext <- 1
    for ( i in 1:length(gff.sub.gene[,1]) ){
      if ( gff.sub.gene$strand[i] == 'plus' ) {
        c1 <- c(gff.sub.gene$start[i], ylim[1]/3)
        c2 <- c(gff.sub.gene$end[i], ylim[1]/3)
      } else {
        c2 <- c(gff.sub.gene$start[i], ylim[1]/3)
        c1 <- c(gff.sub.gene$end[i], ylim[1]/3)
      }
      arrows(c1[1], c1[2], c2[1], c2[2], length=0.125/
             (max(1, (clr.sub$Position[length(clr.sub$Position)] -
                      clr.sub$Position[1])/30000)), lwd=1.5)
      text(min(c(c1[1], c2[1])) + ((max(c(c1[1], c2[1])) - min(c(c1[1], 
                                                                 c2[1])))/2),
           texty[ctext], labels=gff.sub.gene$gene[i], cex=0.6, pos=1)
      if ( ctext < 3 ) { ctext = ctext +1 } else { ctext = 1 }
    }
  }
}
