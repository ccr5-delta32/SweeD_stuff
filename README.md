# SweeD_stuff
a bunch of functions useful to analyze and plot SweeD results

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
##            regions with little space overlap is unavoidable.
