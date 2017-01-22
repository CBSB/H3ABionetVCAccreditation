setwd('~/accreditation/')

## Done in accordance with bedtools protocols - 
# https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md
library(ggplot2)
library(RColorBrewer)

#BP2: Assessing coverage in DNA sequencing experiments (WGS)
{
  ## Plotting Genome-Wide coverage distribution
  {
    cov = read.table('sample.coverage.hist.txt')
    head(cov)
    # extract the genome-wide (i.e., no the per-chromosome) histogram entries
    gcov = cov[cov[,1] == 'genome',]
    head(gcov)
    # plot a density function for the genome-wide coverage
    plot(gcov[1:51,2], gcov[1:51,5], type='h', col='darkgreen', lwd=3, 
         xlab="Depth", ylab="Fraction of genome at depth",
         main = 'Density of genome wide coverage')
    #axis(1,at=c(5,10,15,20,25,30,35,40,45,50))
    plot(gcov[,2], gcov[,5], type='h', col='darkgreen',
         lwd=3, xlab="Depth", ylab="Fraction of genome at depth")
  }
  
  ## Plotting a cumulative density function to assess what fraction of the genome
  ## is covered by more that 1, 10, 20, ...
  {
    gcov_cumul = 1 - cumsum(gcov[,5])
    # Create a plot of the CDF
    plot(gcov[2:51,2], gcov_cumul[1:50], col='darkgreen', type='l', lwd=3, 
         xlab="Depth", ylab="Fraction of genome >= depth", ylim=c(0,1.0),
         main = 'Culmutative Depth distribution')
    #axis(1,at=c(5,10,15,20,25,30,35,40,45,50))
    #axis(2,at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
  }
  
}


#AP2b: Assessing coverage in exome capture experiments 
{
  # bedtools coverage \
  # -hist \
  # -b <sample.exome.bam> \
  # -a <targets.numeric.chroms.bed> \
  # > <sample.exome.coverage.hist.txt>
  
  cov = read.table('sample.exome.coverage.hist.txt',fill = TRUE);
  head(cov)
  names(cov) [c(1,7:10)] = c('Region', 'Depth', 'Count bases at depth', 'Interval size', 'Fraction_bases')
  gcov = cov[cov[,1] == 'all',]
  names(gcov) = c('Region', 'Depth', 'Count bases at depth', 'Interval size', 'Fraction_bases')
  
  ## Plotting coverage distribution:
  {
    # plot a density function for the genome-wide coverage
    g = ggplot(gcov, aes(x = Depth, y =Fraction_bases ))
    g + geom_bar(position=position_stack(), stat="identity", 
                 aes(color='red'), width = .1) + guides(color=FALSE) +
      ggtitle('Coverage Distribution in Captured exonic regions') +
    theme(plot.title = element_text(lineheight=.8, face="bold"))
  }
  
  {
    # Create a cumulative distribution from the "raw" hist 
    # (truncate at depth >=1000)
    Cumultative_Coverage = 1 - cumsum(gcov[,5])
    gcov = cbind(gcov, gcov_cumul)
    # Create a plot of the CDF
    g = ggplot(gcov, aes(x = Depth, y =Cumultative_Coverage ))
    g + geom_line() + scale_color_brewer(palette="Set1") + 
       ggtitle('Cumultative coverage distribution \n in Exomic regions') #+
      # geom_vline(xintercept = 20, linetype=2) +
      # geom_vline(xintercept = 50, linetype=2) +
      # geom_vline(xintercept = 80, linetype=2) +
      # geom_hline(yintercept = .5, linetype=2) +
      # geom_hline(yintercept = .9, linetype=2)
  }
  
  