{if (!require(ggplot2)){
  install.packages(ggplot2)
  library(ggplot2)
}
library(stringr)
library(gridExtra)
}

gatk.allfile = './data/filteration/gatk.all.txt'
gatkall = crunch_annotation(gatk.allfile) %>%
  mutate(caller = 'gatk')
freebayes.allfile = './data/filteration/freebayes.all.txt' 
freebayesall = crunch_annotation(freebayes.allfile) %>%
  mutate(caller = 'freebayes')

  ################################ plot everything:
  {#together
    data = all
    require(reshape2)
    long = melt(data[c(-2,-11,-12)], id.vars= "CHROM",variable.name = 'series')
    ggplot(long, aes(value)) + geom_density() + facet_grid(series ~ .)  
    ggplot(long, aes (value)) +   geom_density() +   facet_wrap(~series)
  }
  
  #################################
annotation = c( 'QD', 'MQ', 'FS', 'SOR')
  g.qual = ggplot(data = gatkall, aes(x=QUAL)) + 
    labs(y="Density") + scale_fill_discrete('Variant caller') + 
    geom_density(aes(fill=caller), alpha=.4) + 
    geom_density(data = freebayesall, aes(fill=caller), alpha=.4) + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  g.dp = ggplot(data = gatkall, aes(x=DP)) +
    labs(y="Density") + scale_fill_discrete('Variant caller') + 
    geom_density(aes(fill=caller), alpha=.4) + 
    geom_density(data = freebayesall, aes(fill=caller), alpha=.4) +
    geom_vline(xintercept = 10, linetype=2) + 
    geom_vline(xintercept = 85, linetype=2) + 
    annotate("text", x=10, y=0.005, label='Minumum recommended DP', 
             hjust=0, vjust=-0.9, angle=90, size=4, fontface=3, family='serif') +
    annotate("text", x=85, y=0.005, label='Maximum recommended DP', 
             hjust=0, vjust=1.5, angle=90, size=4, fontface=3, family='serif') + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  grid.arrange(g.dp, g.qual,ncol = 2)
  
  dev.off()

  g.qd = ggplot(data = gatkall, aes(x=QD)) + 
    labs(y="Density") + scale_fill_discrete('Variant caller') + 
    geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
     theme_bw() + geom_vline(xintercept = 2, linetype=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    annotate("text", x=2, y=0.05, label='Minimum recommended QD', hjust=0, vjust=-0.9,
             angle=90, size=4, fontface=3, family='serif')
  
  g.mq = ggplot(data = gatkall, aes(x=MQ)) + 
    labs(y="Density") + scale_fill_discrete('Variant caller') + 
    geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
    theme_bw() + geom_vline(xintercept = 40, linetype=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    annotate("text", x=40, y=0.05, label='Minimum recommended QD', hjust=0, vjust=-0.9,
             angle=90, size=4, fontface=3, family='serif') 
    
  
  g.fs = ggplot(data = gatkall, aes(x=FS)) + 
    labs(y="Density") + scale_fill_discrete('Variant caller') + 
    geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
    theme_bw() + geom_vline(xintercept = 60, linetype=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    annotate("text", x=60, y=0.05, label='Maximum recommended QD', hjust=0, vjust=1.5,
             angle=90, size=4, fontface=3, family='serif') 
  
  g.sor = ggplot(data = gatkall, aes(x=SOR)) + 
    labs(y="Density") + scale_fill_discrete('Variant caller') + 
    geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
    theme_bw() + geom_vline(xintercept = 3, linetype=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    annotate("text", x=3, y=0.05, label='Maximum recommended QD', hjust=0, vjust=1.5,
             angle=90, size=4, fontface=3, family='serif') 
  
  grid.arrange(g.qd, g.mq, g.fs, g.sor, ncol = 2)