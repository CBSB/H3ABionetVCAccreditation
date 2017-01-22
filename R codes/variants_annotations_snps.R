rm(list=ls())

{if (!require(ggplot2)){
  install.packages(ggplot2)
  library(ggplot2)
}
  library(stringr)
  library(gridExtra)
}
source('functions.R')
snpfile = './data/filteration/CBSB_Khartoum.gatk.tabulted_annotations.snps.txt'
gatksnps = crunch_annotation(snpfile)

#################################
annotation = c( 'QD', 'MQ', 'FS', 'SOR')
################
g.qd = ggplot(data = gatksnps, aes(x=QD)) + 
  labs(y="Density") + scale_fill_discrete('Variant caller') + 
  geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
  theme_bw() + geom_vline(xintercept = 2, linetype=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x=2, y=0.05, label='Minimum recommended QD', hjust=0, vjust=-0.9,
           angle=90, size=4, fontface=3, family='serif')

g.mq = ggplot(data = gatksnps, aes(x=MQ)) + 
  labs(y="Density") + scale_fill_discrete('Variant caller') + 
  geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
  theme_bw() + geom_vline(xintercept = 40, linetype=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x=40, y=0.05, label='Minimum recommended QD', hjust=0, vjust=-0.9,
           angle=90, size=4, fontface=3, family='serif') 

g.fs = ggplot(data = gatksnps, aes(x=FS)) + 
  labs(y="Density") + scale_fill_discrete('Variant caller') + 
  geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
  theme_bw() + geom_vline(xintercept = 60, linetype=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x=60, y=0.05, label='Maximum recommended QD', hjust=0, vjust=1.5,
           angle=90, size=4, fontface=3, family='serif') 

g.sor = ggplot(data = gatksnps, aes(x=SOR)) + 
  labs(y="Density") + scale_fill_discrete('Variant caller') + 
  geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
  theme_bw() + geom_vline(xintercept = 3, linetype=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x=3, y=0.05, label='Maximum recommended QD', hjust=0, vjust=1.5,
           angle=90, size=4, fontface=3, family='serif') 

grid.arrange(g.qd, g.mq, g.fs, g.sor, ncol = 2)
