rm(list=ls())

{if (!require(ggplot2)){
  install.packages(ggplot2)
  library(ggplot2)
}
  library(stringr)
  library(gridExtra)
}
source('functions.R')

indelfile = './data/filteration/CBSB_Khartoum.gatk.tabulted_annotations.indels.txt'
gatkindels = crunch_annotation(indelfile)

################
g.qd = ggplot(data = gatkindels, aes(x=QD)) + 
  labs(y="Density") + scale_fill_discrete('Variant caller') + 
  geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
  theme_bw() + geom_vline(xintercept = 2, linetype=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x=3, y=0.01, label='Minimum recommended QD', hjust=0, vjust=-0.9,
           angle=90, size=4, fontface=3, family='serif')

g.fs = ggplot(data = gatkindels, aes(x=FS)) + 
  labs(y="Density") + scale_fill_discrete('Variant caller') + 
  geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
  theme_bw() + geom_vline(xintercept = 200, linetype=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x=199, y=0.5, label='Maximum recommended QD', hjust=0, vjust=1.5,
           angle=90, size=4, fontface=3, family='serif') 

g.sor = ggplot(data = gatkindels, aes(x=SOR)) + 
  labs(y="Density") + scale_fill_discrete('Variant caller') + 
  geom_density(aes(fill=caller), fill = 'blue', alpha=.4, color='blue') + 
  theme_bw() + geom_vline(xintercept = 10, linetype=2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  annotate("text", x=9.9, y=0.4, label='Maximum recommended QD', hjust=0, vjust=1.5,
           angle=90, size=4, fontface=3, family='serif') 

grid.arrange(g.qd, g.fs, g.sor, ncol = 3)
