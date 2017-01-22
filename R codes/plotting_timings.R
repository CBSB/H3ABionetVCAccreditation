dev.off()
rm(list=ls())
require(gridExtra)
source('functions.R')
setwd("~/Dropbox/University life/PHD quest/Projects_Current/Accreditaion/")
path = paste0(getwd(),'/data/timings')

files=list.files(pattern='timings.*',path=path, full.names = T)

align_file = grep(pattern = 'align', x = files)
  g.align = prep_timing(stage = 'Alignment \n  tool',file = files[align_file],plot = T)

sort_file = grep(pattern = 'sort', x = files)
  g.sort =prep_timing(stage ='Sorting\n tool',file = files[sort_file], plot =T)

index_file = grep(pattern = 'index', x = files)
 g.index = prep_timing(stage = 'Indexing \n tool', file = files[index_file], plot =T)

dedup_file = grep(pattern = 'dedup', x = files)
g.dedup = prep_timing(stage = 'Deduplication \n  tool', file = files[dedup_file], plot =T)
g.dedup
sorting = prep_timing(stage = 'sorting', file = files[sort_file], plot = F, ready_data_flag = F)
indexing = prep_timing(stage = 'indexing', file = files[index_file], plot = F, ready_data_flag = F)
deduplication = prep_timing(stage = 'deduping', file = files[dedup_file], plot = F, ready_data_flag = F)

data = rbind(sorting, indexing, deduplication)
data = data %>% group_by(process, cores) %>%
  summarize(duration = sum(duration))
g.all = prep_timing(stage = 'Sorting, indexing \n& deduplication',ready = T, ready_data = data)

g.align

grid.arrange(arrangeGrob(g.sort,g.dedup+ ylim(0,400)),g.all, ncol=2)
