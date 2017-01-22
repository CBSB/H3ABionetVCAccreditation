### 1: sorting stuff

sort_profile = read.table('./data/profiles/profile.complete_pipeline_sort.log', 
                           skip = 6, header = T,  sep = ',', stringsAsFactors = F) %>%  mutate(time = paste0('2017-', time)) %>% 
  mutate(time = force_tz(ydm_hms(time))) %>% cpu_process()

sort_timing = crunching_times(files [ grep(pattern = 'sort', x = files) ])

sorting_profile_cpu  = intersection(sort_profile, sort_timing)
plotting {
g.novosort = sorting_profile_cpu %>% filter(process == 'novosort') %>% 
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Novosort\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g.picard = sorting_profile_cpu %>% filter(process == 'java') %>% 
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Picard\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g.sambamba = sorting_profile_cpu %>% filter(process == 'sambamba') %>% 
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Sambamba\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g.samtools = sorting_profile_cpu %>% filter(process == 'samtools') %>% 
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Samtools\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


grid.arrange(g.samtools, g.novosort, g.sambamba, ncol = 3)
}
