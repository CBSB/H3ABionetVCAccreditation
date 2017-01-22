### 1: deduplication stuff
dedup_profile = read.table('./data/profiles/profile.complete_pipeline_dedup.log', 
                          skip = 6, header = T,  sep = ',', stringsAsFactors = F) %>%  mutate(time = paste0('2017-', time)) %>% 
  mutate(time = force_tz(ydm_hms(time))) %>% cpu_process()

dedup_timing = crunching_times(files [ grep(pattern = 'dedup', x = files) ])

dedup_profile_cpu  = intersection(dedup_profile, dedup_timing)
plotting {
g.novosort = dedup_profile_cpu %>% filter(process == 'novosort') %>% 
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Novosort\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g.picard = dedup_profile_cpu %>% filter(process == 'java') %>% 
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Picard\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g.sambamba = dedup_profile_cpu %>% filter(process == 'sambamba') %>% 
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Sambamba\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g.samtools = dedup_profile_cpu %>% filter(process == 'samtools') %>% 
  filter(cores > 0 )%>%
  ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
  xlab('Time   (hours)') + ylab('Percentage CPU usage') + 
  scale_color_discrete(name='Samtools\n cores') + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(g.novosort, g.samtools, g.sambamba,g.picard, ncol = 2)
}