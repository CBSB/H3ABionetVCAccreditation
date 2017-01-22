### 1: alignment stuff

align_profile = read.table('./data/profiles/profile.complete_pipeline_align.log', 
                               skip = 6, header = T,  sep = ',', stringsAsFactors = F) %>%  mutate(time = paste0('2017-', time)) %>% 
  mutate(time = force_tz(ydm_hms(time))) %>% cpu_process()

align_timing = crunching_times(files [ grep(pattern = 'align', x = files) ])

alignment_profile_cpu  = intersection(align_profile, align_timing)
 {
  g.novoalign = alignment_profile_cpu %>% filter(process == 'novoalign') %>% 
    ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
    xlab('Time   (seconds)') + ylab('Percentage CPU usage') + 
    scale_color_discrete(name='Novoalign\n cores') + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  g.bwa = alignment_profile_cpu %>% filter(process == 'bwa') %>% 
    ggplot(aes(duration, percentage)) + geom_line(aes(colour = (cores))) +  
    xlab('Time   (seconds)') + ylab('Percentage CPU usage') + 
    scale_color_discrete(name='BWA MEM\n cores') + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  grid.arrange(g.bwa, g.novoalign, ncol = 2)
  
}