crunching = function(file) {
  {require(stringr)
    require(ggplot2)
    require(lubridate)}
  data = read.table(file, sep = '\t', stringsAsFactors = F, header = T)
  data$START_TIME = str_extract(data$START_TIME,'[0-9]+:[0-9]+:[0-9]+')
  data$END_TIME = str_extract(data$END_TIME, '[0-9]+:[0-9]+:[0-9]+')
  data = data %>% 
    separate(PROCESS, into = c('process', 'cores'), sep='_') %>%
    filter(!is.na(cores)) %>%
    mutate (process = as.factor(process)) %>%
    mutate (duration = as.numeric(hms(END_TIME)-hms(START_TIME) )) %>%
    select(process, cores, duration)
}

sum_non_threaded = function(data,clean, non_threaded){
  non_1 = filter(data, process==non_threaded) 
  non_1[1,3] = summarise(non_1, mean(duration))
  clean = rbind(clean, non_1[1,]) %>% na.omit()
}

prep_timing = function(stage, file, non_threaded=c(''), combined = F, 
                       combined_flag = 'novosort', ready_data_flag=F, ready_data, plot = T) {
  
  if (!ready_data_flag) {
    data = crunching(file)
    if (combined)
      data = filter(data, process != combined_flag)
    for (i in seq_along(non_threaded)){
      clean = filter (data, process != non_threaded[i]) 
      data  = sum_non_threaded(data, clean,non_threaded [i])
    }
  } else data = ready_data
  
  if (plot) {
    g = 
      ggplot(data=data, aes(x=cores, y=duration, fill=process)) + 
      geom_bar(colour="black", stat="identity", position='dodge', width = 0.7) +
      ylab('Time   (sec)') + xlab('Number of cores') + scale_fill_discrete(name=stage)
  } else data
}

##########
crunch_annotation = function (file) {
  data = read.table(file, header = FALSE, stringsAsFactors = F) %>%
    filter(V1 == 'chr1')
  header = read.table(file, nrows = 1, comment.char = '', 
                      stringsAsFactors = F) %>% select(-1) %>%  
    str_replace('CBSB_Khartoum:','') %>% str_extract('([A-Za-z])+')
  names(data) = header
  data
}

### Extract columns that are relevant to the profiling of cpu usage of process
cpu_process = function(data){
  cpu = data %>%
    select(time, cpu.process) %>%
    separate(col = 'cpu.process', into = c('process','percentage'),sep = ' / ') %>%
    separate(col = 'percentage', into = c('percentage'), sep = '%') %>%
    mutate(process = as.factor(process)) %>%
    mutate(percentage = as.double(percentage))
}

### Read the timings of processes
crunching_times = function(file) {
  #file=files [ grep(pattern = 'vcall', x = files) ]
  {require(stringr)
    require(ggplot2)
    require(lubridate)}
  data = read.table(file, sep = '\t', stringsAsFactors = F, header = T)
  data$START_TIME = str_extract(data$START_TIME,'[0-9]+:[0-9]+:[0-9]+')
  data$END_TIME = str_extract(data$END_TIME, '[0-9]+:[0-9]+:[0-9]+')
  data = data %>% 
    separate(PROCESS, into = c('process', 'cores'), sep='_') %>%
    mutate (process = as.factor(process)) %>%
    mutate (END_TIME = paste(ymd('2017-01-21'),END_TIME))%>%
    mutate(END_TIME = force_tz(ymd_hms(END_TIME))) %>%
    mutate (START_TIME = paste(ymd('2017-01-15'),START_TIME))%>%
    mutate(START_TIME = force_tz(ymd_hms(START_TIME))) %>%
    mutate(interval = interval(START_TIME,END_TIME))
  #filter(!is.na(cores)) %>%
  data
}

### between features and timing
intersection = function (feature, timing){
  #timing = dedup_timing
  #feature = dedup_profile
  logical_index = feature$time %in% timing$START_TIME
  index = which(logical_index)
  breaks = c(1,index,nrow(feature))
  t_factor = cut(feature$time, breaks = feature[breaks,'time'], include.lowest = TRUE )
  
  t_start_instant = force_tz(ymd_hms(as.character(t_factor)))
  
  feature = cbind(feature, t_factor,t_start_instant) %>%
    mutate(duration = time-t_start_instant) %>%
    mutate(cores = 0) %>%
    mutate(process = as.character(process))
  str(feature)
  head(cbind(feature[,c('time','t_start_instant','duration')]))
  
  for (i in seq_along(1:nrow(feature)))
    for (j in seq_along(1:nrow(timing)))
      #i=508+1;j=1; (feature[i,'time'] %within% timing[j,'interval'])
      if (feature[i,'time'] %within% timing[j,'interval'])
        feature[i,'cores'] = timing[j,"cores"]
  feature
}
  