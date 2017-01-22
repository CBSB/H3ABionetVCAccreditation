rm(list=ls())
library(tidyverse)
library(lubridate)
require(gridExtra)
source('functions.R')
setwd("~/Dropbox/University life/PHD quest/Projects_Current/Accreditaion/")

rawdata = read.table('./data/profiles/profile.complete_pipeline_vcall.log', skip = 6, 
                     header = T,sep = ',', stringsAsFactors = F)
profile_data = rawdata

profile_data[,1] = paste0('2017-', profile_data[,1])
profile_data[,1] = force_tz(ydm_hms(profile_data[,1]))


files = paste0(getwd(),'/data/timings') %>% list.files(pattern='timings.*', full.names=T)
vcall_file = files [ grep(pattern = 'vcall', x = files) ]
vcall_timing = crunching_times(vcall_file)

####################### 3. Most expensive cpu process
{
 stages = as.character(vcall_timing$process)
 interval = vcall_timing$interval
  ####
  for (i in seq_along(stages)){
    p = profile_data %>%
      select(time, cpu.process) %>%
      separate(col = 'cpu.process', into = c('process','percentage'),sep = ' / ') %>%
      separate(col = 'percentage', into = c('percentage'), sep = '%')
    p = p[p$time %within%interval[i],] 
    p[,1] = p[,1]-p[1,1] 
    png(paste0('results/cpu_profile/',stages[i],i))
    p%>%
      mutate(percentage = as.double(percentage))%>%
      ggplot(aes(time, percentage)) + geom_line() + 
      xlab(paste('Time   (sec) running', stages[i]))
    dev.off()
  }

####################### 4. Most expensive i/o process
{
  str(profile_data)
  for (i in seq_along(stages)){ 
    #i = 0; i=i+1
    p = profile_data %>%
      select(time, i.o.process) %>%
      separate(col = 'i.o.process', into = c('process','read_write'),sep = ' / ')  %>%
      separate(col = 'read_write', into = c('read','write'),sep = ':') %>%
      mutate(read = as.double(read))  %>%
      mutate(write = as.double(write)) %>%
      mutate(process = as.factor(process))  %>%
      mutate (read_write = (read + write)/(1024*1024)) 
    p = p[p$time %within%interval[i],] 
    p[,1] = p[,1]-p[1,1] 
    png(paste0('results/io/',stages[i],i))
    p%>%
      ggplot(aes(time, read_write)) + geom_line() + 
      xlab(paste('Time   (sec) running', stages[i])) + ylab('Read & Write, Mega Byte')
    dev.off()
  }
}
####################### 5. Most expensive memory process
{
  str(profile_data)
  for (i in seq_along(stages)){ 
    #i = 0; i=i+1
    p = profile_data %>%
      select(time, memory.process) %>%
      separate(col = 'memory.process', into = c('process','memory'),sep = ' / ')  %>%
      separate(col = 'memory', into = c('memory'),sep = '%') %>%
      mutate(process = as.factor(process)) %>%
      mutate(memory = as.double(memory))  %>%
      mutate(memory = memory / (1024*1024*1024)) 
    p = p[p$time %within%interval[i],] 
    p[,1] = p[,1]-p[1,1] 
    png(paste0('results/memory/',stages[i],i))
    p%>%
      ggplot(aes(time, memory)) + geom_line() + 
      xlab(paste('Time   (sec) running', stages[i])) + ylab('Memory, GB')
    dev.off()
}

####################### 1. Total CPU usage (usr+sys)
{
  profile_data %>% 
  select(time, usr, sys) %>%
  mutate(total_cpu = usr+sys) %>%
  ggplot(aes(time, total_cpu)) + geom_line() + xlab('Time   (sec)')
}
####################### 2. sda utilization
{
  str(profile_data)
  profile_data %>%
    select(time, util) %>%
    ggplot(aes(time, util)) + geom_line() + xlab('Time   (sec)')
}