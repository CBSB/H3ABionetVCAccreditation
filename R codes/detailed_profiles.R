rm(list=ls())
setwd("~/Dropbox/University life/PHD quest/Projects_Current/Accreditaion/")
{
  library(tidyverse)
  library(lubridate)
  require(gridExtra)
}
source('functions.R')

files = paste0(getwd(),'/data/timings') %>% list.files(pattern='timings.*', full.names=T)

####################### 3. Most expensive cpu process

Rscript detailed_profiles_alignment.R
####################### 4. Most expensive i/o process
{
  str(data)
  data %>%
    select(time, i.o.process) %>%
    separate(col = 'i.o.process', into = c('process','read_write'),sep = ' / ')  %>%
    separate(col = 'read_write', into = c('read','write'),sep = ':') %>%
    mutate(read = as.double(read))  %>%
    mutate(write = as.double(write)) %>%
    mutate(process = as.factor(process))  %>%
    mutate (read_write = (read + write)/(1024*1024)) %>%
    ggplot(aes(time, read_write)) + geom_line(aes(colour = process)) + 
    xlab('Time   (sec)') + ylab('Read & Write, Mega Byte')
}
####################### 5. Most expensive memory process
{
  str(data)
  data %>%
    select(time, memory.process) %>%
    separate(col = 'memory.process', into = c('process','memory'),sep = ' / ')  %>%
    separate(col = 'memory', into = c('memory'),sep = '%') %>%
    mutate(process = as.factor(process)) %>%
    mutate(memory = as.double(memory))  %>%
    mutate(memory = memory / (1024*1024*1024)) %>%
    ggplot(aes(time, memory)) + geom_line(aes(colour = process)) + 
    xlab('Time   (sec)') + ylab('Memory, GB')
  
}
# Complete pc plots
{ #Extra_stuff!
  ####################### 1. Total CPU usage (usr+sys)
  {
    d = data %>% 
      mutate(total_cpu = usr+sys)
    d%>%
      ggplot(aes(time, total_cpu), fill = cpu.process ) + geom_line() + xlab('Time   (sec)')
  }
  ####################### 2. sda utilization
  {
    str(data)
    data %>%
      select(time, util) %>%
      ggplot(aes(time, util)) + geom_line() + xlab('Time   (sec)')
  }
}