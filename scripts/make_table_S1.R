


library(tidyverse)
library(readxl)

sample_data <- 
  read_delim("data/processed/pig_filtered_sample_data.tab", delim = "\t")

dat <- read_excel("doc/Manuscript/supplementary data/S1.xlsx", sheet = 2)

meta <- read_excel("doc/Manuscript/supplementary data/S1.xlsx", sheet = 3)

dat %>% 
  mutate(`Tissue ID` = gsub("\\d-", "", Sample)) %>% 
  left_join(meta %>% 
              select(`Tissue ID`, Tissue)) %>% 
  select(Sample, Tissue, everything(), -`Tissue ID`) %>% 
  filter(!complete.cases(.))

samples_in_data <-
  sample_data %>% 
  select(sample_ID, individual) %>% 
  distinct() %>% 
  mutate(sample = paste0(match(individual, letters), "-", gsub("_.$", "", sample_ID)))
  
samples_in_data

dat$Sample[!dat$Sample %in% samples_in_data$sample]
         
savedata <- 
  dat %>% 
  filter(Sample %in% samples_in_data$sample) %>% 
  mutate(`Tissue ID` = gsub("\\d-", "", Sample)) %>% 
  left_join(meta %>% 
              select(`Tissue ID`, Tissue)) %>% 
  select(Sample, Tissue, everything(), -`Tissue ID`) %>% 
  arrange(Sample) 

savedata %>% 
  write_csv("data/processed/sup table S1.csv")


savedata %>% 
  slice(66)

savedata %>% 
  filter(grepl("cer", Sample) & 
           grepl("1|2", Sample)) %>% 
  write_csv("data/processed/cervix_QC.csv")
