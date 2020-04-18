library(tidyverse)
CF = read_delim("G:/gBGC/carl/workflow_fresh/output/CF_collected.tab", 
                           "\t", escape_double = FALSE, trim_ws = TRUE,
                 col_names = c('genospecies', 'unitig', 'file', 'parameter', 'post_mean', 'post_var', 'a_post', 'b_post')) %>% 
    mutate(gene = str_sub(file, 4, str_length(file)-7))
CF %>% View



GC3 = read_delim("G:/gBGC/carl/workflow_fresh/output/GC3.tab", 
                  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    rename(genospecies = gs)




