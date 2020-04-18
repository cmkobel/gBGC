library(tidyverse)
CF_raw = read_delim("G:/gBGC/carl/workflow_fresh/output/CF_collected.tab", 
                           "\t", escape_double = FALSE, trim_ws = TRUE,
                 col_names = c('genospecies', 'unitig', 'file', 'parameter', 'post_mean', 'post_var', 'a_post', 'b_post')) %>% 
    mutate(gene = str_sub(file, 4, str_length(file)-7)) %>% select(-file)


#CF = CF_raw %>% pivot_wider(names_from = parameter, values_from = c(post_mean, post_var, a_post, b_post))
CF = CF_raw %>% filter(parameter == "R/theta") %>% select(-parameter)


GC3_raw = read_delim("G:/gBGC/carl/workflow_fresh/output/GC3.tab", 
                  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    rename(genospecies = gs) %>%
    mutate(sample = str_sub(header, 6),
           gene = str_sub(file, 8, str_length(file)-13)) %>% 
    select(-header, -file)

GC3 = GC3 %>% group_by(genospecies, unitig, gene) %>% 
    summarize(mean_GC3 = mean(GC3), n_samples_GC3 = length(GC3))
    



gff_data <- read_csv("C:/Users/carl/Desktop/xmfa2mfa/3206-3.gff") %>% select(-X1) %>% mutate(mid = start - ((end-start)/2))
gff2_data <- read_delim("C:/Users/carl/Desktop/xmfa2mfa/Gene_function_pop_gene_data.csv", delim = ';')
gff_data = left_join(gff_data, gff2_data %>% select(`Gene group`, `Putative function`) %>% rename(gene_group = `Gene group`)) %>% rename(gene = gene_group)
rm(gff2_data)


data = inner_join(GC3, CF) %>% inner_join(gff_data)


data %>% ggplot(aes(



