library(tidyverse)
rm(list = ls())

mcorr_colnames = c("gene_str", "group", "d_sample", "theta_pool", "phi_pool", "ratio", "fbar", "c", "d_pool", "d_clonal", "theta_s", "phi_s")
setwd("~/genomedk/gBGC/carl/workflow_assemblyalignmentgenerator/")
recomb_data = read_csv("output/Spyo3/Spyo3_final.csv", 
                         col_names = mcorr_colnames) %>% 
    mutate(gene = str_sub(gene_str, 1, str_length(gene_str)-19))



gc_data = read_delim("output/Spyo3/Spyo3_gc.tab", delim = "\t", col_names = c("gene", "GC3"))



data = inner_join(recomb_data, gc_data)

data %>% filter(phi_pool >5) %>% 
    ggplot(aes(GC3, log(phi_pool))) +
    geom_smooth() +
    geom_point()
ggsave("1.1_Spyo.png")

# bin

data %>% drop_na() %>% 
    #filter(log(phi_pool) > 5) %>% 
    mutate(GC3_bin = cut_number(GC3, 20)) %>% 
    group_by(GC3_bin) %>% 
    summarise(mean_GC3 = mean(GC3), median_phi_pool = median(phi_pool)) %>% 
    ggplot(aes(mean_GC3, median_phi_pool)) + geom_point() + geom_smooth()


ggsave("1.2_Spyo_20_bins.png")


