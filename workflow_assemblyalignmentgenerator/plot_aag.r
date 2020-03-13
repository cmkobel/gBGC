library(tidyverse)
rm(list = ls())

mcorr_colnames = c("gene_str", "group", "d_sample", "theta_pool", "phi_pool", "ratio", "fbar", "c", "d_pool", "d_clonal", "theta_s", "phi_s")
setwd("~/genomedk/gBGC/carl/workflow_assemblyalignmentgenerator/")
recomb_data = read_csv("output/Spyo1/Spyo1_recomb.csv", 
                         col_names = mcorr_colnames) %>% 
    mutate(gene = str_sub(gene_str, 1, str_length(gene_str)-22))



gc_data = read_delim("output/Spyo1/Spyo1_gc.tab", delim = "\t", col_names = c("gene", "GC3"))



data = inner_join(recomb_data, gc_data)

data %>% ggplot(aes(GC3, phi_pool)) + geom_point()



