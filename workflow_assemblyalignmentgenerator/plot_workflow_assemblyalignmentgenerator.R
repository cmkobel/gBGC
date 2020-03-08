library(tidyverse)

Spyo1_gc <- read_delim("~/genomedk/gBGC/carl/workflow_assemblyalignmentgenerator/output/Spyo1/Spyo1_gc.tab", 
                       "\t", escape_double = FALSE, col_names = c("gene", "GC3"), 
                       trim_ws = TRUE)


Spyo1_recomb <- read_csv("~/genomedk/ClinicalMicrobio/faststorage/10carl/urecomb/workflow_assemblyalignmentgenerator/output/Spyo1/all_genes_recomb_data.csv", col_names = c("source", "group", "d_sample", "theta_pool", "phi_pool", "ratio", "fbar", "c", "d_pool", "d_clonal", "theta_s", "phi_s")) %>% 
    mutate(gene = str_sub(source, 1, str_length(source)-22)) %>% select(-group, -source)



data = inner_join(Spyo1_gc, Spyo1_recomb)

tail_size_phi_pool = 0.1
data_inliers = data %>%
    filter(phi_pool > quantile(phi_pool, tail_size_phi_pool) & phi_pool < quantile(phi_pool, 1-tail_size_phi_pool)) 


data %>% filter(phi_pool > 50) %>% ggplot(aes(GC3, phi_pool)) +
    geom_point() + 
    #geom_hline(yintercept = 50) +
    scale_y_log10() + 
    geom_smooth(method = "lm") + 
    labs(title = "Spyogenes", subtitle = "64 randomly chosen genomes from genbank", caption = "phi_pool-values below 50 are excluded\nRecombination is inferred with mcorr on each gene at a time in the core genome.")

ggsave("Spyo mcorr.png")
