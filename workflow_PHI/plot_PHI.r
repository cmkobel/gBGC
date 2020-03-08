library(tidyverse)

path = "/home/arabthal/genomedk/gBGC/carl/workflow_PHI/output/export"
recomb_files = list.files(path=path, pattern="*phi_results.tab", full.names=TRUE, recursive=FALSE)
recomb_data = tibble()
for (file in recomb_files) {
    import = read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE, na = c("--"))
    print(paste(file, dim(import)[1]))
    recomb_data = bind_rows(recomb_data, import)
    rm(import)
}

recomb_data = recomb_data %>% mutate(method = paste(method, detail),
       unitig = str_sub(genome, 8, 8),
       genospecies = str_sub(genome, 24, 24)) %>% 
    select(-detail)

alpha = 0.05

recomb_data_tresh = recomb_data %>% 
    group_by(method, genome) %>% 
    mutate(recombine = if_else(pvalue < alpha/length(pvalue), T, F), n_genes = length(pvalue))

# sanity check of pvalues:
recomb_data_tresh %>%
    ggplot(aes(pvalue, fill = recombine)) + 
    geom_histogram() + 
    scale_y_log10() +
    facet_grid(method~genome)




gc_files = list.files(path=path, pattern="*gc.tab", full.names=TRUE, recursive=FALSE)
gc_data = tibble()
for (file in gc_files) {
    import = read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)
    print(paste(file, dim(import)[1]))
    gc_data = bind_rows(gc_data, import)
    rm(import)
}

gc_data_summarized = gc_data %>%
    rename(genome = cli_comment_2) %>% select(-cli_comment_1) %>% 
    group_by(genome, gene) %>% 
    summarize(GC3 = mean(gc_content))




data = inner_join(gc_data_summarized, recomb_data_tresh)
write_csv(data, "gc_and_recomb_joined.csv")


# bin the GC content

data_binned20 = data %>% 
    group_by(genome, method) %>% 
    mutate(GC3_bin = cut_number(GC3, 20)) %>% # View
    group_by(GC3_bin, add = T) %>%
    mutate(mean_GC3 = mean(GC3)) %>% # View
    group_by(genospecies, unitig, genome, method, mean_GC3) %>% # use mean instead of weird bin-range
    count(recombine) %>%
    spread(recombine, n) %>% 
    rename(n_recombining = `TRUE`, n_not_recombining = `FALSE`) %>% 
    mutate(ratio = n_recombining / (n_recombining + n_not_recombining))
write_csv(data_binned20, "binned20.csv")


# Main plot
data_binned20 %>% 
    #filter(method == "PHI (Normal):") %>%
    filter(method == "PHI (Permutation):") %>%
    
    filter(unitig == "0") %>%
    
    ggplot(aes(mean_GC3, ratio)) +
    #ggplot(aes(mean_GC3, ratio)) +
    
    geom_point() +
    facet_wrap(~genome, scales = "free") + 
    geom_smooth(method = "lm")



setwd("~/genomedk/gBGC/carl/workflow_PHI/")
height = 5; width = 8
ggsave("main_n.png", height = height, width = width)
    #summarize(mean_GC3 = mean(GC3))





    