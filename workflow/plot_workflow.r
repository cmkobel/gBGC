# todo: extend this to include phi-results

rm(list = ls())
library(tidyverse)
setwd("~/urecomb/lokal/exports/export1")
#setwd("~/genomedk/gBGC/carl/workflow/output/")


# Chapter 1: import the data, and make some plots.

# import fitted (recombination) parameters
mcorr_files = list.files(path=".", pattern="*_fitpars.csv", full.names=TRUE, recursive=T)
mcorr_colnames = c("gene", "group", "d_sample", "theta_pool", "phi_pool", "ratio", "fbar", "c", "d_pool", "d_clonal", "theta_s", "phi_s", "bin_size", "genome")

mcorr_data = tibble()
i  = 1
for (file in mcorr_files) {
    
    import = read_csv(file, col_names = mcorr_colnames)
    print(paste(i, file, dim(import)[1]))
    mcorr_data = bind_rows(mcorr_data, import)
    rm(import)
    i = i + 1
}


# import gc data
gc_files = list.files(path=".", pattern="*bp_gc.tab", full.names=TRUE, recursive=T)
gc_data = tibble()
i = 1
for (file in gc_files) {
    import <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)
    print(paste(i, file, dim(import)[1]))
    gc_data = bind_rows(gc_data, import)
    #rm(import)
    i = i + 1
}

gc_data_summarised = gc_data %>% group_by(gene, bin_size = cli_comment_1, genome = cli_comment_2) %>% 
    summarise(GC3 = mean(gc_content))

data = inner_join(mcorr_data, gc_data_summarised) %>% mutate(gs = str_sub(genome, str_length(genome)))
#saveRDS(data, "data_workflow.rds")
#data = readRDS("data_workflow.rds")

height = 8
width = 10

# main plot
data %>%  ggplot(aes(GC3, phi_s, color = gs)) +
    geom_point(alpha = 0.5) + 
    scale_y_log10() +
    facet_grid(bin_size~gs)
#ggsave("main_(phi_pool_filtered).png", height = height, width = width)

# various histograms
data %>% ggplot(aes(phi_s, fill = gs)) +
    geom_histogram() + 
    facet_grid(bin_size~gs)




# Chapter 2
#    A) remove outliers (values 1-99th percentiles)



tail_size_phi_pool = 0.1
tail_size_ratio = 0.01
tail_size_c = 0.00
tail_size_phi_s = 0.02

data_inliers = data %>%
    group_by(bin_size, gs) %>% 
    
    filter(phi_pool >= quantile(phi_pool, tail_size_phi_pool) & phi_pool <= quantile(phi_pool, 1-tail_size_phi_pool)) %>% 
    filter(ratio >= quantile(ratio, tail_size_ratio) & ratio <= quantile(ratio, 1-tail_size_ratio)) %>% 
    filter(c >= quantile(c, tail_size_c) & c <= quantile(c, 1-tail_size_c)) %>% 
    
    ungroup()


# sanity check of tail-removal
# various histograms
data_inliers %>% ggplot(aes(phi_pool, fill = gs)) +
    geom_histogram() + 
    facet_grid(bin_size~gs)


# main plot
data_inliers %>%
    #filter(bin_size == 30000) %>% 
    ggplot(aes(GC3, log10(phi_pool), color = gs)) +
    geom_point(alpha = 0.5) + 
    geom_smooth() +
    facet_grid(bin_size~gs)
    
#ggsave("main2_(phi_pool_filtered).png", height = height, width = width)






#    2. B) bin the data, to stretch the variation
# The GC3 varation should be binned with mean, and phi_pool with median.
# About 10-20 bins for each





# Bin phi_pool
data_inliers_binned20 = data_inliers %>% filter(bin_size != 20000) %>% 
    group_by(bin_size, gs) %>% 
    mutate(GC3_bin = cut_number(phi_pool, 20)) %>% 
    
    
    group_by(GC3_bin, add = T) %>% 
    mutate(mean_GC3 = mean(GC3)) %>% # Add the mean GC3, because it looks prettier than bin-coordinates when plotting 
    
    
    group_by(mean_GC3, add = T) %>%
    summarise(median_phi_pool = median(phi_pool),
              median_ratio = median(ratio),
              median_c = median(c),
              median_phi_s = median(phi_s))
    

# Main binned plot
data_inliers_binned20 %>% ggplot(aes(mean_GC3, median_phi_s)) +
    geom_point() +
    scale_y_log10()+
    geom_smooth(method = "lm") +
    facet_grid(bin_size ~ gs, scales = "free") + 
    labs(title = "Inliers, 20 bins")
#ggsave("2B_inliers_binned20_(c).png", height = height, width = width)







