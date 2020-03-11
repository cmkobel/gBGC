# This script compares the recombination values inferred from mcorr and phi. Possibly later also ClonalFrame.


rm(list = ls())
library(tidyverse)
æs = aes
setwd("~/urecomb/lokal/exports/export4compare/")

'
## GC
gather_gc = function(whatever = NULL) {
    gc_files = list.files(path=".", pattern="*bp_gc.tab", full.names=TRUE, recursive=T)
    gc_data = tibble()
    i = 1
    for (file in gc_files) {
        import <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE)
        print(paste(i, file, dim(import)[1]))
        gc_data = bind_rows(gc_data, import)
        rm(import)
        i = i + 1
    }
    #rm(gc_files, gc_data, file, i) # clean up
    gc_data %>% group_by(gene, bin_size = cli_comment_1, genome = cli_comment_2) %>% 
        summarise(GC3 = mean(gc_content))
}
gc_data_summarised = gather_gc()
'

## mcorr    import fitted (recombination) parameters
gather_mcorr = function() {
    mcorr_files = list.files(path=".", pattern="*_fitpars.csv", full.names=TRUE, recursive=T)
    mcorr_colnames = c("gene", "group", "d_sample", "theta_pool", "phi_pool", "ratio", "fbar", "c", "d_pool", "d_clonal", "theta_s", "phi_s", "bin_size", "genome")
    
    mcorr_data = tibble()
    i  = 1
    for (file in mcorr_files) {
        import = read_csv(file, col_names = mcorr_colnames)
        print(paste(i, file, dim(import)[1]))
        mcorr_data = bind_rows(mcorr_data, import)

        i = i + 1
    }
    mcorr_data %>%
        rename(mcorr_ratio = ratio,
               mcorr_c = c) %>% 
        select(-group)
        
}
mcorr_data = gather_mcorr()


## PHI files
gather_PHI = function() {
    phi_files = list.files(path=".", pattern="*phi_results.tab", full.names=TRUE, recursive=T)
    phi_data = tibble()
    i = 1
    for (file in phi_files) {
        import = read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE, na = c("--"))
        print(paste(i, file, dim(import)[1]))
        phi_data = bind_rows(phi_data, import)
        rm(import)
        i = i + 1
    }
    phi_data %>%
        mutate(method = paste(method, detail),
               unitig = str_sub(genome, 8, 8),
               genospecies = str_sub(genome, 24, 24)) %>% 
        select(-detail) %>% 
        spread(method, pvalue) %>% 
        rename(p_maxchisq = 'Max Chi^2:',
               p_nss = 'NSS NA:',
               p_phi_normal = 'PHI (Normal):',
               p_phi_permut = 'PHI (Permutation):')
}
phi_data = gather_PHI() 



# It doen't make sense to bin, because we don't care about GC.


data = inner_join(mcorr_data, phi_data) %>% select(-genome)

# main plot
data %>% ggplot(æs(phi_pool, p_phi_permut)) + 
    geom_point()






## Let's filter the mcorr tails away
tail_size_phi_pool = 0.1
tail_size_ratio = 0.01
tail_size_c = 0.00
tail_size_phi_s = 0.02

data_inliers = data %>%
    group_by(bin_size, genospecies, unitig) %>% 
    filter(phi_pool >= quantile(phi_pool, tail_size_phi_pool) & phi_pool <= quantile(phi_pool, 1-tail_size_phi_pool)) %>% 
    filter(mcorr_ratio >= quantile(mcorr_ratio, tail_size_ratio) & mcorr_ratio <= quantile(mcorr_ratio, 1-tail_size_ratio)) %>% 
    filter(mcorr_c >= quantile(mcorr_c, tail_size_c) & mcorr_c <= quantile(mcorr_c, 1-tail_size_c)) %>% 
    filter(phi_s >= quantile(phi_s, tail_size_c) & phi_s <= quantile(phi_s, 1-tail_size_c)) %>% 
    
    ungroup()


# sanity check of tail-removal
# various histograms
data_inliers %>%
    filter(bin_size == 20000) %>% 
    #pivot_longer(c(phi_pool, mcorr_ratio, mcorr_c, phi_s)) %>% 
    ggplot(aes(phi_pool, fill = genospecies)) +
    geom_histogram() + 
    facet_grid(name~genospecies)


# main plot
data_inliers %>%
    #filter(bin_size == 30000) %>% 
    ggplot(aes(GC3, log10(phi_pool), color = gs)) +
    geom_point(alpha = 0.5) + 
    geom_smooth() +
    facet_grid(bin_size~gs)

