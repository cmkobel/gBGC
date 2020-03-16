# This script compares the recombination values inferred from mcorr and phi. Possibly later also ClonalFrame.


rm(list = ls())
library(tidyverse)
library(ggpmisc)
library(egg)
æs = aes
setwd("~/urecomb/lokal/exports/export4compare_binsize1/")


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
    #gc_data %>%
    #    mutate(unitig = str_sub(genome, 8, 8),
    #           genospecies = str_sub(genome, 24, 24)) %>% 
    #    group_by(gene, bin_size = cli_comment_1, genome = cli_comment_2) %>% 
    #    summarise(GC3 = mean(gc_content))
    gc_data %>% 
        rename(bin_size = cli_comment_1, genome = cli_comment_2) %>% 
        mutate(unitig = str_sub(genome, 8, 8),
               genospecies = str_sub(genome, 24, 24)) %>% 
        group_by(gene, unitig, genospecies, bin_size, genome) %>% 
                 summarise(GC3 = mean(gc_content))
}
gc_data_summarised = gather_gc()

# check distributions of GC
gc_data_summarised %>% filter(unitig == 0) %>% 
    ggplot(aes(GC3)) + 
    geom_histogram() + 
    facet_wrap(~genospecies)


    ## ClonalFrameML
gather_cf = function(whatever = NULL) {
    cf_files = list.files(path="filtered", pattern="*clonalframe.tab", full.names=TRUE, recursive=T)
    cf_data = tibble()
    i = 1
    for (file in cf_files) {
        import <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("gene", "title", "parameter", "post_mean", "post_var", "a_post", "b_post"))
        print(paste(i, file, dim(import)[1]))
        cf_data = bind_rows(cf_data, import)
        rm(import)
        i = i + 1
    }
    #rm(cf_files, cf_data, file, i) # clean up
    #cf_data %>% group_by(gene, bin_size = cli_comment_1, genome = cli_comment_2) %>% 
    #    summarise(GC3 = mean(gc_content))
    cf_data %>%
        mutate(unitig = str_sub(title, 15, 15),
               genospecies = str_sub(title, 31, 31)) %>%
        filter(parameter == "R/theta") %>% select(-title)
}
cf_data = gather_cf()
cf_data$genospecies %>% table


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
    alpha = 0.5
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
        group_by(method, unitig, genospecies) %>% 
   #alpha/length(pvalue), T, F), n_genes = length(pvalue))
        spread(method, pvalue) %>% 
        rename(p_maxchisq = 'Max Chi^2:',
               p_nss = 'NSS NA:',
               p_phi_normal = 'PHI (Normal):',
               p_phi_permut = 'PHI (Permutation):')
        
}
phi_data = gather_PHI() 


# Check distribution of p-values from PHI
phi_data %>% pivot_longer(starts_with("p_")) %>% 
    ggplot(aes(value)) + 
    facet_grid(name~genospecies) + 
    geom_histogram()
ggsave("~/genomedk/gBGC/carl/log/6_pvals_forbinsize1.png")


# It doen't make sense to bin, because we don't care about GC.


data = inner_join(inner_join(gc_data_summarised, phi_data), mcorr_data) %>% select(-genome)


## Let's check that phi measures OK when concatenating the genes. (Remember, we have to concatenate the genes in order to be able to measure anything).

# We have to drop na's, or many bins will fall out.
data = data %>% drop_na()

# bin the GC content
'
data_binned20 = data %>% 
    group_by(genospecies, unitig) %>% 
    mutate(GC3_bin = cut_number(GC3, 20)) %>% 
    group_by(GC3_bin, add = T) %>%
    mutate(mean_GC3 = mean(GC3)) %>% # View
    group_by(genospecies, unitig, mean_GC3) %>%  # use mean instead of weird bin-range
    count(recombine) %>%
    spread(recombine, n, fill = 0) %>% 
    rename(n_recombining = `TRUE`, n_not_recombining = `FALSE`) %>% 
    mutate(ratio = n_recombining / (n_recombining + n_not_recombining))
'

# main plot
data %>% ggplot(aes(log10(phi_pool), (p_phi_normal))) + 
    geom_point() + 
    facet_wrap(~genospecies, scales = "free")+ 
    geom_smooth()
ggsave("~/genomedk/gBGC/carl/log/7_main_compare_logphipool_phinormal_free.png")






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
    #filter(bin_size == 20000) %>% 
    pivot_longer(c(phi_pool, mcorr_ratio, mcorr_c, phi_s)) %>% 
    ggplot(aes(value, fill = genospecies)) +
    geom_histogram() + 
    facet_grid(name~genospecies) + 
    labs(caption = "Only unitig 0 and bin size 20K is included in this dataset.")


# main plot
# mcorr_phi_pool and PHI
data %>% 
    #filter(bin_size == 30000) %>% 
    ggplot(aes(p_phi_permut, (mcorr_ratio), color = genospecies)) +
    geom_point(alpha = 0.5) + 
    geom_smooth(formula = y~x, method = "lm") +
    facet_grid(.~genospecies, scales = "free") +
    theme(axis.text.x = element_text(vjust = 0.5, angle = -90)) +
    labs(title = "")
ggsave("main_compare.pdf")



# Try binning by GC content then compare PHI ratio and median phi_pool





# Sanity check: Check that the relation betmeen PHI-bin p-values and GC3 is consistent.
alpha = 0.05
phi_data_thresh = phi_data %>% 
    group_by(unitig, genospecies) %>% 
    mutate(recombine = if_else(p_phi_permut < alpha/length(p_phi_permut), T, F), n_genes = length(p_phi_permut))


data = inner_join(gc_data_summarised, phi_data_thresh) %>% select(-genome)
data_binned20 = data %>% 
    group_by(genospecies, unitig) %>% 
    mutate(GC3_bin = cut_number(GC3, 20)) %>% 
    group_by(GC3_bin, add = T) %>%
    mutate(mean_GC3 = mean(GC3)) %>%  
    group_by(genospecies, unitig, mean_GC3) %>%  # use mean instead of weird bin-range
    count(recombine) %>% 
    spread(recombine, n, fill = 0) %>% 
    rename(n_recombining = `TRUE`, n_not_recombining = `FALSE`) %>% 
    mutate(ratio = n_recombining / (n_recombining + n_not_recombining), n = n_recombining + n_not_recombining)
    
data_binned20 %>% ggplot(aes(mean_GC3, n_recombining)) + 
    geom_point() + 
    facet_wrap(~genospecies, scales = "free") + 
    geom_smooth(method = "lm")






#### Compare GC3 and clonalframe

cfgc_data = inner_join(cf_data, gc_data_summarised)



# plot directly without any preprocessing
cfgc_data %>% ggplot(aes(GC3, post_mean)) +
    geom_errorbar(aes(ymin = post_mean-sqrt(post_var), ymax = post_mean+sqrt(post_var)), alpha = 0.1) +
    geom_point(alpha = 0.2) +
    labs(y = "R/theta", caption = "Error bars: ± 1 SD") + 
    facet_wrap(~genospecies)
ggsave("~/genomedk/gBGC/carl/log/8_cf_raw.png")


# Let's look at the parameter distributions

cfgc_data %>% pivot_longer(c(post_mean, post_var, a_post, b_post)) %>%
    ggplot(aes(value)) + 
    geom_histogram() + 
    facet_grid(genospecies~name, scales = "free")
ggsave("~/genomedk/gBGC/carl/log/8_cf_parameter_distributions.png")


# Let's bin, so we can compare to Lassalle

cfgc_data %>%
    group_by(unitig, genospecies) %>% 
    mutate(GC3_bin = cut_number(GC3, 20)) %>%
    group_by(GC3_bin, add = T) %>% 
    mutate(mean_GC3 = mean(GC3),
           `median_R/theta` = median(post_mean)) %>% 
    ggplot(aes(mean_GC3, `median_R/theta`)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    stat_poly_eq(formula = y~x, 
                 aes(label = paste(..rr.label..)), 
                 parse = TRUE) +
    facet_wrap(~genospecies)
ggsave("~/genomedk/gBGC/carl/log/8_cd_20_bins_lm.png")



#### Let's compare PHI and ClonalFrame


data = inner_join(cf_data, phi_data)
data %>% ggplot(aes(post_mean, -log(p_phi_normal))) + 
    geom_point() + 
    labs(x = "R/theta") + 
    geom_smooth(method = "lm") + 
    geom_smooth(color = "red")
ggsave("~/genomedk/gBGC/carl/log/9_PHIvsCF.png")
                