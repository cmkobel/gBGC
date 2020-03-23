# This script comparesthe recombination values inferred from mcorr and phi. Possibly later also ClonalFrame.


rm(list = ls())
library(tidyverse)
library(ggpmisc)
library(egg)

x_labels_rotated = theme(axis.text.x = element_text(angle = 90, hjust = 1))

setwd("~/urecomb/lokal/exports/export5_0C_structure/") # could also be on the cluster


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
        group_by(gene, genome) %>% 
                 summarise(GC3 = mean(gc_content))
}
gc_data_summarised = gather_gc()

# check distributions of GC
gc_data_summarised %>% 
    ggplot(aes(GC3)) + 
    geom_histogram() + 
    facet_wrap(~genome)
ggsave("~/genomedk/gBGC/carl/log/10A_0C_GC_distribution.png")


## ClonalFrameML
gather_cf = function(whatever = NULL) {
    cf_files = list.files(path=".", pattern="*clonalframe.tab", full.names=TRUE, recursive=T)
    cf_data = tibble()
    i = 1
    for (file in cf_files) {
        import <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("gene", "title", "parameter", "post_mean", "post_var", "a_post", "b_post"))
        print(paste(i, file, dim(import)[1]))
        cf_data = bind_rows(cf_data, import)
        rm(import)
        i = i + 1
    }
    
    #cf_data %>%
    cf_data %>% 
        mutate(genome = str_sub(title, 8, str_length(title)-2)) %>% 
        filter(parameter == "R/theta") #%>% select(-title)
}
cf_data = gather_cf()


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
        mutate(method = paste(method, detail)) %>% 
        select(-detail)  %>% 
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
    facet_grid(name~genome) + 
    geom_histogram()


#### Investigating difference between CF and PHI
# check the difference between CF and PHI
not_in_phi = anti_join(cf_data, phi_data, by = c("gene", "genome")) %>% transmute(gene, genome, not_in_phi = T)
not_in_cf = anti_join(phi_data, cf_data, by = c("gene", "genome")) %>% transmute(gene, genome, not_in_cf = T)



# Let's start with PHI
not_in_cf$not_in_cf %>% table
phi_gc_data = inner_join(phi_data, gc_data_summarised) %>%
    left_join(not_in_cf) %>% 
    mutate(not_in_cf = if_else(is.na(not_in_cf), F, T))

phi_gc_data %>%
    ggplot(aes(GC3, -log10(p_phi_permut+1e-10), color = not_in_cf)) + 
    geom_point(aes(size = not_in_cf), alpha = 0.3) +
    facet_wrap(~genome) 
ggsave("~/genomedk/gBGC/carl/log/16A_culprit_PHI_not_in_cf.png")

# I guess PHI seems quite honest with only 9 missing genes relative to CF

# And let's look at CF
not_in_phi$not_in_phi %>% table
cf_gc_data = inner_join(cf_data, gc_data_summarised) %>% 
    left_join(not_in_phi) %>% 
    mutate(not_in_phi = if_else(is.na(not_in_phi), F, T))

cf_gc_data %>% 
    ggplot(aes(GC3, log10(post_mean+1e-10), color = not_in_phi)) +
    geom_point(alpha = 0.1) +
    facet_wrap(~genome) + 
    labs(y = "log10(R/theta+1e-10)")
ggsave("~/genomedk/gBGC/carl/log/16B_culprit_CF_not_in_phi.png")






# Sanity check: Check that the relation betmeen PHI-bin p-values and GC3 is consistent.
alpha = 0.05
phi_data_thresh = phi_data %>% 
    group_by(genome) %>% 
    mutate(recombine = if_else(p_phi_permut < alpha/length(p_phi_permut), T, F), n_genes = length(p_phi_permut))


data = inner_join(gc_data_summarised, phi_data_thresh) 
data_binned20 = data %>% 
    group_by(genome) %>% 
    mutate(GC3_bin = cut_number(GC3, 20)) %>% 
    group_by(GC3_bin, add = T) %>%
    mutate(mean_GC3 = mean(GC3)) %>%  
    group_by(genome, mean_GC3) %>%  # use mean instead of weird bin-range
    count(recombine) %>% 
    spread(recombine, n, fill = 0) %>% 
    rename(n_recombining = `TRUE`, n_not_recombining = `FALSE`) %>% 
    mutate(ratio = n_recombining / (n_recombining + n_not_recombining), n = n_recombining + n_not_recombining)
    
data_binned20 %>% ggplot(aes(mean_GC3, n_recombining)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    stat_poly_eq(formula = y~x, 
                 aes(label = paste(..rr.label..)), 
                 parse = TRUE) +
    facet_wrap(~genome, scales = "free") 
ggsave("~/genomedk/gBGC/carl/log/13_phi_permut_bin.png")

###################################################################




#### Compare GC3 and clonalframe

cfgc_data = inner_join(cf_data, gc_data_summarised)


# Let's look at the parameter distributions

cfgc_data %>% pivot_longer(c(post_mean, post_var, a_post, b_post)) %>%
    ggplot(aes(value)) + 
    geom_histogram() + 
    facet_grid(genome~name, scales = "free")# +
    #scale_y_log10()
ggsave("~/genomedk/gBGC/carl/log/10_B_cf_parameter_distributions.png")


# plot directly without any preprocessing
cfgc_data %>% ggplot(aes(GC3, post_mean)) +
    geom_errorbar(aes(ymin = post_mean-sqrt(post_var), ymax = post_mean+sqrt(post_var)), alpha = 0.1) +
    geom_point(alpha = 0.2) +
    labs(y = "R/theta", caption = "Error bars: ± 1 SD") + 
    facet_wrap(~genome)
ggsave("~/genomedk/gBGC/carl/log/10_cf_raw_0C_structure.png")



# Let's bin, so we can compare to Lassalle

cfgc_data %>%
    group_by(genome) %>% 
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
    facet_wrap(~genome, scales = "free")
ggsave("~/genomedk/gBGC/carl/log/14_cf_20_bins_lm.png")





# What is the distribution of GC-content
cfgc_data %>% 
    ggplot(aes(GC3)) + 
    geom_histogram() + 
    facet_wrap(~genome)
# And peek at that from the PHI results
data %>% 
    ggplot(aes(GC3)) +
    geom_histogram() + 
    facet_wrap(~genome)

#



# Let's bin more explicitly with jitter
cfgc_data %>% 
    group_by(genome) %>% 
    mutate(GC3_bin = cut_interval(GC3, 20)) %>% 
    group_by(GC3_bin, add = T) %>% 
    mutate(mean_GC3 = mean(GC3),
           `median_R/theta` = median(post_mean)) %>%
    ungroup %>% 
    ggplot(aes(mean_GC3, `median_R/theta`)) +
    geom_jitter(aes(mean_GC3, post_mean), alpha = 0.05, width = 0.002, height = 0) +
    
    geom_point(aes(mean_GC3, `median_R/theta`), color = "blue") + 
    
    
    geom_smooth(color = "blue", method = "lm") + 
    stat_poly_eq(formula = y~x, 
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
                 parse = TRUE) +
    facet_wrap(~genome)
ggsave("~/genomedk/gBGC/carl/log/15_explicit_with_jitter_interval.png")

















#### Let's compare PHI and ClonalFrame


# compare distributions
data %>% select(post_mean, p_phi_permut, p_phi_normal) %>% 
    mutate(log(post_mean), log(p_phi_permut), log(p_phi_normal)) %>% 
    pivot_longer(everything()) %>% # everything() is all columns
    ggplot(aes((value))) + 
    geom_histogram()+
    facet_wrap(~name, scales = "free")

data = inner_join(cf_data, phi_data)
data %>% ggplot(aes(log(post_mean+1e-10), -log10(p_phi_normal+1e-10))) + 
    geom_point(alpha = 0.5) + 
    #labs(x = "R/theta") + 
    geom_smooth(method = "lm") + 
    geom_smooth(color = "red") +
    facet_wrap(~genospecies, scales = "free")
#ggsave("~/genomedk/gBGC/carl/log/9_PHIvsCF__.png")
                