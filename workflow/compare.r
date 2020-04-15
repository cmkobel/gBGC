# This script compares the recombination values inferred from mcorr and phi. Possibly later also ClonalFrame.



rm(list = ls())
library(tidyverse)
library(ggpmisc)
library(egg)
library(zoo)
#setwd("~/urecomb/lokal/exports/export4compare_binsize1/")
setwd('c:/Users/carl/urecomb/lokal/exports/export4compare_binsize1/')
setwd('c:/Users/carl/urecomb/lokal/exports/export7_cf_all_chromids/')


gff_data <- read_csv("C:/Users/carl/Desktop/xmfa2mfa/3206-3.gff") %>% select(-X1)
gff2_data <- read_delim("C:/Users/carl/Desktop/xmfa2mfa/Gene_function_pop_gene_data.csv", delim = ';')
gff_data = left_join(gff_data, gff2_data %>% select(`Gene group`, `Putative function`) %>% rename(gene_group = `Gene group`))

## GC
gather_gc = function(whatever = NULL) {
    gc_files = list.files(path=".", pattern="*bp_gc.tab", full.names=TRUE, recursive=F)
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
                 summarise(GC3 = mean(gc_content)) %>% 
        ungroup() %>% 
        mutate(gene = as.character(gene))
}
gc_data_summarised = gather_gc()

# check distributions of GC
gc_data_summarised %>% #filter(unitig == 0) %>% 
    ggplot(aes(GC3)) + 
    geom_histogram() + 
    facet_grid(unitig~genospecies, scales = "free_y")
ggsave('g:/gBGC/carl/log/38_all_gc_dist.png')




## ClonalFrameML
gather_cf = function(whatever = NULL) {
    cf_files = list.files(path=".", pattern="*clonalframe.tab", full.names=TRUE, recursive=F)
    print(paste('parsing from', cf_files))
    cf_data = tibble()
    i = 1
    for (file in cf_files) {
        import <- read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("gene", "title", "parameter", "post_mean", "post_var", "a_post", "b_post"))
        print(paste(i, file, dim(import)[1]))
        cf_data = bind_rows(cf_data, import)
        rm(import)
        i = i + 1
    }

    cf_data %>%
        mutate(unitig = str_sub(title, 15, 15),
               genospecies = str_sub(title, 31, 31)) %>%
        filter(parameter == "R/theta") %>% select(-title)
}
cf_data = gather_cf()
cf_data$unitig %>% table; cf_data$genospecies %>% table
cf_data %>% ggplot(aes(post_mean)) +
    geom_histogram() +
    facet_grid(unitig~genospecies, scales = 'free_y')


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
    facet_wrap(~genospecies, scales = "free")
ggsave("~/genomedk/gBGC/carl/log/8_cd_20_bins_lm.png")



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
ggsave("~/genomedk/gBGC/carl/log/9_PHIvsCF__.png")
                










#### Let's look at the distribution of GC3 and recombination over the genome.

# import annotaton 
gff_data$plasmid %>% table
cf_gc_annot = gff_data %>% #filter(plasmid == '3206-3_scaf_3_chromosome-00') %>%
    group_by(plasmid) %>% 
    arrange(start) %>%
    mutate(length = abs(end - start),
           mid = end- (length/2),
           gene = str_sub(gene_group, 6)) %>% 
    inner_join(cf_data %>% mutate(gene = as.character(gene))) %>%   # add CF
    inner_join(gc_data_summarised) %>%                              # add gc
    # add quantile information
    group_by(unitig, genospecies) %>% 
    arrange(post_mean) %>% 
    mutate(rown = row_number(), rank = rown/(length(GC3)))


perimeter = cf_gc_annot %>%
    filter(genospecies == "C") %>%
    filter(unitig == '1')
# Create simple plots for each chromid

X = perimeter %>% 
    ggplot(aes(x = mid)) +
    geom_point(aes(y = post_mean))

Y = perimeter %>% 
    ggplot(aes(x = mid)) +
    geom_point(aes(y = GC3))

ggarrange(X, Y, ncol = 1)




# Compute threshold: #TODO: Try with the median instead.
variable_GC3_threshold = mean(cf_gc_annot %>% filter(genospecies == genospecies) %>% pull(GC3) %>% .^1)
variable_GC3_threshold = cf_gc_annot$GC3 %>% median


# Mark genes above threshold
cf_gc_annot = cf_gc_annot %>%
    mutate(GC3_threshold = if_else(GC3 > variable_GC3_threshold, T, F))




genospecies = 'C'

A = cf_gc_annot %>% filter(genospecies == genospecies) %>%
    ggplot(aes(mid, (post_mean), color = GC3_threshold)) + 
    geom_point(alpha = 0.5, size = 1.2) +
    #geom_errorbar(aes(ymin = post_mean - sqrt(post_var),
    #                  ymax = post_mean + sqrt(post_var)), alpha = 0.15) + 
    labs(#title = paste('Genospecies', genospecies),
         #caption = 'error bars: SD',
         x = '',
         y = 'R/theta') +
    theme(legend.position="none")
    #guides(color=guide_legend(title="Above median GC3"))
A

A_hist = cf_gc_annot %>% filter(genospecies == genospecies) %>% 
    ggplot(aes((post_mean), fill = GC3_threshold)) +
    geom_histogram() +
    theme(legend.position="none") 
    #scale_y_log10()
A_hist


B = cf_gc_annot %>% filter(genospecies == genospecies) %>%
    ggplot(aes(mid, GC3, color = GC3_threshold)) + 
    geom_point(alpha = 0.5, size = 1) +
    labs(title = paste('Genospecies', genospecies), x = 'unitig 0') + 
    geom_hline(yintercept = variable_GC3_threshold) +
    theme(legend.position="none") +
    scale_y_reverse()
B

B_hist = cf_gc_annot %>% filter(genospecies == genospecies) %>% 
    ggplot(aes(GC3, fill = GC3_threshold)) +
    geom_histogram() +
    geom_vline(xintercept = variable_GC3_threshold) +
    theme(legend.position = 'none') +
    scale_x_reverse()
B_hist

ratio = .85
arr = ggarrange(B, B_hist, A, A_hist, ncol = 2, widths = c(ratio, 1-ratio))
#ggarrange(A, B, ncol = 1)
arr

ggsave('c:/Users/carl/repositories/gBGC/log/33_cf_gc.png', plot = arr, height = 8, width = 30, dpi = 400)
ggsave('c:/Users/carl/repositories/gBGC/log/33_cf_gc_small.png', plot = arr)

# Find the top genes:
top1 = cf_gc_annot %>% group_by(genospecies, unitig) %>% 
    filter(post_mean > quantile(post_mean, .99)) %>% 
    ungroup %>% arrange(genospecies, desc(post_mean)) %>% 
    select(genospecies, unitig, gene_group, post_mean, post_var, start, end, length, product, `Putative function`) %>% 
    rename(`R/theta` = post_mean,
           `Var(R/theta)` = post_var,
           `Product (3206-3.gff)` = product,
           `Putative function (pop_gene_data.csv)` = `Putative function`)

write_tsv(top1, 'c:/Users/carl/repositories/gBGC/log/top1.tsv')


    


# Sliding window of many genes.
# also, use three classes of GC3 (for coloring)

selected_genospecies = 'C'
number_of_bins = 50
midpoint_ = mean(cf_gc_annot %>% filter(genospecies == selected_genospecies) %>% pull(GC3) %>% mean)

cf_gc_annot %>%
    ungroup() %>% 
    group_by(unitig, genospecies) %>% 
    mutate(position_bin = cut_number(mid, number_of_bins)) %>% 
    
    group_by(position_bin, add = T) %>% 
    mutate(binned_position = mean(`mid`),
           median_post_mean = median(post_mean),
           mean_GC3 = mean(GC3)) %>% 
    
    ungroup() %>% 
    filter(genospecies == selected_genospecies) %>% 
    arrange(binned_position)  %>% 
    ggplot(aes(binned_position, (median_post_mean), color = mean_GC3)) + 
    geom_line(color = "black", alpha = 0.15)  +
    geom_point(size = 2) +
    #geom_errorbar(aes(ymin = post_mean - sqrt(post_var),
    #                  ymax = post_mean + sqrt(post_var)), alpha = 0.15) + 
    #labs(#title = paste('Genospecies', genospecies),
    #    caption = paste('number of bins:', number_of_bins),#'error bars: SD',
    #    x = '',
    #    y = 'R/theta')  +
    scale_color_gradient2(midpoint = .62, low = "red1", mid = "grey", high = "green4", space = "Lab")
    #theme(legend.position="none")

scaleFUN <- function(x) sprintf("%.2f", x) 

truncFUN <- function(x) x - (x%%10)



annotation_custom2 <- 
    function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
    {
        layer(data = data, stat = StatIdentity, position = PositionIdentity, 
              geom = ggplot2:::GeomCustomAnn,
              inherit.aes = TRUE, params = list(grob = grob, 
                                                xmin = xmin, xmax = xmax, 
                                                ymin = ymin, ymax = ymax))
    }

# 35
# R/theta (color GC3)
# roll instead
i = 1
#p = 0; for (i in c(1, 5, seq(10, 1000, 10), seq(990, 10, -10), 5)) {
for (i in c(1, 5, seq(10, 1000, 10))) {
    p = p + 1
    print(i)
    roll_width = i
    plot_data = cf_gc_annot %>% 
        group_by(unitig, genospecies) %>% 
        arrange(mid) %>% 
        mutate(roll_post_mean = rollapply(post_mean, roll_width, median, fill = NA)) %>%  #View
        mutate(roll_GC3 = rollapply(GC3, roll_width, mean, fill = NA)) %>%
        filter(genospecies == 'C') %>%  #View
        rename(`R/theta` = roll_post_mean,
               `rolling window mean GC3` = roll_GC3)
    plot_data %>% 
        ggplot(aes(mid, `R/theta`+1e-10, color = `rolling window mean GC3`)) + 
        #geom_point() +
        geom_line(size = 1) +
        #geom_vline(xintercept = min(plot_data$mid), color = 'grey') +
        #geom_vline(xintercept = max(plot_data$mid), color = 'grey') +
        

                      
        scale_color_gradientn(colours = c('red1', 'grey', 'green4')) +
        labs(#caption = 'green: high GC3\ngrey: median GC3\nred: low GC3',
            x = 'chromosome position',
            y = 'rolling window median R/theta',
            title = 'Genospecies C',
            subtitle = paste('width', roll_width))  +
        #theme(legend.position = 'none') +
        scale_y_log10(#breaks = seq(0, 4, .02), 
                      labels=scaleFUN) +
        theme(panel.grid.minor = element_blank()) 
    #guides(color=guide_legend(title="Rolling window GC3"))
    #scale_color_continuous(breaks = seq(0.00, 1.00, 0.1)) +
    #ylim(c(0.15, 1.2))
    ggsave(paste0('c:/Users/carl/repositories/gBGC/log/34/34_rolling_recomb_GC3_', str_pad(roll_width, 4, pad = 0), '.png'))
}


#37
# GC3 (color R/theta)
# roll instead
i = 100
#p = 0; for (i in c(1, 5, seq(10, 1000, 10), seq(990, 10, -10), 5)) {
for (i in c(1, 5, seq(10, 1000, 10))) {
        
    print(i)
    p = p + 1
    roll_width = i
    cf_gc_annot %>% 
        group_by(unitig, genospecies) %>% 
        arrange(mid) %>% 
        mutate(roll_post_mean = rollapply(post_mean, roll_width, median, fill = NA)) %>%  #View
        mutate(roll_GC3 = rollapply(GC3, roll_width, mean, fill = NA)) %>%
        filter(genospecies == 'C') %>%  #View
        rename(`R/theta` = roll_post_mean) %>% 
        
        ggplot(aes(mid, roll_GC3, color = log10(`R/theta`+1e-10))) + 
        #geom_point() +
        geom_line(size = 1) +
        
        scale_color_gradientn(colours = c('red1', 'grey', 'green4')) +
        labs(#caption = 'green: high GC3\ngrey: median GC3\nred: low GC3',
            x = 'chromosome position',
            y = 'rolling window mean GC3',
            title = 'Genospecies C',
            subtitle = paste('width', roll_width)) +
        #theme(legend.position = 'none') +
        scale_y_continuous(breaks = seq(0, 4, .01), labels=scaleFUN)
    #theme(panel.grid.minor = element_blank()) 
    #guides(color=guide_legend(title="log10(R/theta + 1e-10)"))
    #scale_color_continuous(breaks = seq(0.00, 1.00, 0.1)) +
    #ylim(c(0.15, 1.2))
    ggsave(paste0('c:/Users/carl/repositories/gBGC/log/36/36_rolling_GC3_recomb_', str_pad(roll_width, 4, pad = 0), '.png'))
}
           