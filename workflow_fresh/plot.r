library(tidyverse)
library(zoo)
library(gganimate)
library(ggpmisc)
#library(ggpubr)
library(cowplot)
setwd("c:/users/carl/urecomb/lokal")

height = 5; width = 8


# ClonalFrameML files
CF_raw = read_delim("G:/gBGC/carl/workflow_fresh/output/CF_collected.tab", 
                           "\t", escape_double = FALSE, trim_ws = TRUE,
                 col_names = c('genospecies', 'unitig', 'file', 'parameter', 'post_mean', 'post_var', 'a_post', 'b_post')) %>% 
    mutate(gene = str_sub(file, 4, str_length(file)-7)) %>% select(-file)


#CF = CF_raw %>% pivot_wider(names_from = parameter, values_from = c(post_mean, post_var, a_post, b_post))
CF = CF_raw %>% filter(parameter == "R/theta") %>% select(-parameter) 



# GC3 files
GC_raw = read_delim("G:/gBGC/carl/workflow_fresh/output/GC3_correct.tab", 
                  "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
    mutate(sample = str_sub(header, 6),
           gene = str_sub(file, 4, str_length(file)-13)) %>% 
    select(-header, -file)


GC3 = GC_raw %>% group_by(genospecies, unitig, gene) %>% 
    summarize(mean_GC = mean(GC), 
              mean_GC1 = mean(GC1), 
              mean_GC2 = mean(GC2), 
              mean_GC3 = mean(GC3), 
              n_samples_GC3 = length(GC3)) %>%
    ungroup %>%
    rename(GC = mean_GC,
           GC1 = mean_GC1,
           GC2 = mean_GC2,
           GC3 = mean_GC3) 

GC_all = GC_raw %>% group_by(genospecies, unitig, gene) %>% 
    summarize(GC = mean(GC),
              GC1 = mean(GC1),
              GC2 = mean(GC2),
              GC3 = mean(GC3)) %>% 
    ungroup

x_min = 0.20
x_max = .90

# Visualize GC1, 2, 3 distributions
A1 = GC_all %>% #pivot_longer(starts_with("GC"), names_to = "positions", values_to = "proportion GC") %>% 
    #filter(positions != "GC") %>% 
    ggplot(aes(GC1)) +
    geom_histogram(fill = "red") +
    #geom_histogram(aes(mean_GC1), alpha = 0.4, fill = "red") +
    #geom_histogram(aes(mean_GC2), alpha = 0.4, fill = "green") +
    #geom_histogram(aes(mean_GC3), alpha = 0.4, fill = "blue")
    
    #facet_grid(.~positions) +
    theme_light() + 
    labs(x = "") +
    coord_flip()  + 
    theme(legend.position = "none") + 
    xlim(c(x_min, x_max))
A1
A2 = GC_all %>% #pivot_longer(starts_with("GC"), names_to = "positions", values_to = "proportion GC") %>% 
    #filter(positions != "GC") %>% 
    ggplot(aes(GC2)) +
    geom_histogram(fill = "green") +
    #geom_histogram(aes(mean_GC1), alpha = 0.4, fill = "red") +
    #geom_histogram(aes(mean_GC2), alpha = 0.4, fill = "green") +
    #geom_histogram(aes(mean_GC3), alpha = 0.4, fill = "blue")
    
    #facet_grid(.~positions) +
    theme_light() + 
    labs(x = "") +
    coord_flip()  + 
    theme(legend.position = "none") +
    xlim(c(x_min, x_max))
A2
A3 = GC_all %>% #pivot_longer(starts_with("GC"), names_to = "positions", values_to = "proportion GC") %>% 
    #filter(positions != "GC") %>% 
    ggplot(aes(GC3)) +
    geom_histogram(fill = "blue") +
    #geom_histogram(aes(mean_GC1), alpha = 0.4, fill = "red") +
    #geom_histogram(aes(mean_GC2), alpha = 0.4, fill = "green") +
    #geom_histogram(aes(mean_GC3), alpha = 0.4, fill = "blue")
    
    #facet_grid(.~positions) +
    theme_light() + 
    labs(x = "") +
    coord_flip()  + 
    theme(legend.position = "none") +
    xlim(c(x_min, x_max))
A3

AA = GC_all %>% pivot_longer(starts_with("GC"), names_to = "position", values_to = "proportion GC") %>%# View
    #mutate(position = factor(position, levels=c("GC", "GC3", "GC1", "GC2"))) %>% 
    filter(position != "GC") %>% 
    ggplot(aes(`proportion GC`, fill = position)) +
    geom_histogram(position = "dodge") +

    
    facet_grid(.~position) +
    theme_light() + 
    labs(x = "") +
    #theme(axis.text.x = element_text(angle = -90, vjust = .4)) +
    coord_flip()  + 
    #theme(legend.position = "none") +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ) +
    theme(axis.text.x = element_text(angle = -90, vjust = .3)) +
    xlim(c(x_min, x_max)) 
AA





B = GC_all %>% pivot_longer(c(GC1, GC2, GC3), names_to = "positions", values_to = "GCn") %>% 
    ggplot(aes(GC, GCn, color = positions)) + 
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm") +
    theme_light() + 
    theme(legend.position = "none") + 
    ylim(c(x_min, x_max)) +
    theme(axis.text.x = element_text(angle = -90, vjust = .5))
B

#ggarrange(B, ggarrange(A1, A2, A3, ncol = 3), ncol = 2)
#ggarrange(B, AA, ncol = 2, labels = c("A", "B"))
plot_grid(BB, A, labels = c("A", "B"), ncol = 2, align = "h", axis = "bt")

#ggsave("final_plots/Y.png", height = height, width = width)



   
# extract statistics from the GC3 measurement
GC_raw %>% group_by(genospecies, unitig) %>% 
    summarize(n_genes = length(unique(gene))) %>% 
    spread(unitig, n_genes) %>% View


## Investigate distribution of GC3
d_mean = GC3$GC3 %>% mean
d_sd = GC3$GC3 %>% sd
d_seq = seq(0, 1, 0.01)
d_d = dnorm(d_seq, mean = d_mean, sd = d_sd)*203.41

ggplot() + 
    geom_histogram(
        data = GC3, 
        mapping = aes(GC3), 
        bins = 100) +
    geom_histogram(
        data = GC3, 
        mapping = aes(GC3), 
        bins = 100) +
    
    
    geom_line(
        data = tibble(x = d_seq, y = d_d), 
        mapping = aes(x,y), linetype = "dashed") +
    geom_area(
        data = tibble(x = d_seq, y = d_d), 
        mapping = aes(x,y), linetype = "dashed", fill = "blue", alpha = 0.15) + 
    theme_light()





## PHI files
gather_PHI = function() {
    alpha = 0.5
    phi_files = list.files(path="c:/users/carl/urecomb/lokal/exports/export7_cf_all_chromids", pattern="*phi_results.tab", full.names=TRUE, recursive=T)
    phi_data = tibble()
    i = 1
    file = phi_files[1]
    for (file in phi_files) {
        import = read_delim(file, "\t", escape_double = FALSE, trim_ws = TRUE, na = c("--")) %>% 
            mutate(genospecies = str_sub(basename(file), 24, 24),
                   unitig = str_sub(basename(file), 8,8),
                   method = paste(method, detail)) %>% 
            select(-detail, -genome)# %>% filter(method == "PHI (Permutation):")
        print(paste(i, file, dim(import)[1]))
        phi_data = bind_rows(phi_data, import)
        rm(import)
        i = i + 1
    }
    phi_data %>% 
        # mutate(method = paste(method, detail),
        #        unitig = str_sub(genome, 8, 8),
        #        genospecies = str_sub(genome, 24, 24)) %>% 
        # select(-detail) %>% 
        # group_by(method, unitig, genospecies) %>% 
        # #alpha/length(pvalue), T, F), n_genes = length(pvalue))
        
        spread(method, pvalue) %>%
        rename(infsites = `infsites `,
               p_maxchisq = 'Max Chi^2:',
               p_nss = 'NSS NA:',
               p_phi_normal = 'PHI (Normal):',
               p_phi_permut = 'PHI (Permutation):')
    
}
phi_data = gather_PHI() %>% select(genospecies, unitig, infsites, gene, p_phi_permut) %>%  mutate(gene = paste0("group", gene))  %>% mutate(unitig = as.numeric(unitig))




gff_data <- read_csv("C:/Users/carl/Desktop/xmfa2mfa/3206-3.gff") %>% select(-X1) %>% mutate(mid = start + ((end-start)/2), length = end-start)
gff2_data <- read_delim("C:/Users/carl/Desktop/xmfa2mfa/Gene_function_pop_gene_data.csv", delim = ';')

# Before I delete gff2_data, I want to investigate a few things.
gff2_data %>% pull

gff_data = left_join(gff_data, gff2_data %>% select(`Gene group`, `Putative function`) %>% rename(gene_group = `Gene group`)) %>% rename(gene = gene_group)
rm(gff2_data)

### Join all data
data = inner_join(GC3, CF) %>% inner_join(phi_data) %>% inner_join(gff_data)
#saveRDS(data, "main_correctGC3_gcall.rds")
data = readRDS("main_correctGC3_gcall.rds")


### informative sites information
### Maybe informative sites can explain R^2
data %>% group_by(genospecies, unitig) %>% 
    ggplot(aes(infsites)) + geom_histogram() +
    facet_grid(genospecies ~ unitig)

infsites_info = data %>% group_by(genospecies, unitig) %>% 
    summarize(sum_infsites = sum(infsites)) %>% 
    mutate(p = paste(genospecies, unitig))

infsites_info %>%     
    ggplot(aes(p, sum_infsites, fill = genospecies)) + 
    geom_col()

infsites_info %>% filter(unitig == 0)
# Groups:   genospecies [5]
# genospecies unitig sum_infsites p    
# <chr>        <dbl>        <dbl> <chr>
# 1 A                0       145718 A 0  
# 2 B                0        41606 B 0  
# 3 C                0       178092 C 0  
# 4 D                0        13112 D 0  
# 5 E                0        51992 E 0 



### Create a plot that shows the distribution PHI
data %>% ggplot(aes(p_phi_permut)) +
    geom_histogram() +
    facet_grid(plasmid ~ genospecies)


# Create a plot that shows the distribution of CFML
data %>% ggplot(aes(post_mean)) +
    geom_histogram() +
    facet_grid(plasmid ~ genospecies)



### Create a plot that shows how PHI and CF correlate
data %>% ggplot(aes(log(post_mean), -p_phi_permut)) + 
    geom_point(alpha = 0.5) + 
    facet_grid(genospecies ~ plasmid) + 
    geom_smooth()
ggsave("g:/gBGC/carl/log/50_A.png")






### Create a plot emulates what lassale did with 20 bins. 
### Here it is grouped by plasmids
data %>% ungroup %>% select(genospecies, plasmid, gene, p_phi_permut, GC3) %>% 
    group_by(genospecies, plasmid) %>% 
    mutate(sig_rec = if_else(p_phi_permut < 0.05/length(p_phi_permut), 1, 0),
           GC3_bin = cut_number(GC3, 20)) %>% 
    group_by(GC3_bin, add = T) %>% 
    mutate(n_sig_rec = sum(sig_rec),
           mean_GC3 = mean(GC3)) %>% 
    
    ggplot(aes(mean_GC3, n_sig_rec)) + 
    geom_point() + 
    facet_grid(plasmid ~ genospecies, scales = "free") + 
    geom_smooth(method = "lm") + 
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..rr.label.., sep = "~~~")), 
                 parse = TRUE)


### This one is the same, but it is grouped by unitig.
data %>% ungroup %>%
    filter(unitig == 0) %>% 
    select(genospecies, unitig, gene, p_phi_permut, GC3) %>% 
    group_by(genospecies, unitig) %>% 
    mutate(sig_rec = if_else(p_phi_permut < 0.05/length(p_phi_permut), 1, 0),
           GC3_bin = cut_number(GC3, 20)) %>% 
    group_by(GC3_bin, add = T) %>% 
    mutate(n_sig_rec = sum(sig_rec),
           mean_GC3 = mean(GC3)) %>% 
    
    ggplot(aes(mean_GC3, n_sig_rec)) + 
    geom_point() + 
    facet_wrap(~ genospecies, scales = "free") + 
    geom_smooth(method = "lm") + 
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..rr.label.., sep = "~~~")), 
                 parse = TRUE)
ggsave("g:/gBGC/carl/log/51_B_.png")





### Now, create a Lassalle-ish plot, but with CF instead of PHI.
data %>% ungroup %>%
    filter(unitig == 0) %>%
    #filter(plasmid %in% c("3206-3_scaf_1_chromosome-01", "3206-3_scaf_2_chromosome-02", "3206-3_scaf_3_chromosome-00")) %>% 
    select(genospecies, unitig, gene, post_mean, GC3) %>% 
    group_by(genospecies, unitig) %>% 
    mutate(#sig_rec = if_else(p_phi_permut < 0.05/length(p_phi_permut), 1, 0),
           GC3_bin = cut_number(GC3, 20)) %>% 
    group_by(GC3_bin, add = T) %>% 
    mutate(#n_sig_rec = sum(sig_rec),
           median_post_mean = median(post_mean),
           mean_GC3 = mean(GC3)) %>% 
    
    ggplot(aes(mean_GC3, median_post_mean)) + 
    geom_point() + 
    facet_wrap(~ genospecies, scales = "free") + 
    geom_smooth(method = "lm") + 
    stat_poly_eq(formula = y ~ x, 
                 aes(label = paste(..rr.label.., sep = "~~~")), 
                 parse = TRUE)
ggsave("g:/gBGC/carl/log/52_B_CF.png")

### Regarding CF, we need to see the data without binning
data %>% ungroup %>% 
    filter(unitig == 0) %>% 
    select(genospecies, unitig, gene, post_mean, GC3) %>% 
    mutate(post_mean = log(post_mean)) %>% 
    ggplot(aes(GC3, post_mean)) + 
    geom_point() + 
    facet_wrap(~genospecies) + 
    geom_smooth()
ggsave("g:/gBGC/carl/log/53_C_CF.png")




### Let's isolate the most extreme values of recombination, and look at GC3
A = data %>% ungroup %>% 
    filter(unitig == 0) %>% 
    group_by(genospecies) %>% 
    filter(post_mean > quantile(post_mean, .95)) %>% 
    ggplot(aes(mid, post_mean, color = GC3)) +
    geom_point()

B = data %>% ungroup %>% 
    filter(unitig == 0) %>% 
    group_by(genospecies) %>% 
    filter(post_mean > quantile(post_mean, .95)) %>% 
    ggplot(aes(mid, GC3)) +
    geom_point()

ggarrange(A, B, ncol = 1)


top_percent = 1
data %>% ungroup %>% 
    filter(unitig == 0) %>% 
    group_by(genospecies) %>% 
    mutate(`recombination class` = if_else(post_mean > quantile(post_mean, 1-(top_percent/100)), paste0("top ", top_percent,"%"), "rest")) %>% 
    
    ggplot(aes(mid, GC3, color = `recombination class`)) + 
    geom_point() +
    facet_wrap(~genospecies)+ 
    labs(x = 'position')
ggsave("g:/gBGC/carl/log/54_position_top10.png")


top_percent = 1
data %>% ungroup %>% 
    filter(unitig == 0) %>% 
    group_by(genospecies) %>% 
    mutate(`recombination class` = if_else(post_mean > quantile(post_mean, 1-(top_percent/100)), paste0("top ", top_percent,"%"), "rest")) %>% 
    
    ggplot(aes(GC3, `recombination class`, fill = `recombination class`)) + 
    geom_boxplot() +
    #facet_grid(.~genospecies)
    facet_grid(genospecies~.)

ggsave("g:/gBGC/carl/log/55_boxplot_top10.png")



### Create data that is used for the shiny app
# First we have to create the sliding window.

#import gene map
gene_map <- read_delim("C:/Users/carl/urecomb/lokal/exonerate/gene_map.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
    rename(genospecies = sample_genospecies,
           unitig = sample_unitig,
           gene = sample_gene) %>% rowid_to_column() %>% group_by(rowid) %>%  mutate(SM3_mid = mean(c(SM3_start, SM3_end)))

data = inner_join(data, gene_map)

# if the unitigs agree, we can remove one of them
data %>% transmute(vs = paste(unitig, SM3_unitig)) %>% table
#   0 0   1 1   2 2   3 3 
# 12535  1301  1176   511
data = data %>% select(-SM3_unitig) 



# How have the genes jumped with the SM3 reference (sanity check)
data %>% mutate(diff = mid - SM3_mid) %>%
    ggplot(aes(diff)) +
    geom_histogram(bins = 100) +
    facet_grid(unitig ~ genospecies) + 
    labs(caption = "# bins = 100") +
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
ggsave("Enquire_Maria.png", height = 8, width = 10)


#test:
data_anim = data %>% 
    select(genospecies, unitig, gene, GC3, post_mean, post_var, plasmid, SM3_mid, length) %>% 
    rename(mid = SM3_mid)
sel_gs = 'C'

#for (roll_width in c(100, 250, 500, 750, 1000)) {
widths = c(1, 100, 1000)
widths = c(1, 5, seq(10, 500, 5))
data_anim_binds = tibble(); for (roll_width in widths) {
    print(roll_width)
    data_anim_binds = bind_rows(data_anim %>%
                                    arrange(mid) %>% 
                                    group_by(genospecies, plasmid) %>% 
                                    
                                    mutate(roll_post_mean = rollapply(post_mean, roll_width, median, fill = NA)) %>%  #View
                                    mutate(roll_GC3 = rollapply(GC3, roll_width, median, fill = NA)) %>%
                                    
                                    mutate(roll_width = roll_width) %>% 
                                    drop_na(),
                                data_anim_binds)
}
saveRDS(data_anim_binds, 'C:/Users/carl/Documents/test2/data_anim_binds_6_median_GC.rds')
                                


render_frame <- function(sel_gs, sel_roll_width) {
    data_anim_binds %>%
        filter(genospecies == sel_gs) %>% 

        filter(roll_width == sel_roll_width) %>%
        ggplot(aes(mid, roll_post_mean, color = roll_GC3)) + 
            geom_point() + 
            facet_grid(plasmid~., scales = "free") + 
            scale_color_gradientn(colours = c('red1', 'grey', 'green4')) + 
        labs(subtitle = paste0("Genospecies ", sel_gs, "\nrolling window width: ", sel_roll_width, " genes"))
}

# test it 
render_frame('C', 100)






data_anim_binds %>% filter(genospecies == sel_gs) %>% 
    filter(roll_width == 100) %>%
    ggplot(aes(mid, roll_post_mean, color = roll_GC3)) + 
        geom_line() + 
        facet_grid(plasmid~., scales = "free") + 
        scale_color_gradientn(colours = c('red1', 'grey', 'green4'))

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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#################

phi_data %>% ggplot(aes(p_phi_permut)) + geom_histogram()
saveRDS(phi_data, "PHI_final.rds")
write_tsv(phi_data, "PHI_final.tsv")


# Check distribution of p-values from PHI
phi_data %>% pivot_longer(starts_with("p_")) %>% 
    ggplot(aes(value)) + 
    facet_grid(name~genospecies) + 
    geom_histogram()

