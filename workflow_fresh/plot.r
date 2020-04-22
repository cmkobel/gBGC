library(tidyverse)
library(zoo)
library(gganimate)


CF_raw = read_delim("G:/gBGC/carl/workflow_fresh/output/CF_collected.tab", 
                           "\t", escape_double = FALSE, trim_ws = TRUE,
                 col_names = c('genospecies', 'unitig', 'file', 'parameter', 'post_mean', 'post_var', 'a_post', 'b_post')) %>% 
    mutate(gene = str_sub(file, 4, str_length(file)-7)) %>% select(-file)


#CF = CF_raw %>% pivot_wider(names_from = parameter, values_from = c(post_mean, post_var, a_post, b_post))
CF = CF_raw %>% filter(parameter == "R/theta") %>% select(-parameter)


GC3_raw = read_delim("G:/gBGC/carl/workflow_fresh/output/GC3.tab", 
                  "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    rename(genospecies = gs) %>%
    mutate(sample = str_sub(header, 6),
           gene = str_sub(file, 8, str_length(file)-13)) %>% 
    select(-header, -file)

GC3 = GC3_raw %>% group_by(genospecies, unitig, gene) %>% 
    summarize(mean_GC3 = mean(GC3), n_samples_GC3 = length(GC3))
    



gff_data <- read_csv("C:/Users/carl/Desktop/xmfa2mfa/3206-3.gff") %>% select(-X1) %>% mutate(mid = start + ((end-start)/2), length = end-start)
gff2_data <- read_delim("C:/Users/carl/Desktop/xmfa2mfa/Gene_function_pop_gene_data.csv", delim = ';')
gff_data = left_join(gff_data, gff2_data %>% select(`Gene group`, `Putative function`) %>% rename(gene_group = `Gene group`)) %>% rename(gene = gene_group)
rm(gff2_data)


data = inner_join(GC3, CF) %>% inner_join(gff_data) %>% inner_join(phi_data %>% rename(gene_group = gene), mutate(unitig = as.character(unitig)))
write_tsv(data, "CF_Rleg_coregenometrees.tab")


data %>% filter(genospecies == 'C') %>%  ggplot(aes(mid, post_mean, color = mean_GC3^4)) + 
    geom_point() + 
    facet_grid(plasmid~., scales = "free") + 
    scale_color_gradientn(colours = c('red1', 'grey', 'green4'))
ggsave("40_cf_coregenetrees.png")



# Create animation

# First we have to create the sliding window.

#test:
data_anim = data %>% 
    select(genospecies, unitig, gene, GC3 = mean_GC3, post_mean, post_var, plasmid, mid, length)
sel_gs = 'C'

#for (roll_width in c(100, 250, 500, 750, 1000)) {
widths = c(1, seq(10, 1000, 10))
widths = c(1, 100, 1000)
data_anim_binds = tibble(); for (roll_width in widths) {
    print(roll_width)
    data_anim_binds = bind_rows(data_anim %>%
                                    arrange(mid) %>% 
                                    group_by(genospecies, plasmid) %>% 
                                    
                                    mutate(roll_post_mean = rollapply(post_mean, roll_width, median, fill = NA)) %>%  #View
                                    mutate(roll_GC3 = rollapply(GC3, roll_width, mean, fill = NA)) %>%
                                    
                                    mutate(roll_width = roll_width) %>% 
                                    drop_na(),
                                data_anim_binds)
}
saveRDS(data_anim_binds, 'C:/Users/carl/Documents/test2/data_anim_binds_2.rds')
                                


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
phi_data = gather_PHI() %>% select(genospecies, unitig, infsites, gene, p_phi_permut) %>%  mutate(gene = paste0("group", gene)) 



phi_data %>% ggplot(aes(p_phi_permut)) + geom_histogram()
saveRDS(phi_data, "PHI_final.rds")
write_tsv(phi_data, "PHI_final.tsv")


# Check distribution of p-values from PHI
phi_data %>% pivot_longer(starts_with("p_")) %>% 
    ggplot(aes(value)) + 
    facet_grid(name~genospecies) + 
    geom_histogram()

