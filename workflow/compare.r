# This script compares the recombination values inferred from mcorr and phi. Possibly later also ClonalFrame.


library(tidyverse)
setwd("~/urecomb/lokal/exports/export4compare/")

## GC
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



## mcorr    import fitted (recombination) parameters
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


## PHI files
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


