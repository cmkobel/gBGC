library(tidyverse)
soiltypes =  read_delim("~/genomedk/gBGC/carl/Rhizobium_project_MICA/Rhizobium_soiltypes_new.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
custom = function(data) {
    paste(list(table(data)))
}

summary = soiltypes %>% 
    group_by(Genospecies) %>% 
    mutate(n = length(`Seq ID`)) %>% 
    group_by(n, add = T) %>% 
    summarize_at(vars(Origin, Origin1, Origin2), custom)

summary %>% unnest(Origin)    


write.csv2(summary, file = 'genospecies_structure.tsv')# delim = '\t')

soiltypes$`Seq ID` %>% sort
                 

#3240