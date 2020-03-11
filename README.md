# gBGC
Investigating GC-Biased Gene Conversion (gBGC) in Rhizobium leguminosarum

I'm developing a number of pipelines that I use to utilize different software in order to measure recombination in bacteria

1. workflow_reads maps reads to a manually picked reference. Recombination is then inferred on the cg-alignment or each gene in the reference.
1. workflow_assemblyalignmentgenerator uses the pipeline developed by kussell-lab.
1. workflow_mcorr uses syntenic cg-alignments. Specifically for Rleguminosarum
1. workflow_PHI uses Bruen's recombination inference tool.




# Log
I want to investigate the relationship between recombination-rate and GC3 content. I found a new method for inferring recombination rate in bacteria [kussell-lab/mcorr](https://github.com/kussell-lab/mcorr) and wanted to use it to indirectly replicate [Lassalle et al. 2015](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004941). In this study they use PHI which I will come back around later. The GC3 content is easy to measure. Just align the reading frame, and count the proportion of G's and C's.

I wanted to have a strong signal, so I chose the species with the strongest relation. in Lassalle et al. 2015 they show that a linear relationship between recombination and GC3 has an R^2 = .68 for Streptococcus pyogenes. Thus I randomly selected ~60 genomes from genbank, aligned the core genomes and measured the recombination for each gene. When I use kussell-lab/mcorr, I get the following results:

![alt text](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png)




First thing I tried running mcorr-xmfa and mcorr-fit on some randomly  a set of syntenic coregenes in 5 genospecies of Rhizobium leguminosarum. The code for this pipeline can be found in [workflow/](https://github.com/cmkobel/gBGC/tree/master/workflow). 
