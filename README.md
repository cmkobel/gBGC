# gBGC
Investigating GC-Biased Gene Conversion (gBGC) in Rhizobium leguminosarum

I'm developing a number of pipelines that I use to utilize different software in order to measure recombination in bacteria

1. workflow_reads maps reads to a manually picked reference. Recombination is then inferred on the cg-alignment or each gene in the reference.
1. workflow_assemblyalignmentgenerator uses the pipeline developed by kussell-lab.
1. workflow_mcorr uses syntenic cg-alignments. Specifically for Rleguminosarum
1. workflow_PHI uses Bruen's recombination inference tool.




# Log
I want to investigate the relationship between recombination-rate and GC3 content. I found a new method for inferring recombination rate in bacteria [kussell-lab/mcorr](https://github.com/kussell-lab/mcorr) and wanted to use it to indirectly replicate [Lassalle et al. 2015](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004941). In this study they use PHI which I will come back around later. The GC3 content is easy to measure. Just align the reading frame, and count the proportion of G's and C's.

I wanted to have a strong signal, so I chose the species with the strongest relation. in Lassalle et al. 2015 they show that a linear relationship between recombination and GC3 has an R^2 = .68 for Streptococcus pyogenes. Thus I randomly selected ~60 S. pyogenes genomes from genbank, aligned the core genomes and measured the recombination for each gene. When I use kussell-lab/mcorr, I get the following results:

![](https://raw.githubusercontent.com/cmkobel/gBGC/master/log/1_spyo1.png)

This figure shows that either there is no relationship, or something is wrong with the method. The problem is that mcorr most often infers that recombination is lacking (close to zero). 

I think the problem is that even though there is recombination, it is hard to infer it from single genes. Recombination is measured using linkage-disequilibrium (LD). If the sequences are too short (average gene length 1000) means that there is not much LD to measure. If we concatenate syntenic genes and inferm recombination from these, we might have enough signal to actually measure recombination.

I found a dataset of core genes from 5 genospecies of Rhizobium leguminosarum. Other studies on this data suggests that the recombination/GC-bias should be existing in this dataset. I concatenated the genes in bins of different sizes, inferred recombination. The code for this pipeline can be found in [workflow/](https://github.com/cmkobel/gBGC/tree/master/workflow). These are the results I got.

![](https://raw.githubusercontent.com/cmkobel/gBGC/master/log/2_Rleg_phi_pool_all.png)

As we can see, again, there is no relationship between the GC3-content and recombination rate. Maybe there is a relationship for genospecies D, bin size 20000. But only if you remove all the points in the bottom.

I tried removing the tails of the distribution, but this doesn't help on the relationship. I'm starting to think that either something is wrong with this method, or I'm using it wrong. As a sanity check, in order the check if there is a relationship in this data or not, using a method that we know works (PHI in Lassalle et al. 2015), we can check that the data we're feeding into mcorr actually has a signal.

I ran PHI for each gene in the core genome of these 5 genospecies of R. leguminosarum. The pipeline for this can be found at [workflow_PHI/](https://github.com/cmkobel/gBGC/tree/master/workflow_PHI), and these are the results I got:

![](https://github.com/cmkobel/gBGC/blob/master/log/3_PHI_ratio.png)

_**Figure 3**: 20 bins with an equal number of genes in each. mean_GC3 is the mean GC3 content in that bin and ratio is the number of significantly recombining genes over the total number of genes in the bin. Significance of recombination is inferred with PHI._



