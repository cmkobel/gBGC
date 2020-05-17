This pipeline uses a syntenic cg-alignment of Rleguminosarum.

It concatenates the syntenic core genes in bins of a given size. This size is given in the `bin-size` variable in `workflow.py`. Then it infers recombination with the Kussell-lab mcorr-package and also with the Bruen PHI-package.

It extracts the core genes in the alignment.

From each gene it measures GC content (GC, GC1, GC2, GC3).
It also measures recombination. This is done with PHI and ClonalFrameML.

ClonalFrameML needs a phylogenetic tree, this is inferred with raxml-ng for each genospecies.
The phylogenetic trees will have to be manually created in advance.

Because the recombination inference jobs are dependent on the name of each gene, the pipeline will have to be run in multiple steps.

First step extracts genes and measures GC-content.
Second step measures recombination with PHI and ClonalFrameML.
Third (last) step collects the parameters from all genes into a single file for each recombination inference tool.
