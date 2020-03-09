This pipeline uses a syntenic cg-alignment of Rleguminosarum.

It concatenates the syntenic core genes in bins of a given size. This size is given in the `bin-size` variable in `workflow.py`. Then it infers recombination with the Kussell-lab mcorr-package and also with the Bruen PHI-package.

Because there will be a job for each gene, it is necessary to run the pipeline several times. First stage bins the genome in the wanted sizes, and splits the genome into concatenated gene-files. The second stage runs mcorr and PHI on each of those concatenated gene-files. The third stage writes the output of each method into a single file, which can then be used for further handling.
