This pipeline isolates the genes from a .xmfa alignment. Subsequently inferring the recombination rate with Bruen's implementation of PHI (Pairwise Homoplasy Index).

If you wish to infer recombination with ClonalFrameML instead, use the workflow in [`../workflow_fresh`](https://github.com/cmkobel/gBGC/tree/master/workflow_fresh).


This pipeline also measures the GC3-content.

Because the recombination inference jobs are dependent on the name of each gene, the pipeline will have to be run in multiple steps.

First step extracts genes and measures GC-content.
Second step measures recombination with PHI.
Third (last) step collects the results for each gene into one tab-file.

Dependencies:
 * python3
 * gwf-org
 * slurm workload manager
