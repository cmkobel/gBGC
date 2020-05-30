This pipeline infers the recombination rate (using ClonalFrameML) on each gene in a set of genes specified in the workflow.

If you wish to infer recombination with PHI instead, use the workflow in [`../workflow`](https://github.com/cmkobel/gBGC/tree/master/workflow).


Each gene must be isolated in advance an put in the specified directory. The phylogenetic tree must also be supplied manually.
There is no reason these things are not automated, other than that the lack of time permits only a limited degree of automation.


![](https://raw.githubusercontent.com/cmkobel/gBGC/master/log/workflow.png)

Dependencies:
 * python3
 * gwf-org
 * slurm workload manager

