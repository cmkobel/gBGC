import gwf
import sys
import glob

def sanify(input):
    """ Makes sure that the name of the gwf target is not illegal. """
    output = []
    
    for i in str(input):
        
        ascii = ord(i)
        if (ascii >= 48 and ascii <= 57) or (ascii >= 65 and ascii <= 90) or (ascii >= 97 and ascii <= 122) or ascii == 95:
            output.append(i)
        else:
            output.append('_')

    return ''.join(output)


gwf = gwf.Workflow(defaults={
    #"mail_user": "kobel@pm.me", "mail_type": "FAIL",
})

base_path = f"/home/cmkobel/gBGC/carl/workflow_fresh"
tree_base = f"data/tabseq/genospecies_core_genomes"

for sel_gs in ['A', 'B', 'C', 'D', 'E']:
	for sel_un in ['0', '1', '2', '3']:
		path = f"data/tabseq/genospecies_unitig_genes/fasta/{sel_gs}/{sel_un}"

		tree_path = tree_base + f"/{sel_gs}/{sel_gs}.fasta.raxml.bestTree"

		for gene in glob.glob(path + '/*.fasta'):

			gene_file = gene.split('/')[-1]
			gene_stem = '.'.join(gene_file.split('.')[:-2]) # I opted to also remove the trivial tabseq part.


			print(gene_file, gene_stem)
			
			gwf.target(sanify(f"CF_{sel_gs}_{sel_un}_{gene_stem}"),
				inputs = [gene],
				outputs = [],
				cores = 1,
				walltime = '2:00:00',
				memory = '4gb',
				account = 'gBGC') << f"""

				# Debug.
				echo {gene_stem}
				ls -l {gene}
				ls -l {tree_path}

				# Make dir.
				cd {path}
				mkdir -p {gene_stem}_cf
				rm -r {gene_stem}_cf
				mkdir {gene_stem}_cf

				# Execute CF.
				ClonalFrameML {base_path}/{tree_path} {base_path}/{gene} {gene_stem}


				# Collect results?
				"""
			break
		break
	break

