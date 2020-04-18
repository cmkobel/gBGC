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

			GC3_output_file = f"{base_path}/output/GC3.tab"

			gwf.target(sanify(F"GC_{sel_gs}_{sel_un}_{gene_stem}"),
				inputs = [gene], 
				outputs = [GC3_output_file],
				cores = 1,
				walltime = '2:00:00',
				memory = '4gb',
				account = 'gBGC') << f"""

				echo {gene_stem}
				ls -l {gene}
				ls -l {tree_path}

				cd data/tabseq/genospecies_unitig_genes/fasta

				# Reset output.
				
				echo -e header\tGC3\tgene\tgs\tunitig\tfile > {GC3_output_file}

				# Loop through all files and output GC3
				for i in {{A..E}}; do
					for j in {{0..3}}; do
						for k in ${{i}}/${{j}}/*.fasta; do
							echo $k;
							cat $k | {base_path}/scripts/fasta_gc.py {gene_stem} {sel_gs} {sel_un} $k >> {GC3_output_file}
						done;
					done;
				done

				"""



			
			gwf.target(sanify(f"CF_{sel_gs}_{sel_un}_{gene_stem}"),
				inputs = [gene],
				outputs = [f"{path}/{gene_stem}_cf/{gene_stem}.em.txt"],
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
				cd {gene_stem}_cf


				# Execute CF.
				ClonalFrameML {base_path}/{tree_path} {base_path}/{gene} {gene_stem}


				# Collect results?
				"""
			break
		break
	break

