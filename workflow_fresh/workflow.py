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


GC3_output_file = f"{base_path}/output/GC3.tab"

gwf.target(sanify(F"GC"),
	inputs = [], 
	outputs = [GC3_output_file],
	cores = 1,
	walltime = '2:00:00',
	memory = '4gb',
	account = 'gBGC') << f"""

		cd data/tabseq/genospecies_unitig_genes/fasta

		# Reset output.
		
		echo -e "header\tGC3\tgenospecies\tunitig\tfile" > {GC3_output_file}
		touch GC_done && rm GC_done

		# Loop through all files and output GC3
		for i in {{A..E}}; do
			for j in {{0..3}}; do
				for k in ${{i}}/${{j}}/*.fasta; do
					echo $k;
					l = $(basename $k)
					cat $k | {base_path}/scripts/fasta_gc.py $i $j $l >> {GC3_output_file}
				done;
			done;
		done

		touch GC_done
	"""

for sel_gs in ['A', 'B', 'C', 'D', 'E']:
	for sel_un in ['0', '1', '2', '3']:
		path = f"data/tabseq/genospecies_unitig_genes/fasta/{sel_gs}/{sel_un}"

		tree_path = tree_base + f"/{sel_gs}/{sel_gs}.fasta.raxml.bestTree"

		for _gene, gene in enumerate(glob.glob(path + '/*.tabseq.fasta', recursive = False)):

			gene_file = gene.split('/')[-1]
			gene_stem = '.'.join(gene_file.split('.')[:-2]) # I opted to also remove the trivial tabseq part.


			#debug
			#print(gene_file, gene_stem)




			
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
					cd {gene_stem}_cf && rm * || echo empty


					# Execute CF.
					ClonalFrameML {base_path}/{tree_path} {base_path}/{gene} {gene_stem}


					# Collect results?

				"""



CF_collect_file = f"{base_path}/output/CF_collected.tab"

gwf.target(sanify(F"ph2_collect_CF"),
	inputs = [], 
	outputs = [CF_collect_file],
	cores = 1,
	walltime = '2:00:00',
	memory = '4gb',
	account = 'gBGC') << f"""

		
		cd data/tabseq/genospecies_unitig_genes/fasta

		# Reset output.
		
		echo "genospecies	unitig	parameter	post_mean	post_var	a_post	b_post" > {GC3_output_file}
		touch CF_done && rm CF_done

		# Loop through all files and output GC3
		for i in {{A..E}}; do
			for j in {{0..3}}; do
				for k in ${{i}}/${{j}}/*cf/*.em.txt; do
					echo $k;
					l=$(basename $k)
					grep -E "(^nu)|(^1/delta)|(^R/theta)" $k | awk -v vari="$i" -v varj="$j" -v vark="$l" '{{ print vari"\t"varj"\t"vark"\t" $0}}' >> {CF_collect_file}
					
				done;
			
			done;
			
		done

		touch CF_done


	"""
		