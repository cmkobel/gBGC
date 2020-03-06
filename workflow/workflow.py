#Author Carl



from gwf import *
import subprocess, sys
import glob
import os.path

#from workflow_templates import *


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


gwf = Workflow(defaults={
    #"mail_user": "kobel@pm.me", "mail_type": "FAIL",
})




""" 
if not path.isdir(f'output/{title}'):
    try:
        subprocess.run(f'mkdir -p output/{title}; AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_fe AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}', shell = True, check = True)
    except subprocess.CalledProcessError as e:
        print(f'\nAn error occured while initializing:\n', e)
        sys.exit()
 """


title_prefix = 'Rlegum' 
genome_dir = '../Rhizobium_project_MICA/corrected_headers'

#bin_sizes = [1, 5000, 8000, 10000, 15000, 20000, 25000, 30000, 40000] 
bin_sizes = [20000]
#bin_sizes = [50000, 70000, 80000, 100000, 120000] #3


for bin_size in bin_sizes:
    #title = title_prefix + str(bin_size)
    #genomes = glob.glob('genomes/corrected_header_2/*')
    genomes = glob.glob(genome_dir + '/*.xmfa')

    for genome in genomes:
        genome_basename = os.path.basename(genome)
        genome_stem = os.path.splitext(genome_basename)[0]

        title = title_prefix + '_' + str(genome_stem) + '_' + str(bin_size) + '_2'
        #title = f"{title_prefix} + '_' + {str(genome_stem)} + '_' + {str(bin_size)}"
        

        print(f"{genome} ({genome_stem})")

        #print(genome_basename)
        #print(genome_stem)
        


        file_binned_xmfa_out = f"{genome_stem}_binned_{int(bin_size)}bp"
        
        gwf.target(sanify('A_bin_' + title),
            inputs = f"{genome}",
            outputs = [f"output/{title}/{file_binned_xmfa_out}.xmfa", f"output/{title}/{file_binned_xmfa_out}_gc.tab"],
            cores = 1,
            walltime = '1:00:00',
            memory = '4gb',
            account = "gBGC") << f"""

        
        mkdir -p output/{title}
        cd output/{title}
        
        # bin
        ../../script/xmfa_bin.py ../../{genome} {int(bin_size)} > {file_binned_xmfa_out}.xmfa


        # compute GC from the binned genome
        ../../script/xmfa_gc.py {file_binned_xmfa_out}.xmfa {bin_size} {genome_stem} > {file_binned_xmfa_out}_gc.tab




        """
        


        gwf.target(sanify('A_split_' + title),
            inputs = [f"output/{title}/{file_binned_xmfa_out}.xmfa"],
            outputs = [f"output/{title}/leg_split.completed"],
            cores = 1,
            walltime = '1:00:00',
            memory = '1gb',
            account = "gBGC") << f"""

        cd output/{title}
        
        mkdir -p split
        cd split


        # split the binned xmfa into single xmfa files for mcorr
        ../../../script/xmfa_split.py ../{file_binned_xmfa_out}.xmfa
        touch ./xmfasplitted

        # split the binned xmfa into single fa files for phi
        ../../../script/xmfa_split_to_fa.py ../../../{file_binned_xmfa_out}.xmfa
        touch ./fasplitted

        touch ../leg_split.completed

        """
        break # debug for the first genome only
        '''

        # When the above targets are done, call it again, to run the rest of the targets below:
        splitted_xmfas = glob.glob(f"output/{title}/split/*.xmfa")
        
        for num, single_gene in enumerate(splitted_xmfas):
            single_gene_basename = os.path.basename(single_gene)
            single_gene_stem = os.path.splitext(single_gene_basename)[0]
            #print(single_gene_basename)


            gwf.target(sanify('B_mcorr_' + title + '_' + str(num) + '_' + single_gene_stem),
                inputs = [single_gene],
                outputs = [f"output/{title}/split/{single_gene_stem}.csv",
                        f"output/{title}/split/{single_gene_stem}_fit_results.csv",
                        f"output/{title}/split/{single_gene_stem}_fitpar.csv"],
                cores = 1,
                walltime = '6:00:00',
                memory = '12gb',
                account = "gBGC") << f"""

    cd output/{title}/split

    mcorr-xmfa {single_gene_basename} {single_gene_stem}

    mcorr-fit {single_gene_stem}.csv {single_gene_stem}

    printf {single_gene_stem}, >> {single_gene_stem}_fitpar.csv
    awk '/^all,/' {single_gene_stem}_fit_results.csv >> {single_gene_stem}_fitpar.csv
    #echo "" >> {single_gene_stem}_fitpar.csv; done


        
        """
            

        # Merge the _fitpar.csv files together, so it can be imported in R later.
        fitpars = glob.glob(f"output/{title}/split/*_fitpar.csv")
        #for num, fitpar in enumerate(fitpars): # I have no idea why this job existed for each 
        
        gwf.target(sanify('C_merge_recomb_' + title),
            inputs = [i for i in fitpars],
            outputs = [f"output/{title}/{file_binned_xmfa_out}_fitpars.csv"],
            cores = 1,
            walltime = '10:00',
            memory = '1gb',
            account = "gBGC") << f"""


cd output/{title}

inputfile="split/*_fitpar.csv"

outputfile="{file_binned_xmfa_out}_fitpars.csv"




if stat -t $inputfile >/dev/null 2>&1; then
    for fitpar in $inputfile; do
        awk '{{print $0, ",{bin_size}, {genome_stem}"}}' $fitpar >> $outputfile
    done
fi


    


            """
        

'''