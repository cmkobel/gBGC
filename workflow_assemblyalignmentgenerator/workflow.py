#Author Carl



from gwf import *
import subprocess, sys

import os.path
from os import path
import glob
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
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
})


title = "Spyo3"


if not path.isdir(f'output/{title}'):
    try:
        subprocess.run(f'mkdir -p output/{title}; AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_fe AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}', shell = True, check = True)
    except subprocess.CalledProcessError as e:
        print(f'\nAn error occured while initializing:\n', e)
        sys.exit()


gwf.target(sanify('aag_xmfa_' + title),
                  inputs = ['accession_list.txt'],
                  outputs = [f'output/{title}/{title}_core.xmfa',
                             f'output/{title}/accession_list.txt'],
                  cores = 8,
                  walltime = '24:00:00',
                  memory = '16gb', account = "clinicalmicrobio") << f"""

mkdir -p output/{title}

cp accession_list.txt output/{title}

#Assemble genomes
AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_be AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}


"""

gwf.target(sanify('aag_split_' + title),
    inputs = [],
    outputs = [f"output/{title}/aag_split.completed"],
    cores = 1,
    walltime = '1:00:00',
    memory = '16gb',
    account = 'clinicalmicrobio') << f"""

# Split into individual genes
mkdir -p output/{title}/single_genes
cd output/{title}/single_genes
../../../scripts/xmfa_split.py ../{title}_core.xmfa

touch ../aag_split.completed
"""

for single_gene in glob.glob(f"output/{title}/single_genes/*.fa"):
    single_gene_basename = os.path.basename(single_gene)
    single_gene_stem = os.path.splitext(single_gene_basename)[0]
    #print(single_gene, single_gene_basename, single_gene_stem)

    gwf.target(sanify('aag_mcorr_' + title + '_' + single_gene_stem),
                    inputs = [f'output/{title}/{title}_core.xmfa'],
                    outputs = [f"output/{title}/single_genes/{single_gene_stem}_mx.csv", f"output/{title}/single_genes/{single_gene_stem}_mx.json",
                               f"output/{title}/single_genes/{single_gene_stem}_mf_fit_results.csv"],
                    cores = 1,
                    memory = '16gb',
                    account = "clinicalmicrobio") << f"""

    cd output/{title}/single_genes


    


    #mcorr-xmfa {single_gene_basename} {single_gene_stem}_mx
    echo mcorr-xmfa done

    mcorr-fit {single_gene_stem}_mx.csv {single_gene_stem}_mf
    echo mcorr_fit done

    #printf "{single_gene_stem}" > {single_gene_stem}_final.csv
    #cat {single_gene_stem}_mf_fit_results.csv | grep "^all" >> {single_gene_stem}_final.csv
    
    """



gwf.target(sanify('aag_gc_' + title),
           inputs = [f'output/{title}/{title}_core.xmfa'],
           outputs = [f'output/{title}/{title}_gc.tab'],
           cores = 1,
           memory = '1gb',
           account = 'clinicalmicrobio') << f"""

cd output/{title}

../../scripts/xmfa_gc.py {title}_core.xmfa > {title}_gc.tab
           
"""

