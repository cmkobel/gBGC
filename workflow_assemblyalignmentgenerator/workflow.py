#Author Carl



from gwf import *
import subprocess, sys

import os.path
from os import path
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


title = "Spyo1"


if not path.isdir(f'output/{title}'):
    try:
        subprocess.run(f'mkdir -p output/{title}; AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_fe AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}', shell = True, check = True)
    except subprocess.CalledProcessError as e:
        print(f'\nAn error occured while initializing:\n', e)
        sys.exit()


gwf.target(sanify('aag_xmfa_' + title),
                  inputs = ['accession_list.txt'],
                  outputs = [f'output/{title}/{title}_core.xmfa', f'output/{title}/accession_list.txt'],
                  cores = 8,
                  walltime = '24:00:00',
                  memory = '16gb', account = "clinicalmicrobio") << f"""

mkdir -p output/{title}

cp accession_list.txt output/{title}

AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_be AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}

"""


gwf.target(sanify('aag_mcorr_xmfa_' + title),
                  inputs = [f'output/{title}/{title}_core.xmfa'],
                  outputs = [f'output/{title}/mcorr-xmfa.completed'],
                  cores = 8,
                  memory = '16gb',
                  account = "clinicalmicrobio") << f"""

cd output/{title}


# Isolate each gene, to infer recombination independently.

mkdir -p single_genes
cd single_genes
../../../scripts/xmfa_split.py ../{title}_core.xmfa



# mcorr-xmfa on each gene
ls *.fa | parallel -j 8 'mcorr-xmfa {{}} {{.}}_mx --num-boot=5'

touch ../mcorr-xmfa.completed
"""


gwf.target(sanify('aag_mcorr_fit_' + title),
                  inputs = [f'output/{title}/mcorr-xmfa.completed'],
                  outputs = [f'output/{title}/mcorr-fit.completed',
                             f'output/{title}/{title}_recomb.csv'],
                  cores = 16,
                  memory = '32gb',
                  walltime = '1:00:00', # 24
                  account = 'clinicalmicrobio') << f"""

cd output/{title}

cd single_genes

echo here0

# Run mcorr-fit
#ls *_mx.csv | parallel -j 16 'mcorr-fit {{}} {{.}}_mf'


echo here1

# extract from txt


echo here2

# Write all fit_results.csv's to a single file
echo "#" > ../all_genes_recomb_data.csv

echo here3



for f in *mx_mf_fit_results.csv; do printf $f, >> ../all_genes_recomb_data.csv; awk '/^all,/' $f >> ../all_genes_recomb_data.csv; echo "" >> ../all_genes_recomb_data.csv; done

echo here4

awk '/.*all,/' ../all_genes_recomb_data.csv > ../{title}_recomb.csv

echo here5

touch ../mcorr-fit.completed

"""



gwf.target(sanify('aag_gc_' + title),
           inputs = [f'output/{title}/{title}_core.xmfa'],
           outputs = [f'output/{title}/{title}_gc.tab'],
           cores = 1,
           memory = '1gb') << f"""

cd output/{title}

../../scripts/xmfa_gc.py {title}_core.xmfa > {title}_gc.tab
           
"""

