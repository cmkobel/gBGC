#Author Carl



from gwf import *
import subprocess, sys
import glob
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

""" 
if not path.isdir(f'output/{title}'):
    try:
        subprocess.run(f'mkdir -p output/{title}; AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_fe AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}', shell = True, check = True)
    except subprocess.CalledProcessError as e:
        print(f'\nAn error occured while initializing:\n', e)
        sys.exit()
 """




genomes = glob.glob('genomes/corrected_header*')

for genome in genomes:
    print(genome)




    gwf.target(sanify('aag_bin_' + title),
                    inputs = [],#'accession_list.txt'],
                    outputs = [f'output/{title}/{title}_core.xmfa', f'output/{title}/accession_list.txt'],
                    cores = 8,
                    walltime = '24:00:00',
                    memory = '16gb', account = "clinicalmicrobio") << f"""

    mkdir -p output/{title}

    cp accession_list.txt output/{title}

    AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_be AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}

    """









    break