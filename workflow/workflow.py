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
    "mail_user": "kobel@pm.me",
    "mail_type": "FAIL",
})




""" 
if not path.isdir(f'output/{title}'):
    try:
        subprocess.run(f'mkdir -p output/{title}; AssemblyAlignmentGenerator/AssemblyAlignmentGenerate_fe AssemblyAlignmentGenerator/assembly_summary_refseq.txt accession_list.txt output/{title} {title}', shell = True, check = True)
    except subprocess.CalledProcessError as e:
        print(f'\nAn error occured while initializing:\n', e)
        sys.exit()
 """


title = 'Rlegum'
bin_size = 5000

genomes = glob.glob('genomes/corrected_header/*')

for genome in genomes:
    genome_basename = os.path.basename(genome)
    genome_stem = os.path.splitext(genome_basename)[0]

    print(genome)
    print(genome_basename)
    print(genome_stem)
    print()


    file_binned_xmfa_out = f"{genome_stem}_binned_{int(bin_size)}bp"
    
    gwf.target(sanify('aag_bin_' + title),
        inputs = f"{genome}",
        outputs = [f"output/{title}/{file_binned_xmfa_out}.xmfa", f"output/{title}/{file_binned_xmfa_out}_gc.tab", f"output/{title}/metadata.tab"],
        cores = 1,
        walltime = '1:00:00',
        memory = '4gb', account = "clinicalmicrobio") << f"""

    
    mkdir -p output/{title}
    cd output/{title}
    
    # bin
    ../../script/xmfa_bin.py ../../{genome} {int(bin_size)} > {file_binned_xmfa_out}.xmfa


    # compute GC
    ../../script/xmfa_gc.py {file_binned_xmfa_out}.xmfa > {file_binned_xmfa_out}_gc.tab


    # write pipeline metadata
    echo -e "{title}\t{bin_size}\t{genome_stem}\t{genome}\t{file_binned_xmfa_out}.xmfa\t{file_binned_xmfa_out}_gc.tab" >> metadata.tab

    """









    