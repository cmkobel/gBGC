#Author Carl


import pandas as pd
from gwf import *
from workflow_templates import *


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

target_dir = '/faststorage/project/ClinicalMicrobio/10carl/urecomb/workflow/output'



### I: Reference genomes and their indexing ###

# group in I and group in II must be equal.

reference_genomes = {'S_mitis_tg': 'mitis_b6'}
                      #'Aactinomycetemcomitans': 'path'}
for group, reference_stem in reference_genomes.items(): # index all genomes
    print(group, reference_stem)

    gwf.target_from_template(sanify('ure_smaltindex_' + group),
                             smalt_index(target_dir, group, reference_stem))




### II: For each sample in the data set ###
input_data = pd.read_csv('input.tab', delimiter = '\t')
print(' ', input_data.columns)

for sample_, sample in input_data.iterrows():

    # Todo: brug et andet variabelnavn end sample. clasher med sample_name

    group = sample['group'].strip()
    sample_name = sample['sample'].strip()
    reference_stem = reference_genomes[group].strip()

    print(' ', group, sample_name, reference_stem)

    # Map reads for each isolate to the reference
    gwf.target_from_template(sanify('ure_smaltmap_' + group[0:5] + '_' + sample_name),
                             smalt_map(group, sample_name, sample['forward'].strip(), sample['reverse'].strip(), reference_stem))

    # infer and correlate recombination
    """  gwf.target_from_template(sanify('ure_mcorr_' + group[0:5] + '_' + sample_name),
                             mcorr_bam_fit(group, sample_name, sample['forward'].strip(), sample['reverse'].strip(), reference_stem))

    # 
    gwf.target_from_template(sanify('ure_consensus_' + group[0:5] + '_' + sample_name),
                             consensus(group, sample_name, reference_stem))
    

    gwf.target_from_template(sanify('ure_extract_' + group[0:5] + '_' + sample_name),
                             extract_fasta_from_bed(group, sample_name, reference_stem))
 """


    gwf.target_from_template(sanify('ure_coverage_' + group[0:5] + '_' + sample_name),
                             coverage(group, sample_name, sample['forward'].strip(), sample['reverse'].strip(), reference_stem))


    

    

    
    
    break




### III: For each group ###