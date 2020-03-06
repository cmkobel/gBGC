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
        # Lower case letters, upper case letters and numbers. And underscore
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


title_prefix = 'Rleg'

genome_dir = '../Rhizobium_project_MICA/corrected_headers/consistent_order'



#title = title_prefix + str(bin_size)
genomes = glob.glob(genome_dir + '/*')

for genome in genomes:
    genome_basename = os.path.basename(genome)
    genome_stem = os.path.splitext(genome_basename)[0]
    title = title_prefix + '_' + str(genome_stem) + '_' + '1'
    print(f"{genome} ({genome_stem})")
    
    gwf.target(sanify('a_gc_and_split_' + title),
        inputs = f"{genome}",
        outputs = [f"output/{title}/{genome_stem}_gc.tab"],
        cores = 1,
        walltime = '0:10:00',
        memory = '1gb',
        account = "gBGC") << f"""

    
    mkdir -p output/{title}
    cd output/{title}

    # compute GC from the alignment
    ../../script/xmfa_gc.py ../../{genome} NA {genome_stem} > {genome_stem}_gc.tab


    # split the genome into genes
    mkdir -p split
    cd split
    ../../../script/xmfa_split_to_fa.py ../../../{genome}

    """

    splitted_fas = glob.glob(f"output/{title}/split/*.fa")
            
    for num, single_gene in enumerate(splitted_fas):
        single_gene_basename = os.path.basename(single_gene)
        single_gene_stem = os.path.splitext(single_gene_basename)[0]
        #print(single_gene_basename)

        gwf.target(sanify('b_phi_' + title + '_' + str(num) + '_' + single_gene_stem),
            inputs = [single_gene],
            outputs = [f"output/{title}/split/{single_gene_stem}_phi.txt"],
            cores = 1,
            walltime = '2:00:00',
            memory = '4gb',
            account = "gBGC") << f"""

    cd output/{title}/split
    
    echo {single_gene_stem}
    mkdir -p "{single_gene_stem}_temp"   
    cd {single_gene_stem}_temp

    ../../../../script/PhiPack/Phi -f ../{single_gene_basename} -p 1000 -o > ../{single_gene_stem}_phi.txt 2> ../{single_gene_stem}_phierr.txt
    # above line fails when there is not enough informative sites.

    #check if stderr is empty
    if [ $(stat -c%s ../{single_gene_stem}_phierr.txt) -eq 0 ]; then
        echo "empty"
        tail -n 5 ../{single_gene_stem}_phi.txt > ../{single_gene_stem}_phiresult.txt
    else
        echo "not empty"
    fi

    cd ..
    rm  {single_gene_stem}_phierr.txt
    rm -r {single_gene_stem}_temp

    # collect results in tab file 
    cat {single_gene_stem}_phiresult.txt | sed 's/NSS/NSS NA/g' | grep -E "^NSS|^Max|^PHI" | awk '{{print $1 "\t" $2 "\t" $3 "\t{genome_stem}\t{single_gene_stem}"}}' > {single_gene_stem}_phiresult.tab

    #printf {single_gene_stem}, >> {single_gene_stem}_fitpar.csv
    #awk '/^all,/' {single_gene_stem}_fit_results.csv >> {single_gene_stem}_fitpar.csv


            """
        #if num >= 10:
        #    break  # first 10 genes
    


    phi_results = glob.glob(f"output/{title}/split/*_phiresult.txt")
    print(len(phi_results), 'phi_results for', genome_stem)    
    gwf.target(sanify('c_collect_results' + title),
        inputs = [i for i in phi_results], # Det smarte her er, at hvis der er en *phiresult.txt fil der bliver opdateret, så bliver resultaterne samlet sammen igen.
        outputs = [f"output/{title}/{genome_stem}_phi_results.tab"],
        cores = 1,
        walltime = '10:00',
        memory = '1gb',
        account = "gBGC") << f"""

    cd output/{title}

    inputfile="split/*_phiresult.tab"
    outputfile="{genome_stem}_phi_results.tab"

    if stat -t $inputfile >/dev/null 2>&1; then
        echo inside

        echo -e "method\tdetail\tpvalue\tgenome\tgene" > $outputfile
        cat $inputfile >> $outputfile
        
        #for resultfile in $inputfile; do
        #    #awk '{{print $0, ",bin_size, {genome_stem}"}}' $resultfile >> $outputfile
        #    cat 
        #done

    fi
        """

        #break # first genome