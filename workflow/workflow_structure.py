#Author Carl



from gwf import *
import subprocess, sys
import glob
import os.path

print('\n\tRemember to activate urecomb2\n')

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





title_prefix = 'Rlegum' 
genome_dir = '../Rhizobium_project_MICA/corrected_headers/consistent_order/0C_origin2'

bin_sizes = [1]


for bin_size in bin_sizes:

    genomes = glob.glob(genome_dir + '/*.xmfa')

    for genome in genomes:
        genome_basename = os.path.basename(genome)
        genome_stem = os.path.splitext(genome_basename)[0]

        title = title_prefix + '_' + str(genome_stem) + '_' + str(bin_size)

        #if not "unitig_0" in title:
        #    continue
                

        print(f"{genome} ({genome_stem})")
        
        




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
        


        gwf.target(sanify('B_split_' + title),
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
            ../../../script/xmfa_split_to_fa.py ../{file_binned_xmfa_out}.xmfa
            touch ./fasplitted

            touch ../leg_split.completed

        """
        
        

        # When the above targets are done, call it again, to run the rest of the targets below:
        splitted_xmfas = glob.glob(f"output/{title}/split/*.xmfa")
        
        for num, single_gene in enumerate(splitted_xmfas):
            single_gene_basename = os.path.basename(single_gene)
            single_gene_stem = os.path.splitext(single_gene_basename)[0]
            #print(single_gene_basename)

            '''
            gwf.target(sanify('C_mcorr_' + title + '_' + str(num) + '_' + single_gene_stem),
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
        '''
        # Run PHI
        splitted_fas = glob.glob(f"output/{title}/split/*.fa")
        
        for num, single_gene in enumerate(splitted_fas):
            single_gene_basename = os.path.basename(single_gene)
            single_gene_stem = os.path.splitext(single_gene_basename)[0]
            #print(single_gene_basename)

            gwf.target(sanify('C_PHI_' + title + '_' + str(num) + '_' + single_gene_stem),
            inputs = [single_gene],
            outputs = [f"output/{title}/split/{single_gene_stem}_phi.txt_FORCEAGAIN"],
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
            echo -e "infsites\\t\\t"$(grep -E "Found [0-9]+ informative sites\\." {single_gene_stem}_phi.txt | grep -oE "[0-9]+")"\\t\\t{single_gene_stem}" >> {single_gene_stem}_phiresult.tab



            """
            break
            
            # Clonalframe also uses .fa files, thus it can use the same single_gene* variables as PHI
            gwf.target(sanify('C_cf_' + title + '_' + str(num) + '_' + single_gene_stem),
                inputs = [single_gene],
                outputs = [f"output/{title}/split/{single_gene_stem}_cf/{single_gene_stem}_cf_final.tab"],
                cores = 1,
                walltime = '2:00:00',
                memory = '4gb',
                account = "gBGC") << f"""

                echo {single_gene_stem}     

                cd output/{title}/split/

                mkdir -p {single_gene_stem}_cf
                rm -r {single_gene_stem}_cf
                mkdir -p {single_gene_stem}_cf
                cd {single_gene_stem}_cf


                # generate guide-tree
                #mkdir -p tree
                #cd tree



                cat ../{single_gene_stem}.fa | ../../../../script/fasta_serialize_headers.py > {single_gene_stem}_serialized.fa

                raxml-ng --msa {single_gene_stem}_serialized.fa --model GTR --threads 1 --redo

                ClonalFrameML {single_gene_stem}_serialized.fa.raxml.bestTree {single_gene_stem}_serialized.fa {single_gene_stem}_cf


                # clean up output
                grep -E "(^nu)|(^1/delta)|(^R/theta)" {single_gene_stem}_cf.em.txt | awk '{{ print "{single_gene_stem}\t{title}\t" $0}}' > {single_gene_stem}_cf_final.tab


                cat {single_gene_stem}_cf_final.tab >> ../../{title}_clonalframe.tab


                # we are in the cf folder


        
                """
            #break # PHI and CF for only one gene

        # Merge the _fitpar.csv files together, so it can be imported in R later.
        fitpars = glob.glob(f"output/{title}/split/*_fitpar.csv")
        #for num, fitpar in enumerate(fitpars): # I have no idea why this job existed for each 
        

        # This target should be run when all the mcorr jobs are done
        # TODO: write this as a cat asterisk job with a checkpoint input instead af manually merging? Doesn't make sense when there is no output.
        '''
        gwf.target(sanify('D_ph2_mcorr_collect' + title),
            inputs = [i for i in fitpars],
            outputs = [f"output/{title}/{file_binned_xmfa_out}_fitpars.csv"],
            cores = 1,
            walltime = '10:00',
            memory = '1gb',
            account = "gBGC") << f"""


cd output/{title}

inputfile="split/*_fitpar.csv"

outputfile="{file_binned_xmfa_out}_fitpars.csv"

echo "" > $outputfile


if stat -t $inputfile >/dev/null 2>&1; then
    for fitpar in $inputfile; do
        awk '{{print $0, ",{bin_size}, {genome_stem}"}}' $fitpar >> $outputfile
    done
fi

            """
        '''


        # Run this target when all the PHI jobs are done.
        # TODO: Consider putting this and the previous job into a post processing script that is run manually.
        phi_results = glob.glob(f"output/{title}/split/*_phiresult.txt")
        print(len(phi_results), 'phi_results for', genome_stem)    
        gwf.target(sanify('D_ph2_PHI_collect' + title),
            inputs = [i for i in phi_results], # Det smarte her er, at hvis der er en *phiresult.txt fil der bliver opdateret, så bliver resultaterne samlet sammen igen.
            outputs = [f"output/{title}/{genome_stem}_phi_results.tab"],
            cores = 1,
            walltime = '10:00',
            memory = '1gb',
            account = "gBGC") << f"""
            cd output/{title}
            inputfile="split/*_phiresult.tab"
            outputfile="{genome_stem}_phi_results.tab"

            # reset outputfile
            echo "" > $outputfile
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

        #break # single genome


