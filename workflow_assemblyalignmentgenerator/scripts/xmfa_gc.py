#!/usr/bin/env python

import sys, json

#Warning! This script only works for single line sequences

input_file = sys.argv[1]


dict = {}

with open(input_file, 'r') as file:
    for line in file:
        line_s = line.strip()
        if line_s == "=":
            #print('skipping')
            del gene, isolate, dna
            continue
        elif line_s[0] == ">":

            gene, isolate = line_s[1:].split("|")
        

            #print(gene, isolate)

            # calculate GC content of DNA string
            
            dna = next(file)
            #print(dna[:60])

            dna_gc = 0.0
            for _i, i in enumerate(dna): # TODO: expand this to use codon table.
                if _i % 3 == 0 and i.upper() in "GC":
                    dna_gc += 1

            gc_content = dna_gc / (float(len(dna))/3)

            if gene in dict:
                dict[gene].append(gc_content)
            else:
                dict[gene] = [gc_content]


#print(json.dumps(dict, indent = 4))

for key, value in dict.items():
    print(key, sum(value)/len(value), sep = '\t')


            
        
            
        

