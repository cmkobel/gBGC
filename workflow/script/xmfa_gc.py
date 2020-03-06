#!/usr/bin/env python

""" 
    This file takes all the genes from an single-line formatted .xmfa-file and calculates
    the mean GC3 content. The output is printed to STDOUT
"""

import sys, json

#Warning! This script only works for single line sequences

input_file = sys.argv[1]
try:
    cli_comment_1 = str(sys.argv[2])
    cli_comment_2 = str(sys.argv[3])
    
except Exception as e:
    print(e)
    cli_comment_1 = ""
    cli_comment2 = ""







def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)

dict = {}

print('gene', 'posA', 'posB', 'comment', 'gc_content', 'cli_comment_1', 'cli_comment_2', sep = '\t')

with open(input_file, 'r') as file:
    for line in file:
        line_s = line.strip()
        if line_s == "=":
            #print('skipping')
            del gene, right, pos, comment, posA, posB, dna, dna_gc, gc_content
            continue
        elif line_s[0] == ">":
            #eprint(line_s)
            
            gene, right = line_s[1:].split(":")
            pos, comment = right.split('+')
            posA, posB = pos.split('-')
            
            #eprint('gene', gene, 'pos', posA, posB,'comment', comment)
            

            #print(gene, isolate)

            # calculate GC content of DNA string
            
            dna = next(file)
            #print(dna[:60])

            dna_gc = 0.0
            for _i, i in enumerate(dna):
                if _i %3 == 0 and i.upper() in "GC":
                    dna_gc += 1

            gc_content = dna_gc / (float(len(dna)) / 3)


            
            print(gene, posA, posB, comment, gc_content, cli_comment_1, cli_comment_2, sep = '\t')


