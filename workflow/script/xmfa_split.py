#!/usr/bin/env python3

"""
    This file splits each gene in an .xmfa-file into individual .xmfa-files.
"""

import sys


input_file = sys.argv[1]

buffer = ""

with open(input_file, 'r') as file:
    for line in file:
        #line = line.strip()
        buffer += line
        
        if line[0] == '=':
            with open(f'{gene.strip()}.fa', 'w') as output_file:
                #buffer += '\n' # apparently necessary for mcorr-xmfa to work. Who would have known..
                output_file.write(buffer)

            del gene
            buffer = ""
            
            continue # try and delete this line

        elif line[0] == '>': # new header
            gene, right = line[1:].split(":")
            pos, comment = right.split('+')
            posA, posB = pos.split('-')

        
        
        
