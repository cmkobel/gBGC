#!/usr/bin/env python3


import sys 
import glob
from Bio import AlignIO
import json
import re

print(sys.argv)
input_file = sys.argv[1]
output_file = sys.argv[2]

print(input_file)
# fix malformed headers
with open(output_file, 'w') as open_output_file:
    with open(input_file, 'r') as open_input_file:
        
        for line in open_input_file:
            if line[0] == '>':
                
                old_spec, group = [i.replace('>', '').strip() for i in line.split("+")]

                seq_num, pos = old_spec.split(':')

                
                new_header = f"> {group.replace('group', '')}:{pos} + {seq_num}\n"
                #print(new_header.strip())
            
                open_output_file.write(new_header)
            else:
                open_output_file.write(line)
            


