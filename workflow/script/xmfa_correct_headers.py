#!/usr/bin/env python3


import sys 
import glob
from Bio import AlignIO
import json
import re

# TODO: write to STDOUT instead of sys.argv[2]

break_xmfa_compatibility_for_the_sake_of_supporting_a_bug_in_mcorr = True

if break_xmfa_compatibility_for_the_sake_of_supporting_a_bug_in_mcorr == True:
    space_or_not = ""
else:
    space_or_not = " "

print(sys.argv)
input_file = sys.argv[1]
output_file = sys.argv[2]

print(input_file)
# fix malformed headers
with open(output_file, 'w') as open_output_file:
    open_output_file.write(f"## Warning: this is not a valid xmfa file, as the sequences within each alignment (=) is not sorted correspondingly.\n")

    with open(input_file, 'r') as open_input_file:
        open_output_file.write(f"## Info: This file was created with xmfa_correct_headers.py using {input_file} as input.\n")

        for line in open_input_file:
            if line[0] == '>':
                
                # Parse content. Consider using regex instead..
                xmfa_prefix, xmfa_comment = [i.replace('>', '').strip() for i in line.strip().split('+')]
                xmfa_seq_num, xmfa_pos = xmfa_prefix.split(':')
                
                # Map old to new format.
                gene = xmfa_comment.replace('group', '')
                pos = xmfa_pos
                isolate = xmfa_seq_num                
                
                new_header = f">{space_or_not}{gene}:{pos} + {isolate}\n"

                # Debug prints:
                # print(line.strip())
                # print(new_header)
                # print()

                open_output_file.write(new_header)
            else:
                # if line[0] == '=': print('=') # debug
                open_output_file.write(line)
            


