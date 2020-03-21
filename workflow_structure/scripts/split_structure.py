#!/usr/bin/env python3

import sys
import pandas as pd
import json
import os
import re


# Help
structure_file = sys.argv[1]
print('This is the structure split script.', file = sys.stderr)
print('Structure file:', structure_file, file = sys.stderr)


# Field variables and functions
CXMFA_HEADER_REGEX = re.compile(r"> *(?P<gene>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) isolate (?P<seqid>\d+)-(?P<strainid>[0-9a-zA-Z]+) bin (?P<bin>\d+)") 

def parse_pop_map(structure_file):
    pop_map = {}
    t = pd.read_table(structure_file)
    t = t.rename(columns = lambda x: x.strip()) # Is this necessary?

    return t

def eprint(string, *args, **kwargs):
    print(string, *args, **kwargs, file = sys.stderr)
    


# Create a dictionary that maps seqids to their respective countries
map = {}
print('This is the map used for structuring the samples:')
df = parse_pop_map(structure_file)[['Seq ID', 'Origin2']]
for i, j in zip(df['Seq ID'], df['Origin2']):
    
    i = str(i)
    j = str(j)
    map[i] = j

eprint(json.dumps(map, indent = 2))



try:
    os.mkdir('split_structure')
except FileExistsError as _:
    pass


print('\nparsing xmfa from stdin:\n')

# Parse header
current_alignment = {}#{j:'' for j in set_origin2}
for line in sys.stdin:
    
    if line[0] == '=': #end of an alignment
        #eprint('\nend of alignment\n')

        #eprint(json.dumps(current_alignment, indent = 2))



        # write the data to each file
        for key, value in current_alignment.items():
            #print(key, str(value)[0:1000]) # For now it just prints out to a single file.
            with open(f"split_structure/{key}.xmfa", 'a') as output_file:
                output_file.write(str(value) + '=\n')







        # Reset current_alignment
        for key in current_alignment:
            current_alignment[key] = ''


    elif line[0] == '>': # New header commences
        m = re.match(CXMFA_HEADER_REGEX, line)
        #eprint(line.strip())

        # These keys are strictly necessary
        current_gene = m.group('gene')
        current_seqid = m.group('seqid')
        current_line = line
        

        current_header = {} 
        for key in ('gene', 'start', 'end', 'strand', 'seqid', 'strainid', 'bin'): 
            try: 
                value = str(m.group(key))
                if key == 'start': 
                    value = int(value) 
                    # Convert to zero based counting 
                    if value > 0: 
                        value -= 1 

                if key == "end": 
                    value = int(value) 
                current_header[key] = value 
            except IndexError: 
                # This will occur if we're asking for a group that 
                # doesn't exist. It's fine as long as it is not the gene, which we already checked.
                pass 



        #eprint(map[current_seqid])
        #eprint(json.dumps(current_header, indent = 2))

        if not map[current_seqid] in current_alignment:
            current_alignment[map[current_seqid]] = ''

        current_alignment[map[current_seqid]] += line



    elif line[0] == '#': # comment, skip
        continue

    else: # DNA
        current_alignment[map[current_seqid]] += line 


# TODO: remember to write the last alignment.