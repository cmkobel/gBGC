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
XMFA_HEADER_REGEX = re.compile(r"> *(?P<id>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<comments>.*)") 
def parse_pop_map(structure_file):
    pop_map = {}
    t = pd.read_table(structure_file)
    t = t.rename(columns = lambda x: x.strip()) # Is this necessary?

    return t#.to_dict()



map = {}
df = parse_pop_map(structure_file)[['Seq ID', 'Origin2']]
for i, j in zip(df['Seq ID'], df['Origin2']):
    
    map[i] = j

print(json.dumps(map, indent = 2))
exit()


set_origin2 = set([i for i in df['Origin2'].values()])


# Create folders for each Origin2
for origin in set_origin2:
    try:
        os.mkdir(origin)
    except FileExistsError as fee:
        continue # not a problem at all.



# test map

map = df['Seq ID', 'Origin2']

'''



# Parse header
current_bin = {}
for line in sys.stdin:
    if line[0] == '=': #end of an alignment
        # write the data to each file
        for key, value in current_bin.items():
            print(key, value[:20], '...')
            current_bin = {}


    elif line[0] == '>': # New header commences
        m = re.match(XMFA_HEADER_REGEX, line)

        current_id = m.group('id')
        current_comments = m.group('comments')
            
        current_header = {} 
        for key in ("start", "end", "id", "strand", "comments", "realname"): 
            try: 
                value = m.group(key) 
                if key == "start": 
                    value = int(value) 
                    # Convert to zero based counting 
                    if value > 0: 
                        value -= 1 

                if key == "end": 
                    value = int(value) 
                current_header[key] = value 
            except IndexError: 
                # This will occur if we're asking for a group that 
                # doesn't exist. It's fine as long as it is not the id, which we already checked.
                pass 



        print(current_id, current_comments)
        print(json.dumps(current_header, indent = 2))




    else: # sequence data
'''
