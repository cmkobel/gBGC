#!/usr/bin/env python3

# Author: Carl Kobel. Some ideas are taken from biopython

import sys 
import glob
import json
import re

target_bin_size = 7 # 1K, 10K, 100K

input_file = sys.argv[1]
#input_file_stem = '.'.join(input_file.split(".")[0:-1])
#input_file_extension = input_file.split(".")[-1]
#output_file = f"{input_file_stem}_binsize_{target_bin_size}.{input_file_extension}


XMFA_HEADER_REGEX = re.compile(r"> (?P<id>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<comments>.*)") 

def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)

concatenated_alignments = {} # could also be named current_bin
current_bin_size = 0

with open(input_file, 'r') as open_input_file:
    for line in open_input_file:
        if line[0] == "=": # The alignment ends.
            #exit()

            # The last current header is irrelevant to the current bin as a whole, and may be deleted now.
            del current_header

            
            # Check the bin size each time an alignment is gone by.
            current_bin_size = [concatenated_alignments[key]['length'] for key in concatenated_alignments.keys()][0]
            eprint(f"current_bin_size (target_bin_size): {current_bin_size} ({target_bin_size})")
            eprint('current_bin_size:', current_bin_size)

            # Når current bin-size er nået (eller overskredet) skal indholdet skrives gemmes som et enkelt alignment, og concatenated_alignments kan slettes, så der kan startes forfra med en ny bin.
            if current_bin_size >= target_bin_size:
                eprint('writing bin!')
                # Skriv bin
                for key in concatenated_alignments.keys():
                    print(f"> {key}\n{concatenated_alignments[key]['sequence']}")
                print('=')

                
                concatenated_alignments = {}
                # Reset the bin size
                #current_bin_size = 0
                    
        elif line[0] == ">": # A new "fasta" header commences.
            m = re.match(XMFA_HEADER_REGEX, line)
            eprint(m)


            # Let's check that the id actually exists.
            parsed_id = m.group('id')
            eprint(parsed_id)

            
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

            eprint(json.dumps(current_header, indent = 4))
        
        else: # DNA sequence data for which we should (hopefully) already have obtained the metadata.
            if not current_header['id'] in concatenated_alignments:
                concatenated_alignments[current_header['id']] = {'sequence': '', 'length': 0}



            concatenated_alignments[current_header['id']]['sequence'] += line.strip()
            concatenated_alignments[current_header['id']]['length'] += len(line.strip())
            
        

eprint(json.dumps(concatenated_alignments))