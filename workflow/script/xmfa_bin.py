#!/usr/bin/env python3

"""
    This file concatenates the genes of a sorted .xmfa-file to the size of whatever is set in 'target_bin_size'.
    Outputs to STDOUT.
"""



# Author: Carl Kobel. Some ideas are taken from biopython

import sys 
import glob
import json
import re

#target_bin_size = 10   # 1K, 10K, 100K

input_file = sys.argv[1]
target_bin_size = int(sys.argv[2])
#input_file_stem = '.'.join(input_file.split(".")[0:-1])
#input_file_extension = input_file.split(".")[-1]
#output_file = f"{input_file_stem}_binsize_{target_bin_size}.{input_file_extension}

def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)

eprint(f"will concatenate sequences of {input_file} in bins of size {target_bin_size}.")
eprint(f"Outputting to STDOUT.")

XMFA_HEADER_REGEX = re.compile(r"> (?P<id>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<comments>.*)") 

current_bin = {} # could also be named current_bin
current_bin_size = 0
n_bins = 1

with open(input_file, 'r') as open_input_file:
    for line in open_input_file:
        if line[0] == "=": # The alignment ends.
            #exit()

            # The last current header is irrelevant to the current bin as a whole, and may be deleted now.
            del current_header

            
            # Check the bin size each time an alignment is gone by.
            current_bin_size = [current_bin[key]['length'] for key in current_bin.keys()][0]
            #eprint(f"current_bin_size (target_bin_size): {current_bin_size} ({target_bin_size})")

            # Når current bin-size er nået (eller overskredet) skal indholdet skrives gemmes som et enkelt alignment, og current_bin kan slettes, så der kan startes forfra med en ny bin.
            if current_bin_size >= target_bin_size:
                eprint('writing bin', n_bins, 'of size', current_bin_size)
                #eprint(json.dumps(current_bin, indent = 4))
                # Skriv bin
                for key in current_bin.keys():
                    print(f"> {key}\n{current_bin[key]['sequence']}")
                print('=')

                # Reset the current bin                
                current_bin = {}
                n_bins += 1
                
                exit()

                    
        elif line[0] == ">": # A new "fasta" header commences.
            m = re.match(XMFA_HEADER_REGEX, line)
            #eprint(m)


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

            #eprint(json.dumps(current_header, indent = 4))
        
        else: # DNA sequence data for which we should (hopefully) already have obtained the metadata.
            if not current_header['id'] in current_bin:
                current_bin[current_header['id']] = {'sequence': '', 'length': 0}



            current_bin[current_header['id']]['sequence'] += line.strip()
            current_bin[current_header['id']]['length'] += len(line.strip())
            


# The last sequence might not necessarily have been written to file.
# Let's do that.
if len(current_bin.keys()) >= 0:
    eprint('writing bin', n_bins, 'of size', current_bin_size)
    for key in current_bin.keys():
        print(f"> {key}\n{current_bin[key]['sequence']}")
        print('=')

