#!/usr/bin/env python3

help_text = """
This file concatenates the genes of a sorted .xmfa-file to the size of whatever is set in 'target_bin_size'.
Outputs to STDOUT.

Example:
    ./xmfa_bin.py genome.xmfa 1000
will bin the genes in genome.xmfa into bins of minimum 1000 basepairs
"""



# Author: Carl Kobel. Some ideas are taken from biopython

import sys 
import glob
import json
import re


def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)

try:
    input_file = sys.argv[1]
    target_bin_size = int(sys.argv[2]) # 1K, 10K, 100K
except Exception as e:
    eprint('Exception:', e)
    eprint(help_text)
    exit()

eprint(f"will concatenate sequences of {input_file} in bins of size {target_bin_size}.")
eprint(f"Outputting to STDOUT.")
eprint()

XMFA_HEADER_REGEX = re.compile(r"> (?P<id>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<comments>.*)") 

current_bin = {} # could also be named current_bin
current_bin_size = 0
n_bins = 0
number = 0
last_end = 1


with open(input_file, 'r') as open_input_file:
    for line in open_input_file:
        if line[0] == "=": # The alignment ends.
            #exit()

            # The last current header is irrelevant to the current bin as a whole, and may be deleted now.
            del current_header
            number = 0

            # Check that all sequences in the current bin have the same length.
            current_bin_sizes = list(set([len(current_bin[key]) for key in current_bin.keys()]))
            if len(current_bin_sizes) > 1:
                raise Exception(f"Error: alignments of {current_id} do not have equal lengths. Observed lengths: {current_bin_sizes}.")
            current_bin_size = current_bin_sizes[0]
            eprint(f"id {current_id}: size {current_bin_size}")


            # Når current bin-size er nået (eller overskredet) skal indholdet skrives gemmes som et enkelt alignment, og current_bin kan slettes, så der kan startes forfra med en ny bin.
            if current_bin_size >= target_bin_size:
                n_bins += 1
                eprint(f"writing bin {n_bins} of size {current_bin_size}\n")
                #eprint(json.dumps(current_bin, indent = 4))
                # Skriv bin
                for key in sorted(current_bin.keys()):
                    print(f"> {current_id}:{last_end}-{last_end+current_bin_size-1} + bin {n_bins}\n{current_bin[key]}")
                print('=')
                
                last_end += current_bin_size

                # Reset 
                current_bin = {}
                del current_id
                
                del current_bin_sizes
                del current_bin_size

                    
        elif line[0] == ">": # A new "fasta" header commences.
            m = re.match(XMFA_HEADER_REGEX, line)
            number += 1
            #eprint(m)


            # Let's check that the id actually exists.
            current_id = m.group('id')
            #eprint('id', current_id)

            # Read metadata
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

            if not number in current_bin:
                current_bin[number] = "" 

            current_bin[number] += line.strip()
            
            

# The last sequence might not necessarily have been written to file.
# Let's do that.
if len(current_bin.keys()) > 0:
    n_bins += 1
    eprint('writing last bin', n_bins, 'of size', current_bin_size)
    #eprint(json.dumps(current_bin, indent = 4))
    # Skriv bin
    for key in sorted(current_bin.keys()):
        #print(f"> {current_id}\n{current_bin[key]}")
        #print(f"> {current_id}:{last_end}-{last_end+current_bin_size-1}\n{current_bin[key]}")
        print(f"> {current_id}:{last_end}-{last_end+current_bin_size-1} + bin {n_bins}\n{current_bin[key]}")
    print('=')














