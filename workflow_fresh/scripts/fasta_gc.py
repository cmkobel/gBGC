#!/usr/bin/env python3

""" 
    This file takes all the genes from a.fasta-file (from stdin) and calculates
    the mean GC3 content. The output is printed to stdout
"""

import sys



#input_file = sys.argv[1]
try:
    cli_comment_1 = str(sys.argv[1])
    cli_comment_2 = str(sys.argv[2])
    cli_comment_3 = str(sys.argv[3])

    
except Exception as e:
    print(e)
    cli_comment_1 = ""
    cli_comment_2 = ""
    cli_comment_3 = ""







def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)



print('header', 'GC3', 'cli_comment_1', 'cli_comment_2', 'cli_comment_3', sep = '\t')


sequence = ''

first_seq = True
for line in sys.stdin:
    line = line.rstrip()


    if line[0] == '>':
        if first_seq == False:

            # Calculate GC3
            temp_GC3 = 0.0
            for _i, i in enumerate(sequence):
                if _i %3 == 0 and i.upper() in "GC":
                    temp_GC3 += 1
            GC3 = temp_GC3 / (float(len(sequence)) / 3)    
            
            # Output .tab format to stdout
            print(header, GC3, cli_comment_1, cli_comment_2, cli_comment_3, sep = '\t')

            # Reset variables.
            del temp_GC3, GC3, header
            sequence = ''
        

        header = line[1:]

        first_seq = False

    else:
        sequence += line
