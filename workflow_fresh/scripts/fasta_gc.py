#!/usr/bin/env python3

""" 
    This file takes all the genes from a.fasta-file (from stdin) and calculates
    the mean GC3 content. The output is printed to stdout
"""

import sys




def eprint(*args, **kwargs):
    # I'm too lazy to write 'file = sys.stderr' manually...
    print(*args, **kwargs, file = sys.stderr)
#input_file = sys.argv[1]

cli_comment_1 = ""
cli_comment_2 = ""
cli_comment_3 = ""
cli_comment_4 = ""


try:
    cli_comment_1 = str(sys.argv[1])
    cli_comment_2 = str(sys.argv[2])
    cli_comment_3 = str(sys.argv[3])
    cli_comment_4 = str(sys.argv[4])

    
except Exception as e:
    pass
    









#print('header', 'GC3', 'cli_comment_1', 'cli_comment_2', 'cli_comment_3', sep = '\t')


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
            print(header, GC3, cli_comment_1, cli_comment_2, cli_comment_3, cli_comment_4, sep = '\t')

            # Reset variables.
            del temp_GC3, GC3, header
            sequence = ''
        

        header = line[1:]

        first_seq = False

    else:
        sequence += line