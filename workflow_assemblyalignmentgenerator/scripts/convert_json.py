#!/usr/bin/env python3

"""
    This file takes the .json-file outputted from mcorr-xmfa and prints the contents
    in a tab-delimited format to STDOUT.
"""

import json, sys


input_file = sys.argv[1]

print('gene\tgroup\tlag\tmean\tvar\tN\ttype')


with open(input_file, 'r') as file:
    for line in file:
        dic = json.loads(line)
        
        


        id = dic['ID'].split('|')
        gene, isolate = id

        for i in dic['Results']:
            print(gene, isolate, end = '\t')
            for j, k in i.items():
                print(k, end = '\t')
            print()
        