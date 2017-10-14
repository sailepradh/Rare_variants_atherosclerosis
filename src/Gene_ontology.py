#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
This is script to find the enrichment of Gene ontology of different proccess as
been indicated in earlier *tsv files. General format could be
as follows

python Gene_ontology.py  ../data/sorted_association.tsv  ../data/tmp.txt  -o ../data/GO_Term.txt

'''
#def count_genes(go):

import sys
import argparse
import re
import numpy as np

## Two function that to interaction and goterm files to dictionary format for
## for efficient manipulation

def go2dict(goterm):
    go = {}
    headerline = next(goterm)
    for lines in goterm:
        line = lines.strip()
        fields = line.split("\t")
        go_name = (fields[4], fields[5])
        genes_names = go.get(go_name, [])
        genes_names.append(fields[2])
        go[go_name] = genes_names
    return(go)

def unique_genes(interaction):
    interaction_genes = []
    for lines in interaction:
        line = lines.strip()
        fields = line.split("\t")
        interaction_genes.append(fields[0])
    unique_interaction_genes = np.unique(interaction_genes)
    return(unique_interaction_genes)

def Main():
    parser = argparse.ArgumentParser(description="Manipulation of interaction files of rare GoTERM")
    parser.add_argument("goterm", help = "goterm downloaded from quickgo")
    parser.add_argument("interactionfile", help = "Interaction outputfile genome and vcfmanipulation script")
    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)
    args = parser.parse_args()

    with open (args.goterm, 'r') as goterm:
        go = go2dict(goterm)
        collector = []
        for i in go.keys():
            unique_1 =len(go[i])
            res = i, np.unique(go[i]), unique_1
            collector.append(res)

    with open (args.interactionfile, 'r') as interaction:
        unique_interaction_genes = unique_genes(interaction)

    ress = {}
    for i in unique_interaction_genes:
        for j in range(len(collector)):
            if i in collector[j][1]:
                test1 = (collector[j][0],collector[j][2])
                collector3 = ress.get(test1, [])
                collector3.append(i)
                ress[test1] = collector3

    if args.output:
        sys.stdout = open(args.output, 'w')
        for i in ress.keys():
            print("\t".join(i[0]), end ="\t")
            print (i[1], end ="\t")
            print (len(ress[i]),end="\t")
            print (",".join(ress[i]), end ="\t")
            print ("")
        sys.stdout.close()

if __name__ == "__main__":
    Main()
