#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' This scripts provides the transcription factor status of the rare and common variants
from public databasess. For this project we have used Chipatlas data that have provided
with the bed files from different experiment. Provided we are looking onto the relevant
celltype or primary cell, the aim is that it would help to identify the functional
enhancer from the datasets.

Usage

python genome_tf.py P_E.txt P_E.bed DNS.CD.HUVEC.bed His_CDV_2.bed Oth.CD.HUVEC.bed -o P_E_Tf.txt


python genome_tf.py ../results/Interactor_variant/Rare_P_E.txt ../data/P_E_Interaction_Rare.bed ../data/DNS.CD.HUVEC.bed ../data/His_CDV_2.bed ../data/Oth.CD.HUVEC.bed  -o ../results/Rare_P_E_TFs.txt

'''

import sys
import argparse
import pybedtools
import numpy as np

def bed2dict(bed):
    dict_intersect = {}
    for line in bed:
        genomic_coordinate = line[0],line[1],line[2]
        dict_intersect_keys = ",".join(genomic_coordinate)
        dict_intersect_values = dict_intersect.get(dict_intersect_keys,[])
        dict_intersect_values.append(line[7])
        dict_intersect[dict_intersect_keys] = dict_intersect_values
    return dict_intersect

def compAndret(list1,list2,list3):
    A = ",".join(list1)
    if A in list2:
        test = list3[A]
        Mark = (" ".join(np.unique(test)))
    else:
        Mark = str(0)
    return(Mark)


def Main():
    parser = argparse.ArgumentParser(description="Manipulation of interaction files of rare GoTERM")
    parser.add_argument("interactionfile", help = "Interaction outputfile from vcfmanipulation script")
    parser.add_argument("interactbedfile", help = "Interaction bedfile vcfmanipulation script")
    parser.add_argument("dnase", help = "goterm downloaded from quickgo")
    parser.add_argument("hist", help = "Interaction outputfile genome and vcfmanipulation script")
    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)
    args = parser.parse_args()

    interact_bed = pybedtools.BedTool(args.interactbedfile)
    hist_bed = pybedtools.BedTool(args.hist)
    Ts_bed = pybedtools.BedTool(args.TF)
    DNAase_bed = pybedtools.BedTool(args.dnase)

    int_bed = interact_bed.intersect(hist_bed, wo = True)
    int_bed_TFs = interact_bed.intersect(Ts_bed, wo = True)
    int_DNAase = interact_bed.intersect(DNAase_bed, wo = True)

    hist_dict = bed2dict(int_bed)
    TFs_dict = bed2dict(int_bed_TFs)
    DNAse_dict = bed2dict(int_DNAase)

    keys_hist = list (hist_dict.keys())
    keys_TFs = list (TFs_dict.keys())
    keys_DNAse = list (DNAse_dict.keys())

    Result = []
    with open (args.interactionfile, 'r') as f:
        for line in f:
            A = line.strip().split("\t")
            coordinates =A[7],A[8],A[9]
            meth_mark = compAndret(coordinates, keys_hist,hist_dict)
            TFs_mark = compAndret(coordinates, keys_TFs,TFs_dict)
            DNAse_mark = compAndret(coordinates, keys_DNAse, DNAse_dict)
            res = line.strip()+"\t"+DNAse_mark+"\t"+meth_mark.strip()+"\t"+TFs_mark.strip()
            Result.append(res)

    if args.output:
        sys.stdout = open(args.output, 'w')
        for element in Result:
            print (element)
        sys.stdout.close()

if __name__ == "__main__":
    Main()
