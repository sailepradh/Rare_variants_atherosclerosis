#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' This scripts overlaps the alleles that are present in rare frequency in
Swed samples.Considering the rare frequncy of < 1% in the swed frequency data
we can estimate the allele count as 1% of 2000 = 20 and removing the private
variants present only in one individual we set the lower threshold as 2.
Usage:

Current usage of script is
python test_2.py [VCFfile] [Interactiondata] (class_of_vriant) -o filename.txt

python vcfmanipulation.py ../data/swegen_20161223/SNP_VCF.gz ../data/BAV_HiCap/BAV_P-E_InteractionDataset.txt  Lowfreq  --output ../results/Interactor_variant/Lowfreq_P_E.txt
python vcfmanipulation.py ../data/swegen_20161223/SNP_VCF.gz ../data/BAV_HiCap/BAV_P-E_InteractionDataset.txt  Common  --output ../results/Interactor_variant/Common_P_E.txt
python vcfmanipulation.py ../data/swegen_20161223/SNP_VCF.gz ../data/BAV_HiCap/BAV_P-E_InteractionDataset.txt  Rare  --output ../results/Interactor_variant/Rare_P_E.txt
'''
import sys
import argparse
import sys
import re
from pysam import VariantFile

def allecount(ac):
    if ac < 20 and ac > 2:
        return ("Rare")
    elif ac< 100 and ac >= 20:
        return ("Lowfreq")
    elif ac > 100:
        return ("Common")

def printformat(dictt):
    for keys, values in dictt.items():
        print (keys, end= "\t")
        counter = 0
        for i in values:
            counter = counter+1
            print (i, end =":")
        print (counter)

pattern = re.compile('chr.')

def Main():
    parser = argparse.ArgumentParser(description="loading vcf and interaction files")
    parser.add_argument("vcfile", help = "Vcf file from swedegen")
    parser.add_argument("interactionfile", help = "Interaction outputfile from HiCap methods")
    parser.add_argument("freq", choices = ['Rare','Lowfreq','Common'],help = "choose the class")
    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)
    args = parser.parse_args()
    result= {}
    Vcfin = VariantFile(args.vcfile)

    ## assert to supply only the vcf files

    with open (args.interactionfile, 'r') as f:
        for line in f:
            line = line.strip().split("\t")
            position1 = [line[9], line[10], line[11]]
            position2 = [line[13], line[14], line[15]]

            if pattern.match(position1[0]):
                chr =  ((line[9])[3:],line[10], line[11])
                fields = (line[0],line[1],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12])
            elif pattern.match(position2[0]):
                chr = ((line[13])[3:], line[14], line[15])
                fields = (line[0],line[1],line[4],line[5],line[6],line[7],line[8],line[13],line[14],line[15],line[18],line[9],line[10],line[16],line[17])

            for rec in Vcfin.fetch (chr[0], int(chr[1]), int(chr[2])):
                for i in range(0, len(rec.info["AC"])):

                    if allecount(rec.info["AC"][i]) == "Rare" and args.freq == "Rare":
                        result_keys = "\t".join(fields)
                        result_values = result.get(result_keys,[])
                        res = [rec.chrom,rec.stop,rec.alleles,(rec.info["AC"])[i]]
                        result_values.append(",".join(str(x) for x in res))
                        result[result_keys] = result_values

                    if allecount(rec.info["AC"][i]) == "Lowfreq" and args.freq == "Lowfreq":
                        result_keys = "\t".join(fields)
                        result_values = result.get(result_keys,[])
                        res = [rec.chrom,rec.stop,rec.alleles,(rec.info["AC"])[i]]
                        result_values.append(",".join(str(x) for x in res))
                        result[result_keys] = result_values

                    elif allecount(rec.info["AC"][i]) == "Common" and args.freq == "Common":
                        result_keys = "\t".join(fields)
                        result_values = result.get(result_keys,[])
                        res = [rec.chrom,rec.stop,rec.alleles,(rec.info["AC"])[i]]
                        result_values.append(",".join(str(x) for x in res))
                        result[result_keys] = result_values

    if args.output:
        sys.stdout = open(args.output, 'w')
        printformat(result)
        sys.stdout.close()

if __name__ == "__main__":
    Main()
