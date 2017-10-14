#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' This scripts overlaps the alleles that are present in rare frequency in Swed samples.
Considering the rare frequncy of < 1% in the swed frequency data we can estimate the allele count as 1% of 2000 = 20 and removing the private variants present only in one individual we set the lower threshold as 2.
'''


import sys
from pysam import VariantFile







Vcfin = VariantFile("/Users/salendrapradh/Documents/Rare_variants_atherosclerosis/data/swegen_20161223/SNP_VCF.gz")

with open ('/Users/salendrapradh/Documents/Rare_variants_atherosclerosis/data/BAV_HiCap/BAV_P-E_InteractionDataset.txt', 'r') as f:
    for line in f:
        line = line.strip().split("\t")
        set = [(line[9])[3:], line[10], line[11]]
        for field in line:
            print (field, end = "\t")

        for rec in Vcfin.fetch (set[0], int(set[1]), int(set[2])):
            for i in range(0, len(rec.info["AC"])):
                if rec.info["AC"][i] > 100 : # Common variant
#                if rec.info["AC"][i] > 20 and (rec.info["AC"])[i] <= 100 : # Low frequency variant
#                if rec.info["AC"][i] > 2 and (rec.info["AC"])[i] <= 20 : # Rare frequency variant
                    print (rec.chrom,rec.start,rec.stop,rec.alleles,(rec.info["AC"])[i] , end ="\t")
        print("")
#        break
