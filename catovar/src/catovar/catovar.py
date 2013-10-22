#!/usr/bin/python
'''
Created on Oct 6, 2013

@author: fuzz
'''

import csv
import sample
import sys

# Define variables
samples = []

# Read in sample definitions, create sample instances, and store in samples[]
with open(sys.argv[1]) as def_file:
    dictreader = csv.DictReader(def_file)
    info_fields = dictreader.fieldnames
    for row in dictreader:
        new_sample = sample.Sample(info_fields, row)
        samples.append(new_sample)
        
# Construct a header list for CSV output
anno_fields = samples[0].get_variant_fields()[0:5]
anno_fields.extend(["Zyg", "GQ", "Freq", "Alt", "% Alt+", "% Alt-", "Ref", "Depth"])
anno_fields.extend(samples[0].get_variant_fields()[5:])
header = samples[0].get_info_fields() + anno_fields

data = [] # list to hold lists of variant data
for sample in samples:
    for variant in sample.get_variant_list():
        v_data = [] # initialize list to line of data for variant
        other = sample.get_vcf_info(variant) # get dict holding GT:GQ:DP... (vcf form)
        v_data.extend(sample.get_info()) # get sample information
        v_data.extend(sample.get_anno(variant)[0:5]) # get first 5 anno fields
        v_data.append(sample.get_zyg(variant)) # insert zyg, qual, and coverage
        v_data.append(sample.get_gq(variant))
        
        # caculate mutation frequency
        fao = other["FAO"].split(",")
        fao_int = 0
        for n in fao:
            fao_int += int(n)
        try:
            freq = fao_int / float(other["FDP"])
            v_data.append(str(freq))
        except ZeroDivisionError:
            v_data.append("NA")
            
        v_data.append(other["FAO"])

        # calculate % of alt reads on + strand
        fsaf = other["FSAF"].split(",")
        fsaf_int = 0
        for n in fsaf:
            fsaf_int += int(n)
        try:
    #        print "fsaf_int " + str(fsaf_int)
    #        print "fao_int " + str(fao_int)
            fsaf_freq = fsaf_int / float(fao_int)
            v_data.append(str(fsaf_freq))
        except ZeroDivisionError:
            v_data.append("NA")
        
        # calculate % of alt reads on - strand
        fsar = other["FSAR"].split(",")
        fsar_int = 0
        for n in fsar:
            fsar_int += int(n)
        try:
            fsar_freq = fsar_int / float(fao_int)
            v_data.append(str(fsar_freq))
        except ZeroDivisionError:
            v_data.append("NA")

        v_data.append(other["FRO"])
        v_data.append(other["FDP"])
        v_data.extend(sample.get_anno(variant)[5:]) # get remaining annotation
        data.append(v_data)

with open("combined_table.csv", 'wb') as out_csv:
    writer = csv.writer(out_csv, dialect='excel')
    writer.writerow(header)
    for row in data:
        writer.writerow(row)
            
global_variants = {}
global_fields = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene", "1000g2012apr_all", "snp132", "cosmic65"]

met_variants = []
nomet_variants = []
monosomy_variants = []
disomy_variants = []
 
# Define global variant lookup table
for sample in samples:
    variants = sample.get_variant_list()

    for variant in variants:
        if variant in global_variants:
            continue
        else:
            anno = sample.get_anno(variant)
            

#     met_status = sample.info_dict["metast"]
#     somy_status = sample.info_dict["somy"]
# 
#     if met_status == "mets":
#         met_variants.extend(variants)
#     elif met_status == "nomets":
#         nomet_variants.extend(variants)
#     else:
#         print("Error: illegal metastasis value.")
# 
#     if somy_status == "monosomy":
#         monosomy_variants.extend(variants)
#     elif somy_status == "disomy":
#         disomy_variants.extend(variants)
#     else:
#         print("Error: illegal somy value.")

# fieldnames = 
# with open('output.csv', 'rb') as outfile:
#     writer = csv.DictWriter(outfile, 

# Various print calls to test output

# for variant in global_variants:
#     for sample in global_variants[variant]:
#         print(samples[sample].get_info())
#         print(samples[sample].get_anno(variant))
            
#print("GLOBAL VARIANTS")
#print set(global_variants)
#print global_variants.keys()

#print("METS VARIANTS")
#print set(met_variants)

#print("NOMETS VARIANTS")
#print set(nomet_variants)

#print("MONOSOMY VARIANTS")
#print set(monosomy_variants)

#print("DISOMY VARIANTS")
#print set(disomy_variants)
