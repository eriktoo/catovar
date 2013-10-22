#!/usr/bin/python
'''
Created on Oct 6, 2013

@author: fuzz
'''

import csv
import sample
import sys

# Define functions

# Write list 'data' into csv called 'filename' with fields from 'header'
def output(filename, header, data):
    with open(filename, 'wb') as out_csv:
        writer = csv.writer(out_csv, dialect='excel')
        writer.writerow(header)
        for row in data:
            writer.writerow(row)

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

output("combined.csv", header, data)
            
global_variants = {}
global_anno_fields = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene", "1000g2012apr_all", "snp132", "cosmic65"]

mets_variants = []
nomets_variants = []
monosomy_variants = []
disomy_variants = []
 
# Define global variant lookup table
for sample in samples:
    variants = sample.get_variant_list()

    for variant in variants:
        if variant in global_variants:
            continue
        else:
            anno = sample.get_anno_dict(variant)
            global_anno = []
            for field in global_anno_fields:
                global_anno.append(anno[field])
            global_variants[variant] = global_anno
            
# Get sample metastasis and somy status 
    met_status = sample.get_mets()
    somy_status = sample.get_somy()

    if met_status == "mets":
        mets_variants.extend(variants)
    elif met_status == "nomets":
        nomets_variants.extend(variants)
    else:
        continue
 
    if somy_status == "monosomy":
        monosomy_variants.extend(variants)
    elif somy_status == "disomy":
        disomy_variants.extend(variants)
    else:
        continue

# Build sets based on mets and somy
mets_set = set(mets_variants)
nomets_set = set(nomets_variants)
monosomy_set = set(monosomy_variants)
disomy_set = set(disomy_variants)

# Build sets of intersections of mets and somy
mets_mono = mets_set & monosomy_set
nomets_mono = nomets_set & monosomy_set
mets_di = mets_set & disomy_set
nomets_di = nomets_set & disomy_set

# Various print calls to test output
# for row in data:
#     print(row)
#             
# print("GLOBAL VARIANTS")
# for key in sorted(global_variants.keys()):
#     print global_variants[key]
# 
# print("METS VARIANTS")
# for variant in mets_set:
#     print variant
#   
# print("NOMETS VARIANTS")
# for variant in nomets_set:
#     print variant
#   
# print("MONOSOMY VARIANTS")
# for variant in monosomy_set:
#     print variant
#   
# print("DISOMY VARIANTS")
# for variant in disomy_set:
#     print variant