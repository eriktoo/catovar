'''
Created on Oct 6, 2013

@author: fuzz
'''
#!/usr/bin/python

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
anno_fields.extend(["Zyg", "GQ", "Freq", "Alt", "Ref", "Depth"])
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
        fao = other["FAO"].split(",")
        fao_int = 0
        for n in fao:
            fao_int += int(n)
        freq = fao_int / float(other["FDP"])
        v_data.append(str(freq)) # to do: mutation frequency
        v_data.append(other["FAO"])
        v_data.append(other["FRO"])
        v_data.append(other["FDP"])
        v_data.extend(sample.get_anno(variant)[5:]) # get remaining annotation
        data.append(v_data)

with open("combined_table.csv", 'wb') as out_csv:
    writer = csv.writer(out_csv, dialect='excel')
    writer.writerow(header)
    for row in data:
        writer.writerow(row)
            
# global_variants = {}
# met_variants = []
# nomet_variants = []
# monosomy_variants = []
# disomy_variants = []
# 
# # Define global list of variables, and variable lists for sample types
# for ident in samples:
#     sample = samples[ident]
#     variants = sample.get_variant_list()
# 
#     for variant in variants:
#         if variant in global_variants:
#             global_variants[variant].append(ident)
#         else:
#             global_variants[variant] = [ident]
# 
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