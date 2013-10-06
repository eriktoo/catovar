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

for sample in samples:
    for variant in sample.get_variant_list():
        print ",".join(sample.get_info() + sample.get_anno(variant))

# global_variants = {}
# met_variants = []
# nomet_variants = []
# monosomy_variants = []
# disomy_variants = []
# 
# # Define global list of variables, and variable lists for sample types
# for ident in samples:
#     sample = samples[ident]
#     variants = sample.get_variants()
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