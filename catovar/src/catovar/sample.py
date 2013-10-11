'''
Created on Oct 6, 2013

@author: fuzz
'''

#!/usr/bin/python

import csv

class Sample:
    def __init__(self, info_fields, info_dict):
        
        # !!! The first info_field MUST be 'filename' 
        
        self.info_fields = info_fields
        self.info_dict = info_dict
        self.variant_fields = []
        self.variant_list = []
        self.anno_dict = {}

        # !!! Add filename checking!
        
        # Read in the multianno file for the sample
        with open(self.info_dict['filename'], 'rb') as variant_file:
            variant_dict_reader = csv.DictReader(variant_file)
            self.variant_fields = variant_dict_reader.fieldnames
            # Store dict describing each variant in anno_dict{}
            for variant in variant_dict_reader:
                # Variants are keyed by tuple of first 5 values
                key_list = []
                for index in self.variant_fields[0:5]:
                    key_list.append(variant[index])
                key = tuple(key_list)
                # otherinfo is a tab delim field of vcf data
                # It must be parse to be useful. The first step is to append
                # the info as a list at the end of the anno list.
                variant["Otherinfo"] = variant["Otherinfo"].split('\t')
                self.anno_dict[key] = variant
                # Keep a list of all the variants present in a  sample, in the
                # order they appear in the multianno file: ascending chr, pos
                self.variant_list.append(key)

    def get_info(self):
        info_list = []
        for key in self.info_fields:
            info_list.append(self.info_dict[key])
        return info_list

    def get_info_fields(self):
        return self.info_fields
    
    def get_variant_fields(self):
        return self.variant_fields[0:-1]

    def get_variant_list(self):
        return self.variant_list
    
    def get_zyg(self, variant):
        return self.anno_dict[variant]["Otherinfo"][0]
    
    def get_gq(self, variant):
        return self.anno_dict[variant]["Otherinfo"][1]
    
    # return dictionary containing vcf info contained in the last 2 fields
    def get_vcf_info(self, variant):
        return dict(zip(self.anno_dict[variant]["Otherinfo"][-2].split(":"), \
                        self.anno_dict[variant]["Otherinfo"][-1].split(":")))
    
    # return annotation list excluding --otherinfo
    def get_anno(self, variant):
        anno_list = []
        for key in self.variant_fields[0:-1]:
            anno_list.append(self.anno_dict[variant][key])
        return anno_list