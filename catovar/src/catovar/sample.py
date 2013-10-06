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
        self.variants = {}

        # !!! Add filename checking!
        
        # Read in the multianno file for the sample
        with open(self.info_dict['filename'], 'rb') as variant_file:
            variant_dict_reader = csv.DictReader(variant_file)
            self.variant_fields = variant_dict_reader.fieldnames
            # Store dict describing each variant in variants{}, keyed by tuple of first 5 values
            for variant in variant_dict_reader:
                key_list = []
                for index in self.variant_fields[0:5]:
                    key_list.append(variant[index])
                key = tuple(key_list)
                self.variants[key] = variant

    def get_anno(self, variant):
        anno_list = []
        for key in self.variant_fields:
            anno_list.append(self.variants[variant][key])
        return anno_list

    def get_variant_fields(self):
        return self.variant_fields

    def get_variants(self):
        return self.variants.keys()

    def get_info_fields(self):
        return self.info_fields

    def get_info(self):
        info_list = []
        for key in self.info_fields:
            info_list.append(self.info_dict[key])
        return info_list