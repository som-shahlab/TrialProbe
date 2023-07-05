
import csv
import json
import os
import math

umls_dir = '/Users/ethanid/Downloads'

icd_map = {}
atc_map = {}

with open(os.path.join(umls_dir, 'MRCONSO.RRF')) as f:
    with open('smaller_mrconso.rrf', 'w') as o:
        for line in f:
            parts = line.split('|')

            source = parts[11]

            if source in ('ICD10CM', 'ICD10', 'ATC', 'MDR'):
                o.write(line)
