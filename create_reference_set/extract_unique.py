import json
import numpy as np
import collections
import argparse
import lxml.etree
import functools

parser = argparse.ArgumentParser(
    prog='extract_unique'
)

parser.add_argument('--require_several_rcts', action=argparse.BooleanOptionalAction)
parser.add_argument('--require_high_quality_rcts', action=argparse.BooleanOptionalAction)

args = parser.parse_args()

@functools.cache
def is_valid_trial(study_path):
    with open(study_path) as f:
        document = lxml.etree.parse(f)
    
    if args.require_high_quality_rcts:
        design = document.find('study_design_info')
        masking = design.find('masking')
        if masking is None:
            return False
        
        if 'Participant' not in masking.text:
                return False
        
        allocation = design.find('allocation')

        if allocation.text != 'Randomized':
            return False

    return True

def is_valid_entry(entry):
    if args.require_several_rcts:
        trials = {i['study'] for i in entry['sub_infos']}
        return len(trials) > 1
    else:
        return True
    
with open('raw_entries.txt') as f:
    lines = f.readlines()
    infos = [json.loads(a) for a in lines]
    infos = [a for a in infos if is_valid_trial(a['study'])]

drug_maps = collections.defaultdict(list)

for i in infos:
    drugs = tuple(sorted([tuple(sorted(a)) for a in i['atc_codes']]))
    side_effects = tuple(sorted(list(i['icd10codes'])))
    drug_maps[(drugs, side_effects)].append(i)
    
new_infos = []

with open('unique_entries.txt', 'w') as f:

    for i, (k, v) in enumerate(drug_maps.items()):
        if i % 1000 == 0:
            print(i, len(drug_maps))
        
        rct_map = collections.defaultdict(dict)
        
        for e in v:
            for arm, drug, row in zip(e['atc_codes'], e['drugs'], e['table']):
                rct_map[(e['study'], tuple(sorted(arm)), drug)][e['side_effect_name']] = (e, row)
        
        drugs = k[0]
        
        def get_table(study):
            table = study['table']
            if drugs != tuple([tuple(sorted(a)) for a in study['atc_codes']]):
                table = [table[1], table[0]]
            return np.array(table)
        
        def collapse_rct(arm_vals):
            arm_vals = list(sorted(arm_vals, key=lambda a: a[1][0], reverse=True))

            final_count = sum(arm_vals[0][1])
            final_events = 0
            
            for e, row in arm_vals:
                if sum(row) == final_count:
                    final_events += row[0]
                
            return [final_events, final_count - final_events]
        
        final_combined_table = np.zeros((2, 2))
        
        for (rct_name, arm, drug), arm_vals in rct_map.items():
            collapsed = collapse_rct(arm_vals.values())
            if arm == drugs[0]:
                final_combined_table[0, :] += collapsed
            else:
                final_combined_table[1, :] += collapsed
        
        ratios = final_combined_table[:, 0] / np.sum(final_combined_table, axis=1)
        if ratios[0] > ratios[1]:
            final_combined_table = np.stack((final_combined_table[1, :], final_combined_table[0, :]))
            drugs = [drugs[1], drugs[0]]
        
        new_info = {
            'sub_infos': v,
            'atc_codes': drugs,
            'icd10codes': k[1],
            'table': final_combined_table.tolist(),
        }

        if is_valid_entry(new_info):
            f.write(json.dumps(new_info) + '\n')

