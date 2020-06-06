import json
import numpy as np
import stats_util
import collections

with open('raw_entries.txt') as f:
    lines = f.readlines()
    infos = [json.loads(a) for a in lines]

drug_maps = collections.defaultdict(list)

for i in infos:
    drugs = tuple(sorted([tuple(sorted(a)) for a in i['atc_codes']]))
    side_effects = tuple(sorted(list(i['icd10codes'])))
    drug_maps[(drugs, side_effects)].append(i)
    
new_infos = []


with open('unique_entries.txt', 'w') as f:

    for i, (k, v) in enumerate(drug_maps.items()):
        # if i != 10092:
        #     continue

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
            # print('----')
            # print(arm_vals)

            final_count = sum(arm_vals[0][1])
            final_events = 0
            
            for e, row in arm_vals:
                if sum(row) == final_count:
                    final_events += row[0]

            # print([final_events, final_count - final_events])
                
            return [final_events, final_count - final_events]
        
        final_combined_table = np.zeros((2, 2))
        
        for (rct_name, arm, drug), arm_vals in rct_map.items():
            # print('Got from rct', rct_name, arm, drug, arm_vals)
            collapsed = collapse_rct(arm_vals.values())
            if arm == drugs[0]:
                final_combined_table[0, :] += collapsed
            else:
                final_combined_table[1, :] += collapsed
        
        ratios = final_combined_table[:, 0] / np.sum(final_combined_table, axis=1)
        if ratios[0] > ratios[1]:
            final_combined_table = np.stack((final_combined_table[1, :], final_combined_table[0, :]))
            drugs = [drugs[1], drugs[0]]
        
        # Compute minimum p values
        def compute_minimim_null_p_value(table):
            total_events = np.sum(table, axis=0)[0]
            counts = np.sum(table, axis=1)
            
            total_p = total_events / sum(counts)

            def helper(counts):
                first_events = min(int(total_p * counts[0] + 0.5), counts[0])
                second_events = total_events - first_events
                
                new_table = [
                    [first_events, counts[0] - first_events],
                    [second_events, counts[1] - second_events],
                ]

                # print(new_table)
                
                return stats_util.get_p_value_null(new_table)

            return min(helper(counts), helper(counts[::-1])) 
        
        def compute_minimim_positive_p_value(table):
            total_events = np.sum(table, axis=0)[0]
            counts = np.sum(table, axis=1)

            def helper(counts):
                most_events = min(total_events, counts[1])
                other_events = total_events - most_events
                
                new_table = [
                    [other_events, counts[0] - other_events],
                    [most_events, counts[1] - most_events],
                ]
                # print(new_table)

                return stats_util.get_p_value_positive(new_table)

            return min(helper(counts), helper(counts[::-1]))

        # print(i, v)
        # print(drugs)
        # print(final_combined_table.tolist())

        new_info = {
            'sub_infos': v,
            'atc_codes': drugs,
            'icd10codes': k[1],
            'table': final_combined_table.tolist(),
            'null_p_value': stats_util.get_p_value_null(final_combined_table),
            'positive_p_value': stats_util.get_p_value_positive(final_combined_table),
            'min_null_p_value': compute_minimim_null_p_value(final_combined_table),
            'min_positive_p_value': compute_minimim_positive_p_value(final_combined_table),
        }

        f.write(json.dumps(new_info) + '\n')

