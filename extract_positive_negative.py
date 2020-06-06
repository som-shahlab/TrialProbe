import json
import numpy as np

with open('unique_entries.txt') as f:
    lines = f.readlines()
    infos = [json.loads(a) for a in lines]

def get_items(alpha, min_k, k, direction):
    filtered_infos = [a for a in infos if a[min_k] < alpha]

    infos_by_p_value_null = list(sorted(filtered_infos, key=lambda i: i[k]))

    baseline = np.array(list(range(1, len(infos_by_p_value_null) + 1))) / len(infos_by_p_value_null) * alpha

    significant = [i[k] < b for i,b in zip(infos_by_p_value_null, baseline)]

    index = len(significant) - 1 - list(reversed(significant)).index(True)
    print('Found ', index + 1, 'for', direction)
    
    for i in range(index+1):
        b = dict(infos_by_p_value_null[i])
        b['direction'] = direction
        yield b

negatives = list(get_items(0.05, 'min_null_p_value', 'null_p_value', '='))
positives = list(get_items(0.025, 'min_positive_p_value', 'positive_p_value', '<'))
        
with open('significant_entries.txt', 'w') as f:
    for i in negatives + positives:
        f.write(json.dumps(i) + '\n')