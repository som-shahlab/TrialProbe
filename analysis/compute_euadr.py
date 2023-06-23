import json
import csv

rows = []

with open('euadr.json') as f:
    for line in f:
        rows.append(json.loads(line))

for name in ['OMOP', 'EUADR']:
    for method in ['unadjusted', 'logistics_match', 'logistics_ipw']:
        num_total = 0
        num_found = 0
        num_correct = 0

        for row in rows:
            if not row['study'].startswith(name):
                continue
            if row['direction'] == '=':
                continue

            data = row['cox'][method]

            num_total += 1

            significant = data['p'] < 0.05
            if significant:
                num_found += 1

                if data['coef'] < 0:
                    num_correct += 1

        print(name, method, num_total, num_correct / num_found, num_correct / num_total)

with open('with_means.csv') as f:
    data = list(csv.DictReader(f))

for name in ['Ours']:
    for method in ['unadjusted', 'match', 'ipw']:
        num_total = 0
        num_found = 0
        num_correct = 0

        for row in data:
            if float(row['p_values']) > 0.05:
                continue 

            p = float(row[method + '_p'])
            coef = float(row[method + '_Z'])

            num_total += 1

            significant = p < 0.05
            if significant:
                num_found += 1

                if coef < 0:
                    num_correct += 1

        print(name, method, num_total, num_correct / num_found, num_correct / num_total)
