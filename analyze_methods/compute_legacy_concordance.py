import json
import csv
from utils import join_reference_set_and_results

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



with open('../reference_set.txt') as f:
    infos = [json.loads(a) for a in f]

with open('observational_results.txt') as f:
    observational = [json.loads(a) for a in f]

print(len(infos))

infos = join_reference_set_and_results(infos, observational)
print(len(infos))
infos = [a for a in infos if 'postmean' in a]
print(len(infos))

infos.sort(key=lambda a:abs(a['postmean']), reverse=True)

for name in ['Ours']:
    for method in ['unadjusted', 'logistic_match', 'logistic_ipw']:
        num_total = 0
        num_found = 0
        num_correct = 0

        for info in infos:
            if info['p'] > 0.05:
                continue 

            p = info["results"]["cox"][method]['p']
            coef = info["results"]["cox"][method]['coef']

            num_total += 1

            significant = p < 0.05
            if significant:
                num_found += 1

                if coef < 0:
                    num_correct += 1

        def h(a):
            return f'{100 * a: 0.0f}'

        print(name, method, num_total, h(num_correct / num_found), h(num_correct / num_total))
