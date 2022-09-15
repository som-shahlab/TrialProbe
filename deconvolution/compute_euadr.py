import json

rows = []

with open('euadr.json') as f:
    for line in f:
        rows.append(json.loads(line))

for method in ['unadjusted', 'logistics_match', 'logistics_ipw']:
    for name in ['OMOP', 'EUADR']:
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
