import csv
import json
import os
import math


icd_map = {}
atc_map = {}

with open('smaller') as f:
    for line in f:
        if False:
            break
        parts = line.split('|')

        source = parts[11]

        if source == 'ICD10CM':
            m = icd_map
        elif source == 'ATC':
            m = atc_map
        else:
            continue

        code = parts[13]
        text = parts[14]

        m[code] = text



with open('results.json') as f:
    infos = [json.loads(a) for a in f]

post_means = {}

with open('with_means.csv') as f:
    reader = csv.DictReader(f)
    for line in reader:
        key = tuple(int(line[k]) for k in ('X1', 'N1', 'X2', 'N2'))
        post_means[key] = float(line['postmeans'])
        
for info in infos:
    key = tuple(int(info['table'][a][b]) for a in range(2) for b in range(2))
    if key not in post_means:
        if key[0] == 0 or key[2] == 0:
            info['postmeans'] = 0
        else:
            print(info, key)
            print(1/ 0)
    else:
        info['postmeans'] = post_means[key]

    

infos.sort(key=lambda a:abs(a['postmeans']), reverse=True)



header = '''

\\begin{center}
\\begin{tabular}{ p{3cm}  p{3cm}  p{3cm} p{2cm} p{2cm} }
 Adverse Event (ICD10) & Drug A (ATC) & Drug B (ATC) & Contingency Table & Denoised odds ratio \\\\
\\hline

'''

footer = '''
\\end{tabular}
\\end{center}
'''

def capitalize(a):
    return ' '.join((b[0].upper()) + b[1:] for b in a.split())

def minimize(codes):
    minimum_set = []
    for c in codes:
        if not any(c[:-1].startswith(d) for d in codes):
            minimum_set.append(c)

    return minimum_set

def get_row(info):
    # print(info)

    event_name = capitalize(info['sub_infos'][0]['side_effect_name'])

    min_codes = minimize(info['icd10codes'])

    adverse_string = f'{event_name} \\newline ({", ".join(min_codes)})'
    
    def get_drug_string(index):
        codes = info['atc_codes'][index]
        names = list({atc_map[c] for c in codes})
        
        if len(names) != 1:
            print('Wat ', info, names)
            print(1/0)
        name = capitalize(names[0])

        return f'{name} \\newline ({", ".join(codes)})'

    def get_value(value):
        if value['coef'] < 0:
            if value["p"] < 0.05:
                color = '\\cellcolor{green!25}'
            else:
                color = ''
            direction = '$B > A$'
        else:
            if value["p"] < 0.05:
                color = '\\cellcolor{red!25}'
            else:
                color = ''
            direction = '$B < A$'

        return f'{color} $p = {value["p"]:.3}$ \\newline {direction}'

    t = [[int(a) for a in b] for b in info['table']]

    space = '\\hspace{0.2cm}'

    table = f'''
        $\\lceil {t[0][0]} {space} {t[0][1]} \\rceil$ \\newline
        $\\lfloor {t[1][0]} {space} {t[1][1]} \\rfloor$
    '''

    denoised_value = f'{math.exp(-info["postmeans"]):.2f}'

    parts = [
        adverse_string,
        get_drug_string(0),
        get_drug_string(1),
        table,
        denoised_value,
    ]

    return ' & '.join(parts) + '\\\\'

print(header)
for info in infos[:10]:
    print(get_row(info))
print(footer)
