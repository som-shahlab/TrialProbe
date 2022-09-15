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

def cell(a, b):
    return '\\makecell[l]{' + a + ' \\\\' + b + ' }'


header = '''
\\begin{tabular}{ | l | l |  l| l | l | l | }
\\hline \\makecell[l]{Adverse Event \\\\ (ICD10)} & \\makecell[l]{Drug A \\\\ (ATC)} & \\makecell[l]{Drug B \\\\ (ATC)} & ''' + cell("Clinical Trial", "NCTs") + ''' & \\makecell[l]{Contingency \\\\ Table} & \\makecell[l]{Denoised \\\\ Odds Ratio} \Tstrut \Bstrut \\\\
\\hline'''

footer = '\\end{tabular}'

def capitalize(a):
    return ' '.join((b[0].upper()) + b[1:] for b in a.split())

def minimize(codes):
    minimum_set = []
    for c in codes:
        if not any(c[:-1].startswith(d) for d in codes):
            minimum_set.append(c)

    return minimum_set

def norm(a):
    s = str(a)
    if len(s) < 3:
        s = '\\phantom{' + '0' * (3 - len(s)) + '}' + s
    return s

def get_row(info):
    # print(info)

    event_name = capitalize(info['sub_infos'][0]['side_effect_name'])
    rcts = list({str(int(a['study'].split('NCT')[-1].split('.')[0])) for a in info['sub_infos']})

    if len(rcts) == 1:
        rct_string = rcts[0]
    elif len(rcts) == 2:
        rct_string = cell(rcts[0], rcts[1])
    else:
        print("Wat", rcts)
        print(1/0)

    min_codes = minimize(info['icd10codes'])

    adverse_string = '\\makecell[l]{' + f'{event_name} \\\\ ({", ".join(min_codes)})' + '}'
    
    def get_drug_string(index):
        codes = info['atc_codes'][index][:1]
        names = list({atc_map[c] for c in codes})
        
        if len(names) != 1:
            print('Wat ', info, names)
            print(1/0)
        name = capitalize(names[0])

        return '\\makecell[l]{' + f'{name} \\\\ ({", ".join(codes)})' + '}'

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

    
    start = '\\makecell[l]{ $ \\begin{bmatrix}'
    end = '\\end{bmatrix} $ } '
    table = f'''{start} {norm(t[0][0])} & {norm(t[0][1])} \\\\
                {norm(t[1][0])} & {norm(t[1][1])} \\\\ {end}'''

    denoised_value = f'{math.exp(-info["postmeans"]):.2f}'

    parts = [
        adverse_string,
        get_drug_string(0),
        get_drug_string(1),
        rct_string,
        table,
        denoised_value,
    ]

    return ' & '.join(parts) + ' \\Tstrut \\Bstrut \\\\ \\hline'

print(header)
for info in infos[:10]:
    print(get_row(info))
print(footer)
