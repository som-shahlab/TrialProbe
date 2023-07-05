import csv
import json
import os
import math
from utils import join_reference_set_and_results


icd_map = {}
atc_map = {}

with open('smaller_mrconso.rrf') as f:
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



with open('../reference_set.txt') as f:
    infos = [json.loads(a) for a in f]

with open('observational_results.txt') as f:
    observational = [json.loads(a) for a in f]


infos = join_reference_set_and_results(infos, observational)
infos = [a for a in infos if 'postmean' in a]

infos.sort(key=lambda a:abs(a['postmean']), reverse=True)


def cell(a, b, color=None):
    if color:
        return '\\Block[fill=' + color + ']{}{' + a + ' \\\\' + b + ' }'
    else:
        return '\\Block{}{' + a + ' \\\\' + b + ' }'
    
header = '''
\\begin{NiceTabular}{ |l| l | l |  l| l | l | }
\\hline''' + f'{cell("Adverse Event", "(ICD10)")} & {cell("Drug A", "(ATC)")} & {cell("Drug B", "(ACT)")} & Unadjusted Cox & PSM Cox & IPSW Cox \Tstrut \Bstrut \\\\ \hline'

footer = '\\end{NiceTabular}'

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

    min_codes = minimize(info['icd10codes'])

    adverse_string = '\\makecell[l]{' + f'{event_name} \\\\ ({", ".join(min_codes)})' + '}'
    
    def get_drug_string(index):
        index = 1 - index
        codes = info['atc_codes'][index][:1]
        names = list({atc_map[c] for c in codes})
        
        if len(names) != 1:
            print('Wat ', info, names)
            print(1/0)
        name = capitalize(names[0])

        return '\\makecell[l]{' + f'{name} \\\\ ({", ".join(codes)})' + '}'

    def fexp(f):
        return int(math.floor(math.log10(abs(f)))) if f != 0 else 0

    def fman(f):
        return f/(10**fexp(f))

    def fformat(num):
        attempt = f'{num:0.3}'
        if 'e' not in attempt:
            return attempt
        e, n = fexp(num), fman(num)
        res = f'{n:0.3}' + '\\times 10^' + '{' + f'{e}' + '}'
        return res

    def get_value(value):
        value["p"] = float(value["p"])
        if value['coef'] < 0:
            if value["p"] < 0.05:
                color = 'green!50'
            else:
                color = None
            direction = '$A > B$'
        else:
            if value["p"] < 0.05:
                color = 'red!50'
            else:
                color = None
            direction = '$A < B$'

        first = f'\\small $p = {fformat(value["p"])}$ \\normalsize'
        return f'{cell(first, direction, color)}'


    t = [[int(a) for a in b] for b in info['table']]

    space = '\\hspace{0.2cm}'

    
    start = '\\makecell[l]{ $ \\begin{bmatrix}'
    end = '\\end{bmatrix} $ } '
    table = f'''{start} {norm(t[0][0])} & {norm(t[0][1])} \\\\
                {norm(t[1][0])} & {norm(t[1][1])} \\\\ {end}'''

    denoised_value = f'{math.exp(-info["postmean"]):.2f}'

    parts = [
        adverse_string,
        get_drug_string(0),
        get_drug_string(1),
        get_value(info["results"]["cox"]["unadjusted"]),
        get_value(info["results"]["cox"]["logistic_match"]),
        get_value(info["results"]["cox"]["logistic_ipw"]),
    ]

    return ' & '.join(parts) + ' \\\\ \\hline'

print(header)
for info in infos[:10]:
    print(get_row(info))
print(footer)
