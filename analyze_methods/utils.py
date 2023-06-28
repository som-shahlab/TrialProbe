import copy

def create_map(entries):
    result = {}
    for entry in entries:
        icdcodes = tuple(sorted(list(entry['icd10codes'])))
        atc_codes = tuple(sorted([tuple(sorted(arm)) for arm in entry['atc_codes']]))
        key = (icdcodes, atc_codes)
        result[key] = entry
    return result

def join_reference_set_and_results(reference_set, observational_results):
    reference_set_map = create_map(reference_set)
    observational_results_map = create_map(observational_results)
    result = []
    for k, v in observational_results_map.items():
        if k not in reference_set_map:
            continue
        entry = copy.deepcopy(reference_set_map[k])
        same_direction = set(entry['atc_codes'][0]) == set(v['atc_codes'][0])
        if same_direction:
            effect_modifier = 1
        else:
            effect_modifier = -1
            
        entry['results'] = {}
        for model_category, model_types in v["results"].items():
            entry['results'][model_category] = {}

            for model_type, model_results in model_types.items():
                entry['results'][model_category][model_type] = {
                    'p': model_results['p'],
                    'se': model_results['se'],
                    'coef': effect_modifier * model_results['coef'],
                }
                
        result.append(entry)
    return result