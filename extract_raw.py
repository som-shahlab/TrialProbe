import collections
import itertools
import json
import lxml.etree
import math
import os

import converters

medra_converter = converters.MeDRAConverter('/labs/shahlab/projects/ethanid/2019AB/')
rxnorm_converter = converters.RXNormConverter('cache.db')

num_experimental_studies = 0
num_possible_entries = 0
num_mappable_entries = 0

with open('raw_entries.txt', 'w') as found_entries:
    with open('nct_with_events.txt') as f:
        for line in f:
            study = line.split(':')[0].replace('nct2/', 'nct/')
            try:
                with open(study) as f:
                    document = lxml.etree.parse(f)

                study_type = document.find('study_type').text

                status = document.find('overall_status').text

                result = document.find('clinical_results')

                if result is None:
                    continue

                enrollment = document.find('enrollment')

                if enrollment is None:
                    continue

                enrollment =  int(enrollment.text)

                events_xml = result.find('reported_events')

                if events_xml is None:
                    continue

                if study_type != 'Interventional':
                    # print('Ignore study that is not interventional?')
                    # print(study)
                    # print(study_type)
                    continue
                if status != 'Completed':
                    # print('Ignore study that is not complete?')
                    # print(study)
                    # print(status)
                    continue

                num_experimental_studies += 1

                if enrollment < 200:
                    continue


                event_items = events_xml.findall('.//event')
                title_to_events = collections.defaultdict(list)

                for event_item in event_items:
                    count_map = {}

                    if event_item.find('sub_title') is None:
                        continue

                    groups = event_item.findall('counts')
                    for group in groups:
                        group_id = group.attrib['group_id']

                        subjects_affected = group.attrib.get('subjects_affected', '')
                        if subjects_affected == '':
                            continue

                        subjects_at_risk = group.attrib.get('subjects_at_risk', '')
                        if subjects_at_risk == '' or subjects_at_risk == '0':
                            continue

                        if int(subjects_at_risk) < 100:
                            continue

                        count_map[group_id] = [
                            int(subjects_affected),
                            int(subjects_at_risk),
                        ]

                    title_to_events[event_item.find('sub_title').text.lower()].append(count_map)

                title_to_event = {}

                for title, events in title_to_events.items():
                    events.sort(key=lambda a: -sum(v[0] for v in a.values()))
                    resulting_map = events[0]

                    for event in events[1:]:
                        for k in event.keys():
                            if k not in resulting_map:
                                resulting_map[k] = event[k]
                            else:
                                if resulting_map[k][1] != event[k][1]:
                                    print('Not aligned', study, events)
                                    continue
                                resulting_map[k][0] += event[k][0]    

                    title_to_event[title] = resulting_map

                groups = events_xml.find('group_list').findall('group')

                num_possible_entries += math.comb(len(groups), 2) * len(title_to_event)

                if len(title_to_event) == 0:
                    continue
                
                group_map = {}
                next_group_index = 0
                drug_list = []
                rxnorm_list = []
                for group in groups:
                    group_id = group.attrib['group_id']
                    group_title_node = group.find('title')
                    if group_title_node is None:
                        continue

                    group_title = group_title_node.text
                    if ' + ' in group_title:
                        continue
                    atc_codes = rxnorm_converter.get_atc_codes(group_title)

                    if atc_codes is None:
                        continue

                    group_map[group_id] = (atc_codes, group_title, next_group_index)
                    next_group_index += 1
                    drug_list.append(group_title)
                    rxnorm_list.append(atc_codes)
                
                for title, event in title_to_event.items():
                    for group_ids in itertools.combinations(list(event.keys()), 2):
                        table = [[0, 0], [0, 0]]
                        drugs = []
                        atc_codes = []
                        ratios = []
                        counts = []
                        if any(group_id not in group_map for group_id in group_ids):
                            continue

                        for i, group_id in enumerate(group_ids):
                            atc, group_title, group_index = group_map[group_id]
                            subjects_affected = event[group_id][0]
                            subjects_at_risk = event[group_id][1]

                            table[i][0] = subjects_affected
                            table[i][1] = subjects_at_risk - subjects_affected
                            ratios.append(subjects_affected / subjects_at_risk)
                            counts.append(subjects_at_risk)
                            drugs.append(group_title)
                            atc_codes.append(sorted(list(atc)))

                        if drugs[0] in drugs[1] or drugs[1] in drugs[0]:
                            continue

                        if atc_codes[0] == atc_codes[1]:
                            continue

                        icd10codes = medra_converter.convert_to_icd10(title)

                        if icd10codes is None:
                            print('Could not find ICD10 codes for', title, study)
                            continue

                        num_mappable_entries += 1

                        study_data = {
                            'study': study,
                            'side_effect_name': title,
                            'icd10codes': list(icd10codes),
                            'table': table,
                            'drugs': drugs,
                            'atc_codes': atc_codes,

                        }
                        found_entries.write(json.dumps(study_data) + '\n')
                   
            except Exception as e:
                print('Got error working with ', study, e)
                raise e
                continue


print(num_experimental_studies)
print(num_possible_entries)
print(num_mappable_entries)
