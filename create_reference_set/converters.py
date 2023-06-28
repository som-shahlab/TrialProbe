import collections
import json
import os
import sqlite3
import urllib.parse
import urllib.request
import xml.dom.minidom

class MeDRAConverter:
    def __init__(self, mrconso_path):
        bad_dra_duplicate = set()
        dra_to_concept = {}
        concept_to_icd10 = collections.defaultdict(set)
        concept2_to_icd10 = collections.defaultdict(set)
        vocabs = set()

        bad_dra_no_mapping = set()

        valid_icd10_codes = set()

        with open(mrconso_path) as f:
            for line in f:
                parts = line.split('|')
                cui = parts[0]
                lang = parts[1]
                vocab = parts[11]
                term = parts[14].lower()
                code = parts[13]

                if lang != 'ENG':
                    continue

                if vocab == 'MDR':
                    previous_cui = dra_to_concept.get(term)
                    if previous_cui is not None and previous_cui != cui:
                        bad_dra_duplicate.add(term)
                    dra_to_concept[term] = cui
                elif vocab == 'ICD10CM':
                    concept_to_icd10[cui].add(code)
                    valid_icd10_codes.add(code)

                if vocab in ('ICD10', 'ICD10CM'):
                    vocabs.add(vocab)
                    concept2_to_icd10[cui].add(code)

        print(vocabs)

        self.dra_to_icd10 = {}

        for dra, cui in dra_to_concept.items():
            if dra in bad_dra_duplicate:
                continue
            icd10codes = concept2_to_icd10[cui] & valid_icd10_codes
            # if concept_to_icd10[cui] != concept2_to_icd10[cui]:
            #     print('Lol', dra, 
            #         concept_to_icd10[cui], 
            #         concept2_to_icd10[cui],
            #         concept2_to_icd10[cui] & valid_icd10_codes)
            if len(icd10codes) == 0:
                bad_dra_no_mapping.add(dra)
            else:
                self.dra_to_icd10[dra] = icd10codes

    def convert_to_icd10(self, dra_text):
        return self.dra_to_icd10.get(dra_text)

class RXNormConverter:
    def __init__(self, cache_db_path):
        self.cache = sqlite3.connect(cache_db_path)
        self.cache.executescript("""
            CREATE TABLE IF NOT EXISTS cached_queries (
                query TEXT PRIMARY KEY,
                result TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS cached_in_queries (
                query TEXT PRIMARY KEY,
                result TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS cached_atc_queries (
                query TEXT PRIMARY KEY,
                result TEXT NOT NULL
            );
        """)


    def run_query(self, a):
        cursor = self.cache.cursor()

        cursor.execute("SELECT result from cached_queries where query = ?", (a,))

        cached = list(cursor)
        if cached:
            return cached[0][0]
        else:
            url = 'https://rxnav.nlm.nih.gov/REST/approximateTerm?' + urllib.parse.urlencode({
                'term': a,
                'option': '1'
            })

            with urllib.request.urlopen(url) as f:
                result = f.read().decode('utf-8')
                cursor.execute("INSERT INTO cached_queries (query, result) VALUES (?, ?)", (a, result))

            cursor.close()
            self.cache.commit()

            return result

    def run_in_query(self, a):
        cursor = self.cache.cursor()

        cursor.execute("SELECT result from cached_in_queries where query = ?", (a,))

        cached = list(cursor)
        if cached:
            return cached[0][0]
        else:
            url = f'https://rxnav.nlm.nih.gov/REST/rxcui/{a}/related?tty=IN'
            with urllib.request.urlopen(url) as f:
                result = f.read().decode('utf-8')
                cursor.execute("INSERT INTO cached_in_queries (query, result) VALUES (?, ?)", (a, result))

            cursor.close()
            self.cache.commit()

            return result

    def run_atc_query(self, a):
        cursor = self.cache.cursor()

        cursor.execute("SELECT result from cached_atc_queries where query = ?", (a,))

        cached = list(cursor)
        if cached:
            return cached[0][0]
        else:
            url = f'https://rxnav.nlm.nih.gov/REST/rxcui/{a}/allProperties.json?prop=codes'
            with urllib.request.urlopen(url) as f:
                result = f.read().decode('utf-8')
                cursor.execute("INSERT INTO cached_atc_queries (query, result) VALUES (?, ?)", (a, result))

            cursor.close()
            self.cache.commit()

            return result

    def get_atc(self, a):
        info = json.loads(self.run_atc_query(a))

        atc_codes = set()

        for entry in info["propConceptGroup"]["propConcept"]:
            if entry['propName'] == 'ATC':
                atc_codes.add(entry['propValue'])

        return atc_codes

    def getrxnorm_in(self, a):
        document = xml.dom.minidom.parseString(self.run_in_query(a))
        properties = document.getElementsByTagName('conceptProperties')

        result = []

        for propertyEntry in properties:
            rxcui = propertyEntry.getElementsByTagName('rxcui')

            if len(rxcui) != 1:
                print('Wat2', a, properties.toxml())
                exit(-1)

            result.append(int(rxcui[0].childNodes[0].data))

        return result

    def get_atc_codes(self, a):
        document = xml.dom.minidom.parseString(self.run_query(a))
        
        found = set()
        result = []

        candidates = document.getElementsByTagName('candidate')
        for candidate in candidates:
            cui = int(candidate.getElementsByTagName('rxcui')[0].childNodes[0].data)
            score = int(candidate.getElementsByTagName('score')[0].childNodes[0].data)
            if score <= 50:
                continue

            ingredients = self.getrxnorm_in(cui)

            for ingredient in ingredients:
                if ingredient not in found:
                    result.append(ingredient)
                    found.add(ingredient)

        if len(result) != 1:
            print(f'More or less than one ingredient for: {a}, got {result}')
            return None

        ingredient = result[0]

        atc_codes = self.get_atc(ingredient)

        if len(atc_codes) == 0:
            print(f'Could not find atc code for: {a} and in particular {ingredient}')
            return None

        return atc_codes

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.cache.close()
