import json
import simplejson

with open('results.json') as f:
    with open('clean_results.json', 'w') as o:
        for line in f:
            data = json.loads(line)
            o.write(simplejson.dumps(data, ignore_nan=True) + '\n')