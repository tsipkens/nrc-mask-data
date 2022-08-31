
import json

def json_read(fn):
    with open(fn, 'r') as f:
        data = json.load(f)
    return data

s = json_read('sd22/data.json')
