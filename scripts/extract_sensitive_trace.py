#!/usr/bin/env python3
"""Extract native score/path fixtures from trace_stampy.py JSONL and optional BAQ evidence."""
import argparse
import json
import gzip
import io
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('trace')
parser.add_argument('output')
parser.add_argument('--stride', type=int, default=1)
parser.add_argument('--confidence', action='store_true', help='Keep modes 4/5 and their complete confidence result')
args = parser.parse_args()
if args.stride < 1:
    parser.error('stride must be positive')
with io.StringIO() as out:
    out.write('# Native sensitive trace: event read quals ref prior open extend mean sd max band mode score read_aln ref_aln offset [alt_read alt_ref]\n')
    index = 0
    for line in open(args.trace):
        e = json.loads(line)
        if e['stage'] != 'sensitive' or (args.confidence and e['parameters'][-1] < 4):
            continue
        if index % args.stride == 0:
            fields = [e['event'], e['sequence'], ','.join(map(str, e['qualities'])), e['reference']]
            fields += [int(p) for p in e['parameters']]
            fields += e['result'] if args.confidence else e['result'][:6]
            out.write('\t'.join(map(str, fields)) + '\n')
        index += 1
    data = out.getvalue().encode()
    Path(args.output).write_bytes(gzip.compress(data, mtime=0) if args.output.endswith(".gz") else data)
