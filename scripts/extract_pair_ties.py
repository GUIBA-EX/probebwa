#!/usr/bin/env python3
"""Extract equal-best pairs in original Stampy tuple order from trace_stampy.py."""
import argparse
import gzip
from pathlib import Path
import json

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('trace', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
lines = [
    '# Original Stampy realign_pairs() trace; equal-best pairs in returned order.\n',
    '# case, joint score, mate1(score/read/ref/signed locus), mate2(score/read/ref/signed locus).\n',
]
cases = 0
with args.trace.open() as source:
    for line in source:
        event = json.loads(line)
        if event['stage'] != 'pair_candidates' or not event['tried']:
            continue
        candidates = event['candidates']
        best = min(candidate[0] for candidate in candidates)
        tied = [candidate for candidate in candidates if candidate[0] == best]
        if len(tied) < 2:
            continue
        cases += 1
        for candidate in tied:
            fields = [cases, candidate[0]] + candidate[1][:4] + candidate[2][:4]
            lines.append('\t'.join(map(str, fields)) + '\n')
data = ''.join(lines).encode()
with args.output.open('wb') as target:
    if args.output.suffix == '.gz':
        with gzip.GzipFile(filename='', fileobj=target, mode='wb', mtime=0) as compressed:
            compressed.write(data)
    else:
        target.write(data)
print('{} cases, {} pairs'.format(cases, len(lines) - 2))
