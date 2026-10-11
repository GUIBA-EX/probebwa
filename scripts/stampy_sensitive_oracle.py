#!/usr/bin/env python2
"""Emit unbanded native Viterbi fixtures; does not claim banded/BAQ parity.

python2 stampy_sensitive_oracle.py /path/to/stampy output.tsv [trace.jsonl]
The optional trace is replayed with band=-1, mode=2 and no score bailout.
All original reference windows and scoring parameters are retained.
"""
from __future__ import print_function
import json
import argparse
import os
import random
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('stampy')
parser.add_argument('output')
parser.add_argument('trace', nargs='?')
parser.add_argument('--band', type=int, default=-1, choices=[-1, 8, 15])
args = parser.parse_args()
sys.path.insert(0, os.path.abspath(args.stampy))
import maptools

rng = random.Random(20261012)
cases = []
for i in range(300):
    n = rng.choice([10, 19, 35, 50, 100])
    seq = ''.join(rng.choice('ACGT') for _ in range(n))
    ref = list(seq)
    for j in range(rng.randrange(1, 7)):
        k = rng.randrange(len(ref))
        op = rng.randrange(3)
        if op == 0: ref[k] = rng.choice('ACGTN')
        elif op == 1: ref.insert(k, rng.choice('ACGT'))
        elif len(ref) > 1: del ref[k]
    if i % 5 == 0:
        seq = seq[:n//2] + 'N' + seq[n//2+1:]
    ref = 'GATC' + ''.join(ref) + 'CCAGT'
    qualities = [rng.choice([0, 3, 10, 27, 30, 40, 60]) for _ in seq]
    cases.append((seq, qualities, ref, [5, rng.choice([20, 40, 42]), 3,
                                         rng.choice([0, 4, 10]), rng.choice([5, 10000])]))
# Explicit homopolymer/ambiguous placement and low-quality terminal insertion.
for seq, ref in [('AAAAAAA', 'AAAAAAAAAA'), ('ACGTAAAAACGT', 'ACGTAAAACGT'),
                 ('ACGTAAAACGT', 'ACGTAAAAACGT'), ('NNNNNN', 'ACGTACGT'),
                 ('TTACGTACGT', 'ACGTACGT')]:
    cases.append((seq, [30]*len(seq), ref, [5, 42, 3, 0, 10000]))
if args.trace:
    with open(args.trace) as trace:
        for line in trace:
            e = json.loads(line)
            if e['stage'] == 'sensitive':
                cases.append((e['sequence'].encode('ascii'), e['qualities'],
                              e['reference'].encode('ascii'), e['parameters'][:5]))
with open(args.output, 'w') as out:
    out.write('# Stampy 1.0.32 native Viterbi: id read quals ref ins open extend mean sd score primary_read primary_ref offset alt_read alt_ref band\n')
    for i, (seq, qualities, ref, p) in enumerate(cases):
        # Cython converts these native arguments to C integers.
        p = [int(v) for v in p]
        band = -1 if args.band < 0 else max(args.band, len(ref) - len(seq))
        if band > 15:
            raise ValueError('This oracle only covers the fixed core, not findband')
        result = maptools.alignsensitive(seq, qualities, ref, *(p + [999999, band, 2]))
        row = [i, seq, ','.join(map(str, qualities)), ref] + p + list(result[:6]) + [band]
        out.write('\t'.join(map(str, row)) + '\n')
print('Wrote', len(cases), 'native cases with requested band', args.band)
