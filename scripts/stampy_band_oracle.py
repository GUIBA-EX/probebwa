#!/usr/bin/env python2
"""Call original exported C coarse aligner to regenerate scalar-port fixtures.
Usage: python2 stampy_band_oracle.py /path/to/stampy output.tsv
"""
from __future__ import print_function
import ctypes as c
import os
import random
import sys

native = c.CDLL(os.path.join(sys.argv[1], 'maptools.so'))
fn = native.traceback_sse2
fn.argtypes = [c.c_char_p, c.c_char_p, c.POINTER(c.c_int), c.c_char_p, c.c_char_p,
               c.c_int, c.c_int, c.c_int, c.c_int, c.c_int, c.POINTER(c.c_int), c.c_int]
fn.restype = c.c_int
rng = random.Random(20261013)
with open(sys.argv[2], 'w') as out:
    out.write('# Stampy 1.0.32 traceback_sse2: id read quals ref open extend prior score read_aligned ref_aligned offset\n')
    for i in range(600):
        n = rng.choice([10, 19, 35, 50, 100, 151])
        seq = ''.join(rng.choice('ACGTN') for _ in range(n))
        padding = rng.randrange(16)
        ref = [rng.choice('ACGT') for _ in range(padding)] + list(seq)
        for _ in range(rng.randrange(12)):
            k = rng.randrange(len(ref))
            op = rng.randrange(3)
            if op == 0: ref[k] = rng.choice('ACGTN')
            elif op == 1: ref.insert(k, rng.choice('ACGT'))
            else: del ref[k]
        ref = ''.join(ref + ['G']*(n+15))[:n+15]
        qualities = [rng.choice([0, 3, 10, 27, 30, 40, 60]) for _ in seq]
        go, ge, prior = rng.choice([20, 40, 42]), rng.choice([1, 3, 5]), rng.choice([3, 5, 7])
        qal = (c.c_int*n)(*qualities)
        ar, ag = c.create_string_buffer(2*n+64), c.create_string_buffer(2*n+64)
        offset = c.c_int()
        score = fn(ref, seq, qal, ag, ar, len(ref), n, go, ge, prior, c.byref(offset), 0x3c00)
        row = [i, seq, ','.join(map(str, qualities)), ref, go, ge, prior, score, ar.value, ag.value, offset.value]
        out.write('\t'.join(map(str, row))+'\n')
print('Wrote 600 native coarse cases')
