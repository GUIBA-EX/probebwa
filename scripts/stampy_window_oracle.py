#!/usr/bin/env python2
"""Replay observed sensitive calls through original findband and emit fixtures.
Usage: python2 stampy_window_oracle.py /path/to/stampy trace.jsonl output.tsv
"""
from __future__ import print_function
import ctypes as c
import json
import os
import sys
native = c.CDLL(os.path.join(sys.argv[1], 'maptools.so'))
fn = native.findband
fn.argtypes = [c.c_char_p, c.c_char_p, c.c_int, c.c_int, c.c_char_p] + [c.c_int]*6 + [c.POINTER(c.c_int)]*3
fn.restype = None
with open(sys.argv[3], 'w') as out:
    out.write('# Stampy native findband: event read quals reference open extend prior band start length width\n')
    for line in open(sys.argv[2]):
        e = json.loads(line)
        if e['stage'] != 'sensitive': continue
        seq, ref = e['sequence'].encode('ascii'), e['reference'].encode('ascii')
        p = [int(v) for v in e['parameters']]
        band = max(p[6], len(ref)-len(seq))
        if band <= 15: continue
        start, length, width = c.c_int(), c.c_int(), c.c_int()
        fn(seq, ref, len(seq), len(ref), ''.join(chr(q+33) for q in e['qualities']),
           p[0], p[1], p[2], band, 150, 25, c.byref(start), c.byref(length), c.byref(width))
        row = [e['event'], seq, ','.join(map(str,e['qualities'])), ref, p[1], p[2], p[0], band, start.value, length.value, width.value]
        out.write('\t'.join(map(str,row))+'\n')
