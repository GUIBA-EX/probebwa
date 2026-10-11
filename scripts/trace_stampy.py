#!/usr/bin/env python2
"""Observe original Stampy calls without substituting their results.

python2 trace_stampy.py --trace trace.jsonl /path/to/stampy.py -- [Stampy args]
Coordinates in events are Stampy's signed genome coordinates. Alignment
strings and reference slices are retained to avoid silently mixing coordinate
systems or pretending the two implementations' score units are equivalent.
Run with one Stampy worker so event order follows read order.
"""
from __future__ import print_function
import argparse
import json
import os
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--trace', required=True)
parser.add_argument('--pair-candidates', action='store_true', help='record native slow-path pair tuple ordering')
parser.add_argument('--candidates-only', action='store_true', help='omit sensitive alignment events')
parser.add_argument('--tie-counter', action='store_true', help='read the native ELF placement counter (requires nm symbols)')
parser.add_argument('stampy')
parser.add_argument('arguments', nargs=argparse.REMAINDER)
args = parser.parse_args()
sys.path.insert(0, os.path.dirname(os.path.abspath(args.stampy)))
import maptools
from Stampy import singleendedmapper

tie_counter = None
if args.tie_counter:
    import ctypes
    import subprocess
    symbols = {}
    for line in subprocess.check_output(['nm', maptools.__file__]).splitlines():
        fields = line.split()
        if len(fields) == 3:
            symbols[fields[2]] = int(fields[0], 16)
    counter_name = next(name for name in symbols if name.startswith('pseudo_random.'))
    native = ctypes.CDLL(maptools.__file__)
    address = ctypes.cast(native.unique_alignments_ptr, ctypes.c_void_p).value
    address += symbols[counter_name] - symbols['unique_alignments_ptr']
    # Read only: the counter remains owned and advanced by original Stampy.
    tie_counter = ctypes.c_int.from_address(address)

stream = open(args.trace, 'w')
event = [0]


def emit(stage, payload):
    event[0] += 1
    payload.update(stage=stage, event=event[0])
    stream.write(json.dumps(payload, separators=(',', ':')) + '\n')


locate_original = maptools.locatealignandscore
align_original = maptools.alignsensitive
report_original = singleendedmapper.PairedEndMapper.report
pairs_original = singleendedmapper.PairedEndMapper.realign_pairs


def locate(*a, **kw):
    before = tie_counter.value if tie_counter is not None else None
    result = locate_original(*a, **kw)
    payload = dict(sequence=a[0], repeat_mask=result[0], candidates=result[1], posterior=result[2])
    if tie_counter is not None:
        payload.update(tie_pass_before=before, tie_pass_after=tie_counter.value)
    emit('candidates', payload)
    return result


def align(*a, **kw):
    result = align_original(*a, **kw)
    emit('sensitive', dict(sequence=a[0], qualities=a[1], reference=a[2],
                          parameters=list(a[3:]), result=result))
    return result


def report(self, *a, **kw):
    result = report_original(self, *a, **kw)
    emit('pair_report', dict(read_id=a[2], sequence=a[0],
                            evidence=list(a[5:12]), candidate=a[4],
                            contig=result[6], position=result[7], strand=result[4],
                            read_alignment=result[2], reference_alignment=result[3],
                            single_phred=result[11], pair_phred=result[12],
                            proper_pair=result[27]))
    return result


def realign_pairs(self, *a, **kw):
    result = pairs_original(self, *a, **kw)
    emit('pair_candidates', dict(candidates=result[0], posterior=result[1], tried=result[2]))
    return result


maptools.locatealignandscore = locate
if not args.candidates_only:
    maptools.alignsensitive = align
singleendedmapper.PairedEndMapper.report = report
if args.pair_candidates:
    singleendedmapper.PairedEndMapper.realign_pairs = realign_pairs
native_args = args.arguments
if native_args[:1] == ['--']:
    native_args = native_args[1:]
sys.argv = [args.stampy] + native_args
try:
    with open(args.stampy, 'rb') as source:
        code = compile(source.read(), args.stampy, 'exec', 0, True)
    exec(code, {'__name__': '__main__', '__file__': args.stampy})
finally:
    stream.close()
