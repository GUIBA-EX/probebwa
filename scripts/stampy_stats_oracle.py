#!/usr/bin/env python
"""Regenerate native Stampy statistics fixtures using its Python 2 environment.

Usage: python2 scripts/stampy_stats_oracle.py /path/to/stampy tests/fixtures
Only imports the supplied reference implementation; never substitutes a model.
"""
from __future__ import print_function
import os
import random
import sys

sys.path.insert(0, os.path.abspath(sys.argv[1]))
import maptools

rng = random.Random(20261011)
outdir = sys.argv[2]
if not os.path.isdir(outdir):
    os.makedirs(outdir)

with open(os.path.join(outdir, 'stampy_mutation.tsv'), 'w') as out:
    out.write('# Stampy 1.0.32 native calcmutprobs: id qualities repeat_phred variants miss_phred expectedQ\n')
    index = 0
    for length in [14, 15, 18, 19, 20, 29, 34, 35, 40, 49, 50, 52, 75, 100, 151, 250]:
        for quality in [0, 3, 10, 20, 27, 30, 40, 79, 80, 93]:
            for mode in [-1, 1, 2, 3]:
                q = [quality] * length
                mask = [1000] * length
                result = maptools.calcmutprobs(q, mask, mode)
                out.write('%d\t%s\t%s\t%d\t%d\t%d\n' %
                          (index, ','.join(map(str, q)), ','.join(map(str, mask)), mode, result[0], result[1]))
                index += 1
    for index in range(index, index + 240):
        length = rng.choice([19, 34, 35, 49, 50, 100, 151])
        q = [rng.randrange(0, 94) for _ in range(length)]
        mask = [rng.choice([0, 3, 10, 30, 79, 80, 1000]) for _ in range(length)]
        mode = rng.choice([-1, 1, 2, 3])
        result = maptools.calcmutprobs(q, mask, mode)
        out.write('%d\t%s\t%s\t%d\t%d\t%d\n' %
                  (index, ','.join(map(str, q)), ','.join(map(str, mask)), mode, result[0], result[1]))

with open(os.path.join(outdir, 'stampy_entropy.tsv'), 'w') as out:
    out.write('# Stampy 1.0.32 native logentropy: id read_aligned ref_aligned qualities searchspace min_q score\n')
    for index in range(400):
        length = rng.choice([15, 25, 40, 50, 100, 151])
        read = []
        reference = []
        qualities = []
        for i in range(length):
            r = rng.choice('ACGTN')
            op = rng.randrange(12)
            if op == 0:
                read.append('-'); reference.append(r)
            else:
                read.append(r)
                reference.append('-' if op == 1 else (rng.choice('ACGTN') if op < 5 else r))
                qualities.append(rng.choice([0, 3, 9, 10, 20, 30, 40]))
        read = ''.join(read)
        reference = ''.join(reference)
        space = rng.choice([400, 4641652, 3000000000])
        minq = rng.choice([0, 3, 10])
        score = maptools.logentropy(read, len(read), reference, qualities, space, minq)
        out.write('%d\t%s\t%s\t%s\t%d\t%d\t%d\n' %
                  (index, read, reference, ','.join(map(str, qualities)), space, minq, score))
# Exercise the original bytecode's update methods, rather than a Python
# restatement of the Rust implementation.
import new
from Stampy.singleendedmapper import PairedEndMapper
model = new.instance(PairedEndMapper)
model.autosensing = True
model.sd2 = -1.0
model.as_sd_correction = 0.9707
model.as_maxz = 3.0
model.as_multiprocessing = False
model.as_x = [10.0, 2500.0, 10.0 * (60.0 * 60.0 + 250.0 * 250.0), 0.0, 0.0, 0.0]
model.calc_avgdist_sd()
insert_rng = random.Random(20261011)
observations = [10000.0, -100.0, 380.0, 430.0, 250.0, 251.0]
observations += [insert_rng.gauss(400, 60) for _ in range(194)]
with open(os.path.join(outdir, 'stampy_insert.tsv'), 'w') as out:
    out.write('# Stampy 1.0.32 PairedEndMapper native updates: observed mean sd accepted_pairs\n')
    for observed in observations:
        model.update_pair_stats(observed)
        model.update_pair_stats(observed)
        model.calc_avgdist_sd()
        out.write('%.17g\t%.17g\t%.17g\t%d\n' %
                  (observed, model.averagedistance, model.sd, model.as_x[0] - 10))
print('Wrote native Stampy fixtures to', outdir)
