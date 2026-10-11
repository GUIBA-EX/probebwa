#!/usr/bin/env python3
"""Fail unless SAM records and non-program headers match the reference.

Only RNEXT '=' versus the explicit reference name and optional-tag ordering
are normalized. @PG and @CO are reported separately as run metadata. No
coordinate tolerance, MAPQ tolerance, or exclusion of unmapped reads is used.
"""
import argparse
import collections
import json
from pathlib import Path

FIELDS = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
          'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'TAGS']


def load(path):
    records = collections.defaultdict(list)
    order, headers, metadata = [], [], []
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith('@'):
                (metadata if line[:3] in ('@PG', '@CO') else headers).append(line.rstrip('\n'))
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 11:
                raise ValueError('%s:%d: incomplete SAM record' % (path, number))
            flag = int(cols[1])
            key = (cols[0], flag & 0xC0)
            order.append(key)
            if cols[6] == '=':
                cols[6] = cols[2]
            records[key].append(cols[:11] + [tuple(sorted(cols[11:]))])
    return records, order, headers, metadata


def compare(reference, actual):
    left, lo, lh, lm = load(reference)
    right, ro, rh, rm = load(actual)
    mismatches = collections.Counter({name: 0 for name in FIELDS})
    equal_core, equal_all, compared, missing, extra = 0, 0, 0, 0, 0
    examples = []
    for key in sorted(left.keys() | right.keys()):
        a, b = left[key], right[key]
        missing += max(0, len(a) - len(b))
        extra += max(0, len(b) - len(a))
        for index, (x, y) in enumerate(zip(a, b)):
            compared += 1
            equal_core += x[:11] == y[:11]
            equal_all += x == y
            fields = [name for name, xv, yv in zip(FIELDS, x, y) if xv != yv]
            mismatches.update(fields)
            if fields and len(examples) < 20:
                examples.append(dict(read=key, occurrence=index, fields=fields,
                                     reference=x, actual=y))
    passed = (not missing and not extra and compared == equal_all and lo == ro and lh == rh)
    return dict(passed=passed, reference_records=len(lo), actual_records=len(ro),
                compared=compared, equal_core_records=equal_core, equal_all_records=equal_all,
                missing_records=missing, extra_records=extra,
                record_order_equal=lo == ro, headers_equal=lh == rh,
                field_mismatches=dict(mismatches), examples=examples,
                reference_headers=lh, actual_headers=rh,
                reference_run_metadata=lm, actual_run_metadata=rm)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('actual', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.reference, args.actual)
    args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ('examples', 'reference_headers', 'actual_headers',
                      'reference_run_metadata', 'actual_run_metadata')}, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
