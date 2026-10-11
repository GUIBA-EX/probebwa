#!/usr/bin/env python3
"""Summarize observable candidate-to-final changes in a Stampy trace.

This attributes changes to the combined pairing/realignment stages. It does
not infer that an absent candidate was lost specifically by the hash index.
"""
import argparse
import collections
import json
from pathlib import Path


def summarize(path):
    latest = {}
    counts = collections.Counter()
    examples = []
    with open(path) as stream:
        for line in stream:
            event = json.loads(line)
            counts[event['stage']] += 1
            if event['stage'] == 'candidates':
                latest[event['sequence']] = event
            elif event['stage'] == 'pair_report':
                preliminary = latest.get(event['sequence'])
                if preliminary is None or not preliminary['candidates']:
                    counts['final_without_initial_candidate'] += 1
                    continue
                first = preliminary['candidates'][0]
                final = event['candidate']
                changed = []
                for name, index in [('score', 0), ('read_alignment', 1),
                                    ('reference_alignment', 2), ('signed_locus', 3)]:
                    if first[index] != final[index]:
                        counts[name + '_changed'] += 1
                        changed.append(name)
                if changed and len(examples) < 20:
                    examples.append(dict(read_id=event['read_id'], changed=changed,
                                         preliminary=first[:4], final=final[:4]))
    return dict(counts=dict(counts), examples=examples)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trace', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.trace)
    args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['counts'], indent=2))
