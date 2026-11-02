#!/usr/bin/env python3
"""
Work out what the spacers that found no genomic match are actually targeting.

Checks them against any extra FASTA (the AAV/reporter plasmid, for instance),
both strands and exact, then characterises whatever is still unexplained so the
leftovers are described rather than just counted.

Usage: check_unmapped.py spacers.tsv raw_hits.tsv extra.fa [extra2.fa ...]
"""
import re
import sys
from collections import Counter

import pandas as pd

SPACERS, HITS = sys.argv[1], sys.argv[2]
EXTRA = sys.argv[3:]
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


sp = pd.read_csv(SPACERS, sep="\t")
hit_names = set(pd.read_csv(HITS, sep="\t")["name"]) if HITS != "-" else set()
un = sp[~sp["name"].isin(hit_names)].copy()
print(f"{len(un)} of {len(sp)} spacers had no genomic match", file=sys.stderr)


def load(path):
    name, seqs, cur = None, {}, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(cur).upper()
            name, cur = line[1:].strip(), []
        else:
            cur.append(line.strip())
    if name:
        seqs[name] = "".join(cur).upper()
    return seqs


explained = {}
for path in EXTRA:
    for rec, seq in load(path).items():
        short = rec.split()[0]
        for nm, s in zip(un["name"], un.spacer):
            if nm in explained:
                continue
            if s in seq:
                explained[nm] = (short, "+", seq.find(s) + 1)
            elif rc(s) in seq:
                explained[nm] = (short, "-", seq.find(rc(s)) + 1)

un["match"] = un["name"].map(lambda n: explained.get(n, ("", "", ""))[0])
un["strand"] = un["name"].map(lambda n: explained.get(n, ("", "", ""))[1])
un["pos"] = un["name"].map(lambda n: explained.get(n, ("", "", ""))[2])

print(f"explained by the extra sequences: {len(explained)}", file=sys.stderr)
if explained:
    print(Counter(v[0] for v in explained.values()), file=sys.stderr)

rest = un[un["match"] == ""]
print(f"\nstill unexplained: {len(rest)}", file=sys.stderr)
if len(rest):
    gc = rest.spacer.map(lambda s: (s.count("G") + s.count("C")) / len(s))
    print(f"  GC content: median {gc.median():.2f}, range "
          f"{gc.min():.2f}-{gc.max():.2f}", file=sys.stderr)
    runs = rest.spacer.map(lambda s: bool(re.search(r"(.)\1{3,}", s)))
    print(f"  with a homopolymer run of 4+: {int(runs.sum())}", file=sys.stderr)
    print(f"  first base: {dict(Counter(rest.spacer.str[0]))}", file=sys.stderr)
    print("  examples:", file=sys.stderr)
    for r in rest.head(8).itertuples():
        print(f"    {r.name:22s} {r.spacer}", file=sys.stderr)

un.to_csv(sys.stdout, sep="\t", index=False)
