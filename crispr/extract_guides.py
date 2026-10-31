#!/usr/bin/env python3
"""
Pull 20 bp spacers out of oligo sequences in a spreadsheet.

Layout-agnostic on purpose: every cell of every sheet is scanned for the U6
sgRNA context, so the script does not care which column the oligos sit in.

    ...TATATATCTTGTGGAAAGGACGAAACACC [ G? ] NNNNNNNNNNNNNNNNNNNN GTTT...

The vector usually carries a transcription-start G between the 5' flank and the
spacer, so the captured insert is 20 or 21 nt. The script does not assume which:
it records the raw insert and its length, trims a leading G only when the insert
is 21 nt, and prints the length distribution so the call is made on the data.

Usage: extract_guides.py INPUT.xlsx OUT.tsv
"""
import re
import sys
from collections import Counter

import pandas as pd

IN, OUT = sys.argv[1], sys.argv[2]

FLANK5 = "TATATATCTTGTGGAAAGGACGAAACACC"
FLANK3 = "GTTT"
# capture everything between the flanks, 18-22 nt, so a wrong-length oligo is
# surfaced rather than silently skipped
PAT = re.compile(FLANK5 + r"([ACGT]{18,22}?)" + FLANK3, re.I)
DNA_ONLY = re.compile(r"^[ACGTNacgtn]+$")

sheets = pd.read_excel(IN, sheet_name=None, header=None, dtype=str)
rows, scanned, dna_cells = [], 0, 0
for sheet, df in sheets.items():
    for r in range(df.shape[0]):
        for c in range(df.shape[1]):
            v = df.iat[r, c]
            if not isinstance(v, str):
                continue
            scanned += 1
            s = v.strip().replace(" ", "").replace("\n", "")
            if DNA_ONLY.match(s):
                dna_cells += 1
            m = PAT.search(s)
            if not m:
                continue
            insert = m.group(1).upper()
            spacer = insert[1:] if len(insert) == 21 and insert[0] == "G" else insert
            # the label is usually the first non-sequence cell on the same row
            label = ""
            for cc in range(df.shape[1]):
                o = df.iat[r, cc]
                if isinstance(o, str) and o.strip() and not DNA_ONLY.match(
                        o.strip().replace(" ", "")):
                    label = o.strip()
                    break
            rows.append({"sheet": sheet, "row": r + 1, "col": c + 1,
                         "label": label, "insert": insert,
                         "insert_len": len(insert), "spacer": spacer,
                         "spacer_len": len(spacer), "oligo": s})

if not rows:
    sys.exit(f"no oligo matched the expected context in {IN}\n"
             f"  cells scanned: {scanned}, pure-DNA cells: {dna_cells}\n"
             f"  expected: {FLANK5} [G?] <20 nt> {FLANK3}...")

out = pd.DataFrame(rows)
print("insert length distribution:", dict(Counter(out.insert_len)), file=sys.stderr)
print("spacer length distribution:", dict(Counter(out.spacer_len)), file=sys.stderr)
dups = out.spacer.duplicated(keep=False)
print(f"{len(out)} oligos matched; {out.spacer.nunique()} unique spacers; "
      f"{int(dups.sum())} rows share a spacer", file=sys.stderr)
out.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}", file=sys.stderr)
