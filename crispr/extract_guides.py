#!/usr/bin/env python3
"""
Pull 20 nt sgRNA spacers out of full-length oligos in a spreadsheet.

Construct in this library:

    ...CACTTGCTCGTACGACGTATCTTGTGGAAAGGACGAAACACC [G?] <20 nt spacer>
    GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG...   (modified sgRNA scaffold)

The U6 promoter needs a G at the transcription start, so an oligo whose spacer
does not already begin with G carries an extra one. The insert between the 5'
constant and the scaffold is therefore 20 or 21 nt. Matching is greedy over
{20,21}: a 21 nt match is only possible when the extra G is really there,
because a 20 nt match would leave "TTTAAGAGC..." rather than "GTTTAAGAGC..."
where the scaffold is expected. The spacer is then the 20 nt immediately
preceding the scaffold, which is correct either way.

Usage: extract_guides.py INPUT.xlsx OUT.tsv
"""
import re
import sys
from collections import Counter

import pandas as pd

IN, OUT = sys.argv[1], sys.argv[2]

FLANK5 = "GGAAAGGACGAAACACC"
SCAFFOLD = "GTTTAAGAGCTATGCTGG"
PAT = re.compile(FLANK5 + r"([ACGT]{20,21})" + SCAFFOLD)
DNA_ONLY = re.compile(r"^[ACGTNacgtn]+$")

sheets = pd.read_excel(IN, sheet_name=None, header=None, dtype=str)
rows, unmatched = [], []
for sheet, df in sheets.items():
    for r in range(df.shape[0]):
        cells = [df.iat[r, c] for c in range(df.shape[1])]
        text = [str(v).strip() for v in cells if isinstance(v, str) and str(v).strip()]
        seqs = [t.replace(" ", "").replace("\n", "").upper()
                for t in text if DNA_ONLY.match(t.replace(" ", "").replace("\n", ""))]
        if not seqs:
            continue
        label = next((t for t in text
                      if not DNA_ONLY.match(t.replace(" ", ""))), "")
        for s in seqs:
            m = PAT.search(s)
            if not m:
                unmatched.append((sheet, r + 1, label, s[:60]))
                continue
            insert = m.group(1)
            rows.append({
                "sheet": sheet, "row": r + 1, "name": label,
                "insert": insert, "insert_len": len(insert),
                "extra_G": int(len(insert) == 21),
                "spacer": insert[-20:], "oligo_len": len(s),
            })

if not rows:
    sys.exit(f"no oligo matched {FLANK5} [G?] <20 nt> {SCAFFOLD} in {IN}")

out = pd.DataFrame(rows)
print(f"matched {len(out)} oligos; {len(unmatched)} unmatched", file=sys.stderr)
print("insert length:", dict(Counter(out.insert_len)), file=sys.stderr)
print(f"spacers: {out.spacer.nunique()} unique of {len(out)}", file=sys.stderr)
print("first base of spacer:", dict(Counter(out.spacer.str[0])), file=sys.stderr)
if unmatched:
    for u in unmatched[:5]:
        print("  UNMATCHED", u, file=sys.stderr)
dup = out[out.spacer.duplicated(keep=False)]
if len(dup):
    print(f"NOTE {len(dup)} oligos share a spacer with another oligo",
          file=sys.stderr)
out.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}", file=sys.stderr)
