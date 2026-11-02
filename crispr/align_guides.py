#!/usr/bin/env python3
"""
Map 20 nt sgRNA spacers to the genome by exact match, on both strands, and
report every hit with its PAM and the gene it falls in.

Exact matching rather than a mismatch-tolerant aligner: a guide is expected to
match its intended target perfectly, and at 20 nt an exact search over both
strands is both complete and unambiguous. All hits are kept, so a guide that
matches more than one locus is visible instead of being assigned arbitrarily.

The Cas9 PAM is checked but never used to filter: a 20 nt match with no NGG
next to it is a coincidental genomic match, not a target site, and saying so is
more useful than dropping it silently.

Usage: align_guides.py spacers.tsv genome.fa.gz genes.gtf.gz OUT_PREFIX
"""
import gzip
import sys
from collections import defaultdict

import ahocorasick
import pandas as pd

SPACERS, GENOME, GTF, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

CHUNK = 8 << 20          # 8 MB of sequence at a time
PAD = 32                 # overlap between chunks; also room for PAM context
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


# ------------------------------------------------------------ patterns -------
sp = pd.read_csv(SPACERS, sep="\t")
auto = ahocorasick.Automaton()
pat = defaultdict(list)                      # sequence -> [(name, strand), ...]
for name, spacer in zip(sp["name"], sp.spacer):
    pat[spacer].append((name, "+"))
    pat[rc(spacer)].append((name, "-"))
for seq, hits in pat.items():
    auto.add_word(seq, (seq, hits))
auto.make_automaton()
print(f"{len(sp)} spacers -> {len(pat)} search patterns (both strands)",
      file=sys.stderr)


# -------------------------------------------------------------- scan ---------
def scan(path):
    """Yield (chrom, end_offset, matched_seq, hits, local_seq, local_off)."""
    chrom, buf, base = None, "", 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if chrom:
                    yield from flush(chrom, buf, base, final=True)
                chrom = line[1:].split()[0]
                buf, base = "", 0
                print(f"  scanning {chrom}", file=sys.stderr)
                continue
            buf += line.strip().upper()
            if len(buf) >= CHUNK + PAD:
                yield from flush(chrom, buf, base)
                keep = PAD
                base += len(buf) - keep
                buf = buf[-keep:]
    if chrom:
        yield from flush(chrom, buf, base, final=True)


def flush(chrom, buf, base, final=False):
    for end, (seq, hits) in auto.iter(buf):
        start = end - len(seq) + 1
        # a hit inside the retained overlap was already reported last chunk
        if not final and start < 0:
            continue
        yield chrom, base + start, seq, hits, buf, start


records = []
seen = set()
for chrom, gpos, seq, hits, buf, local in scan(GENOME):
    key = (chrom, gpos, seq)
    if key in seen:
        continue
    seen.add(key)
    for name, strand in hits:
        if strand == "+":
            pam = buf[local + 20:local + 23]
        else:
            pam = rc(buf[max(local - 3, 0):local])
        records.append({"name": name, "chrom": chrom, "start": gpos + 1,
                        "end": gpos + 20, "strand": strand, "spacer": seq
                        if strand == "+" else rc(seq), "PAM": pam,
                        "PAM_ok": int(len(pam) == 3 and pam[1:] == "GG")})

hits = pd.DataFrame(records)
print(f"{len(hits)} genomic matches for {hits['name'].nunique() if len(hits) else 0} "
      f"of {len(sp)} spacers", file=sys.stderr)


# --------------------------------------------------------- annotation --------
genes = []
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene":
            continue
        a = f[8]

        def attr(k):
            i = a.find(k + ' "')
            if i < 0:
                return ""
            j = a.index('"', i + len(k) + 2)
            return a[i + len(k) + 2:j]

        genes.append((f[0], int(f[3]), int(f[4]), f[6], attr("gene_id"),
                      attr("gene_name"), attr("gene_biotype")))
gtf = pd.DataFrame(genes, columns=["chrom", "start", "end", "strand",
                                   "gene_id", "gene_name", "gene_biotype"])
print(f"{len(gtf)} genes from the GTF", file=sys.stderr)

by_chrom = {c: d.sort_values("start").reset_index(drop=True)
            for c, d in gtf.groupby("chrom")}


def overlapping(chrom, s, e):
    d = by_chrom.get(chrom)
    if d is None:
        return []
    m = d[(d.start <= e) & (d.end >= s)]
    return list(zip(m.gene_id, m.gene_name, m.gene_biotype, m.strand,
                    m.start, m.end))


rows = []
for r in hits.itertuples(index=False):
    ov = overlapping(r.chrom, r.start, r.end)
    if not ov:
        rows.append({**r._asdict(), "gene_id": "", "gene_name": "",
                     "gene_biotype": "intergenic", "tss_dist": ""})
        continue
    for gid, gname, gtype, gstrand, gs, ge in ov:
        tss = gs if gstrand == "+" else ge
        rows.append({**r._asdict(), "gene_id": gid, "gene_name": gname,
                     "gene_biotype": gtype,
                     "tss_dist": int((r.start + r.end) // 2 - tss)
                     * (1 if gstrand == "+" else -1)})

ann = pd.DataFrame(rows)
ann.to_csv(f"{OUT}_hits.tsv", sep="\t", index=False)

# one row per guide, collapsing its hits
def collapse(g):
    coding = g[g.gene_biotype == "protein_coding"]
    pick = coding if len(coding) else g
    names = sorted({x for x in pick.gene_name if x})
    ids = sorted({x for x in pick.gene_id if x})
    return pd.Series({
        "n_hits": len(g), "n_PAM_ok": int(g.PAM_ok.sum()),
        "gene_id": ";".join(ids), "gene_name": ";".join(names),
        "gene_biotype": ";".join(sorted({x for x in pick.gene_biotype if x})),
        "loci": ";".join(f"{r.chrom}:{r.start}-{r.end}({r.strand})"
                         for r in g.itertuples()),
        "PAM": ";".join(sorted(set(g.PAM))),
    })


per_guide = ann.groupby("name").apply(collapse, include_groups=False).reset_index()
out = sp[["name", "spacer", "extra_G"]].merge(per_guide, on="name", how="left")
out["n_hits"] = out.n_hits.fillna(0).astype(int)
out.to_csv(f"{OUT}_per_guide.tsv", sep="\t", index=False)

print(f"\nguides with no genomic match: {(out.n_hits == 0).sum()}", file=sys.stderr)
print(f"guides matching exactly one locus: {(out.n_hits == 1).sum()}", file=sys.stderr)
print(f"guides matching >1 locus: {(out.n_hits > 1).sum()}", file=sys.stderr)
print(f"wrote {OUT}_hits.tsv and {OUT}_per_guide.tsv", file=sys.stderr)
