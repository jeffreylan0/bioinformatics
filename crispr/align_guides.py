#!/usr/bin/env python3
"""
Map 20 nt sgRNA spacers to the genome and report every hit with its PAM and the
gene it falls in.

The search is on positions 2-20, not the full 20 nt. U6 needs a G at the
transcription start, and this library gets one by SUBSTITUTING the first base
of the spacer rather than prepending a base: at the Rpl10 locus the genome
reads TCTGGCGCGGTGAGTGAGTT where the oligo carries GCTGGCGCGGTGAGTGAGTT. An
exact 20 nt search therefore misses every guide whose target does not already
begin with G - two thirds of this library. Searching the trailing 19 nt and
recovering the real 5' base from the genome finds them, and the substitution is
reported per hit so a genuine G-start can be told from a substituted one.

A 19 nt query is still highly specific: about 0.02 occurrences expected by
chance in a 2.7 Gb genome over both strands. All hits are kept, so a guide that
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
pat = defaultdict(list)          # 19 nt query -> [(name, full spacer, strand)]
for name, spacer in zip(sp["name"], sp.spacer):
    q = spacer[1:]                               # drop the possibly substituted base
    pat[q].append((name, spacer, "+"))
    pat[rc(q)].append((name, spacer, "-"))
for seq, hits in pat.items():
    auto.add_word(seq, (seq, hits))
auto.make_automaton()
print(f"{len(sp)} spacers -> {len(pat)} search patterns of 19 nt (both strands)",
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


# Chunks overlap by PAD, so a match near a boundary is seen twice: once with
# its PAM context truncated by the chunk edge and once complete. Deduplicate on
# the locus and keep the copy that has a full 3 nt PAM, rather than whichever
# chunk happened to report it first.
best = {}
for chrom, gpos, seq, hits_, buf, local in scan(GENOME):
    for name, spacer, strand in hits_:
        if strand == "+":
            # 19 nt query sits at [local, local+19); the real 5' base precedes it
            if local < 1 or local + 19 + 3 > len(buf):
                continue
            proto = buf[local - 1:local + 19]
            pam = buf[local + 19:local + 22]
            start = gpos - 1
        else:
            # rc(query) sits at [local, local+19); the 5' base follows it
            if local + 20 > len(buf) or local < 3:
                continue
            proto = rc(buf[local:local + 20])
            pam = rc(buf[local - 3:local])
            start = gpos
        if len(proto) != 20 or "N" in proto:
            continue
        key = (chrom, start, strand, name)
        prev = best.get(key)
        if prev is not None and len(prev["PAM"]) >= len(pam):
            continue
        best[key] = {"name": name, "chrom": chrom, "start": start + 1,
                     "end": start + 20, "strand": strand,
                     "spacer": spacer, "protospacer": proto,
                     "g_substituted": int(proto[0] != spacer[0]),
                     "PAM": pam,
                     "PAM_ok": int(len(pam) == 3 and pam[1:] == "GG")}

records = list(best.values())
hits = pd.DataFrame(records)
truncated = sum(1 for r in records if len(r["PAM"]) != 3)
if truncated:
    print(f"WARNING {truncated} hits still lack a full PAM context",
          file=sys.stderr)

# The genome scan is the expensive half, so persist it before touching the GTF;
# an annotation failure should never cost a rescan.
hits.to_csv(f"{OUT}_raw_hits.tsv", sep="\t", index=False)
print(f"wrote {OUT}_raw_hits.tsv", file=sys.stderr)

unmapped = sorted(set(sp["name"]) - set(hits["name"]))
pd.Series(unmapped, name="name").to_csv(f"{OUT}_unmapped.txt", index=False,
                                        header=False)
print(f"{len(unmapped)} spacers with no genomic match -> {OUT}_unmapped.txt",
      file=sys.stderr)
print(f"{len(hits)} genomic matches for {hits['name'].nunique() if len(hits) else 0} "
      f"of {len(sp)} spacers", file=sys.stderr)
if len(hits):
    ok = hits[hits.PAM_ok == 1]
    print(f"  with a canonical NGG PAM: {len(ok)} hits, "
          f"{ok['name'].nunique()} spacers", file=sys.stderr)
    print(f"  5' base substituted to G: {int(hits.g_substituted.sum())} hits",
          file=sys.stderr)


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
        "protospacer": ";".join(sorted(set(g.protospacer))),
        "g_substituted": int(g.g_substituted.max()),
    })


per_guide = ann.groupby("name").apply(collapse, include_groups=False).reset_index()
out = sp[["name", "spacer", "extra_G"]].merge(per_guide, on="name", how="left")
out["n_hits"] = out.n_hits.fillna(0).astype(int)
out.to_csv(f"{OUT}_per_guide.tsv", sep="\t", index=False)

print(f"\nguides with no genomic match: {(out.n_hits == 0).sum()}", file=sys.stderr)
print(f"guides matching exactly one locus: {(out.n_hits == 1).sum()}", file=sys.stderr)
print(f"guides matching >1 locus: {(out.n_hits > 1).sum()}", file=sys.stderr)
print(f"wrote {OUT}_hits.tsv and {OUT}_per_guide.tsv", file=sys.stderr)
