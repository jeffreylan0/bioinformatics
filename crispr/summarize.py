#!/usr/bin/env python3
"""
Final gene assignment: one row per guide, one row per gene.

Only hits carrying a canonical NGG PAM are used to assign genes. A 19 nt query
can land somewhere that is not a Cas9 target site; such a hit is real sequence
but not a place the guide can cut, so it must not name a gene. The hits file
keeps those rows regardless, so nothing is hidden.

Guides with no genomic hit are reported too, as controls rather than failures.

Usage: summarize.py annotated_hits.tsv spacers.tsv OUT_PREFIX [gene_table.tsv]
"""
import re
import sys
from collections import Counter

import pandas as pd

HITS, SPACERS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

ann = pd.read_csv(HITS, sep="\t")
sp = pd.read_csv(SPACERS, sep="\t")

ok = ann[ann.PAM_ok == 1].copy()
print(f"{len(ann)} hits -> {len(ok)} with a canonical NGG PAM", file=sys.stderr)

# A guide that cuts between annotated genes still belongs to a target. Resolve
# those to the nearest protein-coding gene and record the distance, so they are
# reported as "nearest" rather than left blank -- a promoter-proximal guide
# often falls just outside the gene body it is meant to silence.
GTF_GENES = sys.argv[4] if len(sys.argv) > 4 else ""
if GTF_GENES:
    gt = pd.read_csv(GTF_GENES, sep="\t")
    gt = gt[gt.gene_biotype == "protein_coding"]
    inter = ok[ok.gene_biotype == "intergenic"]
    fixes = {}
    for chrom, grp in inter.groupby("chrom"):
        cand = gt[gt.chrom.astype(str) == str(chrom)]
        if not len(cand):
            continue
        for r in grp.itertuples():
            mid = (r.start + r.end) // 2
            d = ((cand.start - mid).abs()).combine(
                (cand.end - mid).abs(), min)
            inside = (cand.start <= mid) & (cand.end >= mid)
            d = d.where(~inside, 0)
            j = d.idxmin()
            fixes[r.Index] = (cand.at[j, "gene_id"], cand.at[j, "gene_name"],
                              int(d.loc[j]))
    for idx, (gid, gname, dist) in fixes.items():
        ok.at[idx, "gene_id"] = gid
        ok.at[idx, "gene_name"] = gname
        ok.at[idx, "gene_biotype"] = "protein_coding"
        ok.at[idx, "nearest_dist"] = dist
    print(f"resolved {len(fixes)} intergenic hits to a nearest gene",
          file=sys.stderr)
if "nearest_dist" not in ok:
    ok["nearest_dist"] = pd.NA


def join(series):
    return ";".join(dict.fromkeys(str(x) for x in series if str(x) not in ("", "nan")))


# ---------------------------------------------------------------- per guide --
coding = ok[ok.gene_biotype == "protein_coding"]
rows = []
for name, g in ok.groupby("name"):
    c = coding[coding.name == name]
    pick = c if len(c) else g
    rows.append({
        "guide": name,
        "n_loci": g[["chrom", "start", "strand"]].drop_duplicates().shape[0],
        "gene_id": join(pick.gene_id),
        "gene_name": join(pick.gene_name),
        "gene_biotype": join(pick.gene_biotype),
        "locus": join(f"{r.chrom}:{r.start}-{r.end}({r.strand})"
                      for r in g.drop_duplicates(["chrom", "start", "strand"])
                      .itertuples()),
        "protospacer": join(g.protospacer),
        "PAM": join(g.PAM),
        "g_substituted": int(g.g_substituted.max()),
        "nearest_gene_dist": join(g.nearest_dist.dropna().astype(int)),
    })
per_guide = pd.DataFrame(rows)

out = sp[["name", "spacer"]].rename(columns={"name": "guide"}).merge(
    per_guide, on="guide", how="left")
out["n_loci"] = out.n_loci.fillna(0).astype(int)
out["status"] = [
    "non-targeting control" if n == 0 else
    ("unique" if n == 1 else f"multi-locus ({n})")
    for n in out.n_loci]
out = out.fillna("")
out.to_csv(f"{OUT}_guides_to_genes.tsv", sep="\t", index=False)

# ------------------------------------------------- library block structure ---
# The oligos are ordered in fixed blocks of five guides per gene (verified: 331
# blocks, 324 of them naming one gene unanimously). Where a guide's cut site
# also overlaps a second, usually predicted, gene the per-guide call can differ
# from its neighbours; the block majority resolves that without overriding the
# alignment, and both are reported.
def classify(n):
    m = re.match(r"^(.*?)_(?:(pos|neg)_)?(\d+)$", str(n))
    if not m:
        return "target", None
    kind = {"pos": "positive control", "neg": "non-targeting control",
            None: "target"}[m.group(2)]
    return kind, int(m.group(3))


out["class"], out["index"] = zip(*out.guide.map(classify))
# Only the candidate-gene guides are laid out in blocks of five. The positive
# controls cover several genes with a couple of guides each, so blocking them
# would merge distinct genes; they keep their per-guide call.
targets = out[out["class"] == "target"].copy().sort_values("index")
contiguous = list(targets["index"]) == list(range(1, len(targets) + 1))
targets["block"] = (targets["index"] - 1) // 5
sizes = set(targets.groupby("block").size())
print(f"{len(targets)} candidate-gene guides, contiguous index: {contiguous}; "
      f"{targets.block.nunique()} blocks, block sizes {sizes}", file=sys.stderr)


guide_block = dict(zip(targets.guide, targets.block))
ok_blocks = ok[ok.name.isin(guide_block)].copy()
ok_blocks["block"] = ok_blocks.name.map(guide_block)


def consensus(grp):
    """Majority gene over the guides in a block.

    Counted per guide, not per hit, so a guide whose cut site overlaps two
    annotations cannot outvote its neighbours. Protein-coding genes win when
    any guide in the block hits one: these blocks target patient candidate
    genes, and several sit inside a lncRNA that overlaps the real target - at
    block 73 all five guides fall in a 73 bp window shared by Glra4 and the
    lncRNA BC065397.
    """
    for subset in (grp[grp.gene_biotype == "protein_coding"], grp):
        c = Counter()
        for gene, guides in subset.groupby("gene_name").name:
            c[gene] = len(set(guides))
        if c:
            return c.most_common(1)[0]
    return ("", 0)


block_gene = {b: consensus(grp) for b, grp in ok_blocks.groupby("block")}
out["block"] = out.guide.map(dict(zip(targets.guide, targets.block)))
out["block_consensus_gene"] = out.block.map(
    lambda b: block_gene.get(b, ("", 0))[0] if pd.notna(b) else "")
out["block_votes"] = out.block.map(
    lambda b: block_gene.get(b, ("", 0))[1] if pd.notna(b) else "")
out["agrees_with_block"] = [
    "" if not bc else int(bc in str(gn).split(";"))
    for bc, gn in zip(out.block_consensus_gene, out.gene_name)]

cols = ["guide", "class", "spacer", "protospacer", "g_substituted", "status",
        "n_loci", "locus", "PAM", "gene_id", "gene_name", "gene_biotype",
        "nearest_gene_dist", "block", "block_consensus_gene", "block_votes",
        "agrees_with_block"]
out = out[[c for c in cols if c in out.columns]]
out.to_csv(f"{OUT}_guides_to_genes.tsv", sep="\t", index=False)

# ----------------------------------------------------------------- per gene --
rows = []
for b, (gene, votes) in sorted(block_gene.items()):
    grp = targets[targets.block == b]
    gid = join(ok[ok.gene_name == gene].gene_id) if gene else ""
    rows.append({"gene_name": gene, "gene_id": gid.split(";")[0] if gid else "",
                 "block": b, "n_guides": len(grp),
                 "consensus_votes": f"{votes}/{len(grp)}",
                 "class": "target", "guides": ";".join(grp.guide)})
# positive controls, one row per gene, straight from the alignment
pos = out[out["class"] == "positive control"]
for gene, grp in pos.groupby(pos.gene_name.map(lambda s_: str(s_).split(";")[0])):
    if not gene:
        continue
    rows.append({"gene_name": gene,
                 "gene_id": join(grp.gene_id).split(";")[0],
                 "block": "", "n_guides": len(grp),
                 "consensus_votes": f"{len(grp)}/{len(grp)}",
                 "class": "positive control", "guides": ";".join(grp.guide)})
blocks = pd.DataFrame(rows)
per_gene = (blocks.groupby(["gene_name", "gene_id", "class"], as_index=False)
            .agg(n_blocks=("block", "count"), n_guides=("n_guides", "sum"),
                 guides=("guides", lambda v: ";".join(v))))
per_gene = per_gene.sort_values(["class", "gene_name"])
per_gene.to_csv(f"{OUT}_genes.tsv", sep="\t", index=False)
blocks.to_csv(f"{OUT}_blocks.tsv", sep="\t", index=False)

print(f"\n{len(out)} guides: " + ", ".join(
    f"{k} {v}" for k, v in out["class"].value_counts().items()), file=sys.stderr)
cand = per_gene[per_gene["class"] == "target"]
print(f"candidate genes recovered: {len(cand)}", file=sys.stderr)
print(f"positive-control genes: "
      f"{sorted(per_gene[per_gene['class'] == 'positive control'].gene_name)}",
      file=sys.stderr)
frac = blocks.consensus_votes.str.split("/", expand=True).astype(int)
weak = blocks[(frac[0] / frac[1]) < 0.8]
if len(weak):
    print(f"blocks whose gene call is below 4/5: {len(weak)}", file=sys.stderr)
    for r in weak.itertuples():
        print(f"    block {r.block} -> {r.gene_name} ({r.consensus_votes})",
              file=sys.stderr)
print(f"wrote {OUT}_guides_to_genes.tsv, {OUT}_genes.tsv, {OUT}_blocks.tsv",
      file=sys.stderr)
