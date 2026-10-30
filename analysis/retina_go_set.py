#!/usr/bin/env python3
"""
Write the set of retina / vision GO terms: the curated roots below plus every
descendant of them in the GO DAG (is_a + part_of).

Defining the set from the ontology rather than by matching words in term names
means a term like "9+0 non-motile cilium" is excluded (it is not under an eye
root) while "photoreceptor connecting cilium" is included, and nothing depends
on the wording a term happens to use.

Usage: retina_go_set.py go-basic.obo out.tsv
"""
import sys
from collections import defaultdict

OBO, OUT = sys.argv[1], sys.argv[2]

# Curated entry points into the eye/vision part of the ontology. Deliberately
# excludes "sensory organ development" (GO:0007423), which also covers ear,
# nose and lateral line, and "sensory perception" generally.
ROOTS = [
    # development / maintenance
    "GO:0001654",  # eye development
    "GO:0043010",  # camera-type eye development
    "GO:0150063",  # visual system development
    "GO:0060041",  # retina development in camera-type eye
    "GO:0048592",  # eye morphogenesis
    "GO:0048593",  # camera-type eye morphogenesis
    "GO:0001895",  # retina homeostasis
    "GO:0045494",  # photoreceptor cell maintenance
    "GO:0046530",  # photoreceptor cell differentiation
    "GO:0042462",  # eye photoreceptor cell development
    # function / signalling
    "GO:0007601",  # visual perception
    "GO:0050953",  # sensory perception of light stimulus
    "GO:0007602",  # phototransduction
    "GO:0016056",  # G protein-coupled opsin signaling pathway
    "GO:0016038",  # absorption of visible light
    # retinoid / visual cycle
    "GO:0001523",  # retinoid metabolic process
    "GO:0042572",  # retinol metabolic process
    "GO:0042574",  # retinal metabolic process
    "GO:0016918",  # retinal binding
    # structures
    "GO:0001750",  # photoreceptor outer segment
    "GO:0001917",  # photoreceptor inner segment
    "GO:0097733",  # photoreceptor cell cilium
]


def parse(path):
    parents, names, nsp, alt = defaultdict(set), {}, {}, {}
    cur = None
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {"id": None}
            continue
        if line.startswith("["):
            cur = None
            continue
        if cur is None or not line:
            continue
        k, _, v = line.partition(": ")
        if k == "id":
            cur["id"] = v
        elif k == "name":
            names[cur["id"]] = v
        elif k == "namespace":
            nsp[cur["id"]] = v
        elif k == "alt_id":
            alt[v] = cur["id"]
        elif k == "is_a":
            parents[cur["id"]].add(v.split(" ! ")[0].strip())
        elif k == "relationship":
            rel, _, tgt = v.partition(" ")
            if rel == "part_of":
                parents[cur["id"]].add(tgt.split(" ! ")[0].strip())
    return parents, names, nsp, alt


parents, names, nsp, alt = parse(OBO)
children = defaultdict(set)
for c, ps in parents.items():
    for p in ps:
        children[p].add(c)

roots = [alt.get(r, r) for r in ROOTS]
missing = [r for r in roots if r not in names]
if missing:
    sys.exit(f"root terms absent from this ontology release: {missing}")

keep, stack = set(roots), list(roots)
while stack:
    t = stack.pop()
    for c in children[t]:
        if c not in keep:
            keep.add(c)
            stack.append(c)

NS_SHORT = {"biological_process": "BP", "molecular_function": "MF",
            "cellular_component": "CC"}
with open(OUT, "w") as fh:
    fh.write("GO_ID\tontology\tterm\tis_root\n")
    for t in sorted(keep):
        fh.write(f"{t}\t{NS_SHORT.get(nsp.get(t), '?')}\t{names.get(t, '')}"
                 f"\t{int(t in set(roots))}\n")
print(f"{len(roots)} roots -> {len(keep)} retina/vision GO terms", file=sys.stderr)
