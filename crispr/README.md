# Reverse-engineering the patient-candidate knockdown gRNA library

Maps the oligos in `2_patient_cand_051224.xlsx` back to the genes they target.
Run order: `extract_guides.py` → `align_guides.py` → `summarize.py`.

## Result

| | |
|---|---|
| Oligos in the sheet | 1675 (all matched, all spacers unique) |
| Candidate-gene guides | 1645, in 329 blocks of 5 |
| **Candidate genes recovered** | **307** (285 with 5 guides, 22 with 10) |
| Positive controls | 10 guides over Cep290, Crx, Nmnat1, Otx2, Rpgrip1 |
| Non-targeting controls | 20 (`..._neg_1656`–`1675`), correctly unmappable |
| Guides mapped to the genome | 1655 of 1655 targeting + control guides |
| Hits with a canonical NGG PAM | every hit used for a gene call |

Genome: Ensembl GRCm39, release 116. Species comes from the data, not an
assumption: the sheet is `gRNA_mm_ptcand_AAV_GFP_bc_05102` and every oligo is
`mm_ptcand_AAV_*`.

Outputs in `results/`:

| File | Contents |
|---|---|
| `mm_ptcand_guides_to_genes.tsv` | one row per guide: spacer, genomic protospacer, locus, PAM, gene ID/name, block call |
| `mm_ptcand_genes.tsv` | one row per gene: gene ID, guide count, guide names |
| `mm_ptcand_blocks.tsv` | one row per block of five, with the consensus vote |

## Two things the construct does that break a naive extraction

**The stated flanks are not the ones in the file.** The 5' flank given as
`TATATATCTTGTGGAAAGGACGAAACACC` appears in **0 of 1675** oligos; the real
constant is `...CACTTGCTCGTACGACGTATCTTGTGGAAAGGACGAAACACC`. The 3' side is the
modified scaffold `GTTTAAGAGCTATGCTGG`, not the classic
`GTTTTAGAGCTAGAAATAGCAAG`, which is absent. Anchoring on the real pair matches
every oligo.

**The 5' base is substituted, not prepended.** U6 needs a G at the
transcription start. This library obtains it by replacing the spacer's first
base, so an exact 20 nt search finds only guides whose target already begins
with G — 561 of 1675. At the *Rpl10* locus the genome reads
`TCTGGCGCGGTGAGTGAGTT` while the oligo carries `GCTGGCGCGGTGAGTGAGTT`: one
mismatch, at position 1. Searching positions 2–20 and recovering the real 5'
base from the genome lifts this to 1655 of 1655. 1101 guides turn out to carry
a substituted base, and each row records which.

The giveaway was the PAM. The 561 guides that did match were plainly genuine —
563 of 577 hits carried an NGG and sat in protein-coding genes — but they
averaged two guides per gene where the library clearly held five. A real
library with a systematic hole, not a weak one. Recording the PAM instead of
filtering on it is what made that visible.

## Method notes

- **Exact match on positions 2–20, both strands.** No aligner binary is
  available and none is needed: a guide matches its target exactly. A 19 nt
  query is still specific — about 0.02 chance occurrences in a 2.7 Gb genome
  across both strands. Implemented as a single Aho-Corasick pass.
- **PAM is recorded, never used as a filter for reporting.** A 20 nt match with
  no NGG beside it is a coincidental genomic match rather than a target site.
  Gene calls use only NGG hits; the hits file keeps the rest.
- **Every locus is kept.** 1634 guides hit exactly one locus, 21 hit more; those
  are reported with all loci rather than assigned to one arbitrarily.
- **Block consensus.** The oligos are ordered in fixed blocks of five per gene
  (verified: 329 blocks, contiguous indices, 5 guides each). Where a cut site
  overlaps two annotations the block majority resolves it, counted per guide so
  one ambiguous guide cannot outvote its neighbours, preferring the
  protein-coding gene. Only 2 of 329 blocks fall below 4/5, both where a
  protein-coding gene and a lncRNA overlap (`Glra4`/`BC065397` within a 73 bp
  window, `Erich2`/`Gad1os`); both are flagged in `mm_ptcand_blocks.tsv`.
- **Intergenic hits** are resolved to the nearest protein-coding gene with the
  distance recorded, rather than left blank — 43 hits.

## Consistent with a knockdown library

The five guides of a block sit within a median of 210 bp of each other
(325 of 328 within 1 kb), i.e. tightly clustered on one small window per gene
rather than spread across coding exons. That is what TSS-directed knockdown
design looks like, and matches the "KD gRNA" description.

## Caveat

Gene assignment is by genomic position, so it says where each guide cuts, not
how well it silences. Guides are assigned to overlapping genes where the
annotation overlaps (105 guides name more than one gene); for a knockdown
readout the block consensus column is the one to use.
