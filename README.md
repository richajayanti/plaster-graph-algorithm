# plaster (graph edition)

A modification of **plaster** — the Treangen Lab (Rice University)'s tool for
comparing the DNA of related organisms — that summarizes a collection of genomes
as a **graph** (a network of connected nodes) rather than as a long list of
sequences.

In plain terms: a genome is just a long string of DNA letters, and a group of
related organisms will share many stretches of that string while differing in
others. The older approach lines the genomes up against one reference and reports
which pieces are present or missing — but that flattens everything into a single
straight line and hides *how* the genomes are rearranged relative to each other.
This version instead builds a network where each distinct stretch of DNA is a
node and each genome is a route through those nodes. Shared stretches become
nodes that many routes pass through; differences show up as the routes splitting
apart and rejoining. The result lets you *see*, at a glance, what the genomes
have in common and where they diverge.

## Biological background and motivation

A **pangenome** is the full set of sequences found across a group of related
genomes (typically isolates of the same species). It is conventionally split
into two compartments:

- **Core genome** — sequence present in *every* genome in the set. These regions
  tend to encode essential, conserved functions.
- **Accessory genome** — sequence present in *some but not all* genomes. This is
  where strain-to-strain diversity lives: acquired genes, mobile elements,
  insertions, and deletions.

Traditional pangenome analysis compares genomes against a single linear
reference and reports gene *presence/absence*. That framing loses **structural
variation** — rearrangements, duplications, and alternative paths — because a
linear coordinate system can only describe one layout at a time. A **graph**
representation keeps every layout at once: shared sequence becomes a single node
that multiple genomes pass through, and structural differences appear as
branches and merges rather than being flattened away.

Plaster was originally built to construct pangenomes from **plasmids** — mobile,
often circular DNA elements that are a major vehicle for horizontal gene
transfer (e.g. antibiotic-resistance cassettes) and are consequently rich in the
kind of structural variation a graph captures well.

## What the graph means

Each FASTA record is treated as a sequence **fragment**:

- **Nodes** are unique fragments. A node stores its DNA sequence and the set of
  genomes it appears in. A fragment shared by several genomes collapses into a
  single node — a **core-like** fragment. A fragment seen in only one genome is
  an **accessory-like** fragment.
- **Edges** are directed and encode *fragment order within a genome*: an edge
  `frag_A → frag_B` means some genome traverses A immediately before B. Each edge
  records which genomes support it, so an edge shared by many genomes is a
  conserved junction, while a lone edge marks a genome-specific rearrangement.
- **Paths** — each genome's ordered walk through the graph — are preserved and
  exported so tools like Bandage can highlight one isolate's route through the
  shared structure.

Because identical fragments merge, the output stops being a set of parallel
linear paths and becomes a real graph with **branch points** (one fragment
leading to different successors in different genomes) and **merge points** (one
fragment reached from different predecessors).

## Example

Running on three toy genomes (`genome1-3.fasta`) produces this graph — red nodes
are fragments shared across genomes (core-like), gray nodes are unique to one
genome (accessory-like):

![Pangenome graph for genome1-3](real_output_graph.png)

`frag_2` (`ATGCGTACGTTAG`) is shared by genome1 and genome3 and branches to two
different next fragments — a point of structural divergence. `frag_1`
(`GGGCTTAGACCAA`) is a merge point reached by both genome1 and genome2.

A larger, more realistic run on `exp_train1-6.fasta`:

![Pangenome graph for exp_train1-6](exp_output_graph.png)

## Requirements

- Python 3.8+
- [Biopython](https://biopython.org/), [networkx](https://networkx.org/), [tqdm](https://github.com/tqdm/tqdm)
- [matplotlib](https://matplotlib.org/) (only for the visualization scripts)

```bash
pip install -r requirements.txt
```

## Usage

Build a graph from one or more FASTA files (each file is one genome; each record
is one fragment):

```bash
python plaster genome1.fasta genome2.fasta genome3.fasta -o real_output
```

This writes three representations of the same graph:

| File | Format | Use |
|------|--------|-----|
| `real_output.gml`      | GML          | Human-readable, easy to inspect |
| `real_output.graphml`  | GraphML      | Opens in Cytoscape / Gephi |
| `real_output.gfa`      | GFA 1.0 (rGFA-style tags) | Opens in [Bandage](https://rrwick.github.io/Bandage/) and other pangenome viewers |

The GFA export carries rGFA-style segment tags (`SN` stable name, `SO` stable
offset, `SR` rank) so the first-processed genome acts as the reference backbone
and later genomes are placed relative to it — the same convention used by
reference-graph pangenome tools.

Common options (`python plaster --help` for the full list):

- `-o, --output` &nbsp; output file stem (no extension)
- `-t, --template` &nbsp; seed genome to start the pangenome from (chosen as the
  reference backbone; by default the largest input is used, a proxy for the
  longest genome)
- `--align` &nbsp; merge fragments by **homology** instead of exact match
  (requires MUMmer; see [Two ways to merge fragments](#two-ways-to-merge-fragments))
- `-l, --length` &nbsp; minimum alignment length passed to `nucmer` (only used
  with `--align`; must be **shorter than your fragments** or nothing aligns)
- `-p, --threads` &nbsp; number of threads
- `-v, --verbose` &nbsp; print each `record -> fragment` assignment

## Visualizing

`visualize_graph.py` renders `real_output.graphml` with matplotlib, coloring
shared (core-like) fragments red and unique (accessory-like) fragments gray:

```bash
python visualize_graph.py
```

## How it works

- `plaster` — command-line entry point. Orders inputs by size (proxy for genome
  length), seeds the pangenome from a reference genome, reads each FASTA,
  deduplicates fragments by sequence, and grows the graph one genome at a time.
- `graph_construction.py` — the graph itself:
  - `initialize_graph()` — create the empty directed graph and the per-genome
    path tracker.
  - `add_fragment_node()` — add a fragment node and record which genomes carry
    it (accumulated as a comma-joined list so it survives GML/GraphML export).
  - `update_genome_path()` — append the fragment to the current genome's walk
    and add the directed edge from the previous fragment, tracking edge support.
  - `export_graph()` / `export_gfa()` — write the GML, GraphML, and GFA outputs.
  - `parse_coords_file()` — parse `show-coords -qT` output into aligned query
    intervals, grouped by record (used by `--align`).
  - `split_query_record()` — cut one query sequence into aligned/unaligned
    fragments given those intervals (used by `--align`).

## Two ways to merge fragments

The whole point of a pangenome graph is that shared sequence becomes a *single
shared node*. plaster can decide "shared" in two ways.

### 1. Exact match (default — no external tools)

By default, two fragments merge into one node only when their sequences are
**byte-for-byte identical**. This needs nothing but Python and works on any
platform (including Windows).

- ✅ Correct and fast for **pre-fragmented input**, where shared regions are
  already emitted as identical records (this is the case for the bundled
  `genome*` and `exp_train*` examples).
- ⚠️ On **whole-genome input** it under-connects: two related genomes that differ
  by even one base won't merge, so the graph degenerates toward parallel linear
  paths. A merge always means true identity, so any shared node you see is real —
  but *absence* of sharing can't be trusted.

### 2. Homology (`--align` — requires MUMmer)

With `--align`, fragments merge by **alignment** instead of exact identity, so
near-identical regions (SNPs, small indels, shifted boundaries) still collapse
into one node. This is closer to how real pangenome graphs are built.

```bash
python plaster genomeA.fasta genomeB.fasta genomeC.fasta -o out --align -l 20
```

How it works (`build_graph_with_alignment` in `plaster`): the largest genome
seeds a reference; each later genome is aligned to the growing reference with
`nucmer`, and every record is cut by `split_query_record()` into aligned pieces
(merged into the reference fragment they hit) and unaligned pieces (added as new,
genome-specific fragments that later genomes can in turn align to).

**Requirements & caveats:**

- Needs [MUMmer](https://github.com/mummer4/mummer) (`nucmer`, `dnadiff`,
  `show-coords`) on your `PATH`. MUMmer is **Linux/macOS only** — on Windows use
  WSL. Without it, `--align` exits with a clear error and changes nothing.
- Pass `-l` **smaller than your fragment lengths** (e.g. `-l 20` for ~45 bp
  fragments); `nucmer`'s default cluster size will otherwise discard short hits.
- v1 reuses a reference fragment when a query aligns to it, but does **not** yet
  sub-split a reference fragment that's only partially covered, so fragments with
  very different boundaries may not merge perfectly. Edge weights (alignment
  identity) are also not recorded yet.

> **Testing note:** the fragment-cutting and coords-parsing logic
> (`split_query_record`, `parse_coords_file`) and the end-to-end merge logic in
> `build_graph_with_alignment` are covered by unit/simulation tests that run
> without MUMmer. The live `nucmer` invocation itself has **not** been run on
> this (Windows) machine — validate the first real `--align` run on a Linux/WSL
> box with MUMmer installed.


## Acknowledgements

Developed by Richa Jayanti under the mentorship of Dr. Todd
Treangen. The first iteration of Plaster was originally created by members of the 
Treangen Lab — Qi Wang, Bryce Kille, Tian Rui Liu, R. A. Leo Elworth, and Todd Treangen.
