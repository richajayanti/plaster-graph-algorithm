# Plaster (graph edition)

A modification of **plaster** (the Treangen Lab (Rice University)'s plasmid pangenome tool) that
builds a **pangenome graph** instead of a linear concatenated pangenome (a list of sequences). Instead
of reading a flat FASTA, you can *see* which stretches of sequence are shared
between genomes and where individual genomes diverge.

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
  - `parse_coords_file()` — parse `show-coords` alignment output into aligned
    query intervals (used by the alignment-based roadmap below).

## Current limitation: exact-match merging

In this build, fragments merge into one node only when their sequences are
**exactly identical**. That works for cleanly pre-fragmented inputs, but real
biological homology is usually *partial* — two fragments sharing a common prefix,
or aligning at 98% identity rather than 100%. Capturing that requires a sequence
aligner.

The original Plaster pipeline handles this with
[MUMmer](https://github.com/mummer4/mummer): align each genome to the growing
reference with `nucmer`/`dnadiff`, use `show-coords` to split each sequence into
aligned and unaligned fragments, and only then build nodes. That machinery is
present here (`run_nucmer`, `parse_coords_file`, and the `--align-only` /
`--realign` modes) but is **not yet wired into graph construction** — MUMmer is
Linux-only and is not invoked during the default graph build. To use it, run
under Linux/WSL with MUMmer installed, or swap in a pure-Python aligner.


## Acknowledgements

Developed by Richa Jayanti under the mentorship of Dr. Todd
Treangen. The first iteration of Plaster was originally created by members of the 
Treangen Lab — Qi Wang, Bryce Kille, Tian Rui Liu, R. A. Leo Elworth, and Todd Treangen.
