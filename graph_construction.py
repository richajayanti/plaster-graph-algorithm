# Richa Jayanti, Spring 2025
# Last modified: 4/13/2025

import networkx as nx

def initialize_graph():
    """
    Initializes and returns the data structures used to build a pangenome graph.

    Returns:
    - G: a directed graph (DiGraph) from networkx where:
        • Each node is a DNA fragment
        • Each edge indicates the sequential order of fragments within a genome
    - genome_paths: a dictionary mapping genome IDs to an ordered list of fragment IDs,
      representing each genome’s unique walk through the graph
    """
    G = nx.DiGraph()         # Directed: fragment order matters (e.g., genome1: frag_1 → frag_2)
    genome_paths = {}        # Maps each genome ID to the list of fragments it passes through
    return G, genome_paths


def add_fragment_node(G, frag_id, sequence, genome_id):
    """
    Adds a fragment node to the graph G. If the sequence already exists,
    the genome ID is appended to that node's genome tracking field.

    Parameters:
    - G: the graph object
    - frag_id: a unique ID string for this fragment (e.g., "frag_4")
    - sequence: the DNA string of this fragment
    - genome_id: which genome this fragment belongs to

    Note:
    - This function assumes that fragments are uniquely identified by frag_id.
    - In your deduplication logic, identical sequences across genomes should share frag_id.
    """
    if frag_id not in G:
        # Create a new node with the sequence and origin genome.
        # 'genomes' is stored as a comma-joined string so it survives GML/
        # GraphML export (those formats can't store Python lists) and matches
        # what visualize_graph.py expects.
        G.add_node(frag_id, sequence=sequence, genomes=genome_id)
    else:
        # Fragment already exists in the graph.
        # Add genome_id to its 'genomes' field if not already recorded.
        genomes = G.nodes[frag_id]["genomes"].split(",")
        if genome_id not in genomes:
            genomes.append(genome_id)
            G.nodes[frag_id]["genomes"] = ",".join(genomes)


def update_genome_path(G, genome_paths, frag_id, genome_id):
    """
    Updates the graph with genome-specific path information:
    - Appends this fragment to the genome's current path
    - Adds a directed edge from the previous fragment to the current one

    Parameters:
    - G: the graph
    - genome_paths: maps each genome to its ordered list of fragments
    - frag_id: the current fragment being added
    - genome_id: the genome being processed

    Example:
    If genome1 has frag_1 followed by frag_2, this creates an edge:
    frag_1 → frag_2 (genome=genome1)
    """
    if genome_id not in genome_paths:
        genome_paths[genome_id] = []

    path = genome_paths[genome_id]

    # If not the first fragment, link previous to current.
    # A single edge may be traversed by several genomes, so accumulate the
    # supporting genome ids in a comma-joined string (same convention as nodes).
    if path:
        prev_id = path[-1]
        if G.has_edge(prev_id, frag_id):
            genomes = G.edges[prev_id, frag_id]["genome"].split(",")
            if genome_id not in genomes:
                genomes.append(genome_id)
                G.edges[prev_id, frag_id]["genome"] = ",".join(genomes)
        else:
            G.add_edge(prev_id, frag_id, genome=genome_id)

    # Add current fragment to genome's walk
    path.append(frag_id)


def export_graph(G, output_file_stem, genome_paths=None):
    """
    Exports the graph to disk in multiple formats for downstream analysis.

    Parameters:
    - G: the constructed graph object
    - output_file_stem: file path prefix (e.g., 'real_output') for output files
    - genome_paths: optional {genome_id: [frag_id, ...]} mapping. When given,
      a GFA file (with per-genome P-lines) is also written for visual tools
      such as Bandage.

    Output files:
    - .gml: readable text format, good for inspection/debugging
    - .graphml: XML-based format, compatible with Cytoscape and Gephi
    - .gfa: pangenome graph format (only when genome_paths is provided)
    """
    nx.write_gml(G, f"{output_file_stem}.gml")
    nx.write_graphml(G, f"{output_file_stem}.graphml")
    if genome_paths is not None:
        export_gfa(G, output_file_stem, genome_paths)


def export_gfa(G, output_file_stem, genome_paths):
    """
    Writes the pangenome graph as GFA 1.0 with rGFA-style segment tags.

    Each fragment becomes an S(egment) line, each graph edge an L(ink) line,
    and each genome an ordered P(ath) line. Segments carry rGFA tags:
    - SN:Z  stable name  = the genome that first introduced the fragment
    - SO:i  stable offset = base offset within that genome's walk
    - SR:i  rank          = 0 for the first-processed genome, else 1

    The result opens directly in Bandage and other GFA viewers so you can
    see shared fragments as branch/merge points instead of a linear string.
    """
    # Assign rGFA tags by walking genomes in processing order. The first
    # genome to introduce a fragment "owns" its stable name/offset.
    reference = next(iter(genome_paths), None)
    seg_tags = {}
    for genome_id, path in genome_paths.items():
        offset = 0
        for frag_id in path:
            seq_len = len(G.nodes[frag_id]["sequence"])
            if frag_id not in seg_tags:
                seg_tags[frag_id] = {
                    "SN": genome_id,
                    "SO": offset,
                    "SR": 0 if genome_id == reference else 1,
                }
            offset += seq_len

    with open(f"{output_file_stem}.gfa", "w") as gfa:
        gfa.write("H\tVN:Z:1.0\n")
        for frag_id, data in G.nodes(data=True):
            seq = data["sequence"] or "*"
            tags = seg_tags.get(frag_id, {"SN": data["genomes"].split(",")[0], "SO": 0, "SR": 0})
            gfa.write(
                "S\t{}\t{}\tLN:i:{}\tSN:Z:{}\tSO:i:{}\tSR:i:{}\n".format(
                    frag_id, seq, len(data["sequence"]),
                    tags["SN"], tags["SO"], tags["SR"]))
        for src, dst in G.edges():
            # All fragments are stored on the forward strand, so links are +/+
            # with no overlap (0M).
            gfa.write("L\t{}\t+\t{}\t+\t0M\n".format(src, dst))
        for genome_id, path in genome_paths.items():
            if path:
                walk = ",".join(f"{frag_id}+" for frag_id in path)
                gfa.write("P\t{}\t{}\t*\n".format(genome_id, walk))


def parse_coords_file(coords_file):
    """
    Parse a `show-coords -qT` output file into aligned intervals, grouped by the
    query sequence (FASTA record) they belong to.

    In tab-delimited (-T) show-coords output each alignment row is:
        [S1] [E1] [S2] [E2] [LEN1] [LEN2] [%IDY] ... [REF_TAG] [QUERY_TAG]
    where S1/E1 are 1-based reference coordinates, S2/E2 are 1-based query
    coordinates, and the last two columns are the reference/query sequence names.
    Header lines (the input paths, "NUCMER", the column titles) do not start with
    an integer, so they are skipped automatically.

    Parameters:
    - coords_file: path to a `.coords_out` file produced by show-coords -qT

    Returns:
    - by_query: dict mapping each query record id to a list of
      (q_start, q_end, ref_tag) tuples, 1-based inclusive, sorted by q_start.
      q_start <= q_end always (reverse-strand hits are normalized). ref_tag is
      the reference fragment the query region aligned to, which is what lets the
      caller reuse (merge into) an existing reference node.
    """
    by_query = {}
    with open(coords_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            # Data rows begin with four integer coordinates; anything else
            # (headers, blank lines) fails this check and is skipped.
            try:
                q_start, q_end = int(parts[2]), int(parts[3])
            except ValueError:
                continue
            ref_tag = parts[-2]
            query_tag = parts[-1]
            lo, hi = sorted((q_start, q_end))
            by_query.setdefault(query_tag, []).append((lo, hi, ref_tag))
    for tag in by_query:
        by_query[tag].sort()
    return by_query


def split_query_record(sequence, intervals):
    """
    Cut one query sequence into aligned and unaligned fragments.

    This is the step that turns an alignment into graph fragments: the regions
    that aligned to the reference become fragments that will be *merged* into the
    matching reference node, and the gaps between them become new,
    genome-specific fragments.

    Parameters:
    - sequence: the query record's DNA string
    - intervals: list of (q_start, q_end, ref_tag) as returned by
      parse_coords_file for this record (1-based inclusive query coordinates)

    Returns:
    - fragments: ordered list of (subsequence, is_aligned, ref_tag) tuples whose
      subsequences concatenate back to `sequence`. For unaligned pieces
      is_aligned is False and ref_tag is None. Overlapping alignments are
      resolved greedily left-to-right so every base is emitted exactly once.
    """
    n = len(sequence)
    ivs = sorted(
        ((max(1, s), min(n, e), tag) for (s, e, tag) in intervals if s <= e),
        key=lambda x: (x[0], x[1]),
    )

    fragments = []
    cursor = 1  # 1-based index of the next base not yet emitted
    for start, end, ref_tag in ivs:
        if end < cursor:
            continue           # already covered by an earlier alignment
        if start < cursor:
            start = cursor     # trim the overlap with the previous alignment
        if start > cursor:
            fragments.append((sequence[cursor - 1:start - 1], False, None))
        fragments.append((sequence[start - 1:end], True, ref_tag))
        cursor = end + 1

    if cursor <= n:
        fragments.append((sequence[cursor - 1:], False, None))
    return fragments
