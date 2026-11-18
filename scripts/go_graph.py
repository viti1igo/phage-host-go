"""Merge labels into unique genes, build each gene's GO subtree, and join them into one graph per group.

A label such as gene=dnaA_2 is one copy of gene=dnaA; Prokka adds the _N suffix when a
genome has several copies. Labels that differ only in that suffix become one gene, whose GO
terms are the union of theirs. Labels written as protein=... are cut-off product names
(protein=putative, protein=26), so they are not merged and stay as they are.

A gene's subtree is its listed terms plus every ancestor under --relations. The graph of a
group is the union of its genes' subtrees, with the GO edges between those terms. Run from
the repo root, one graph per group:

  python scripts/go_graph.py data/go_data/go-basic.obo code/富集分析data_process.ipynb data/go_data/graph \\
      --group host=model/host_key_gene_0.00004.csv --group phage=model/phage_key_gene_0.0001.csv

For Bakta runs (one folder per genome), --bakta reads the GO terms and gene names straight
from the Bakta tables. A gene is the Bakta gene name (column Gene); a feature without one is
keyed by its UniRef90 cluster, and failing that by its own ID. Copies of a gene, in one genome
or across genomes, merge into one gene. Writes the group "host" (or --name) for all genomes together, plus
one graph per genome:

  python scripts/go_graph.py data/go_data/go-basic.obo data/go_data/gaborieau/host/bakta_filled data/go_data/gaborieau/host/graph \\
      --bakta --genomes ../coli_phage_interactions_2023/data/interactions/interaction_matrix.csv

  hosts/GENOME.graphml  the graph of one genome's genes (node counts are genes of that genome)
  hosts.tsv             per genome: features with GO, genes, graph terms and edges

To draw a chosen set of proteins from many genomes instead, --select FILE names them, one per
line: genome|locus tag (one protein), a gene key such as ompC or UniRef90_... (that gene in every
genome), or genome|gene key. It writes one graph, NAME, and no per-genome graphs:

  python scripts/go_graph.py data/go_data/go-basic.obo data/go_data/gaborieau/host/bakta_filled out \\
      --bakta --select receptors.txt --name receptors

Writes to out_folder, for every group NAME:
  genes_NAME.tsv     one row per unique gene: labels merged, listed, deepest and subtree terms, and a
                     QuickGO link that draws the gene's deepest terms (QuickGO's GO release and edges,
                     including regulates, so it can show more than the subtree here)
  subtrees_NAME.tsv  one row per gene and subtree term, marking the gene's deepest terms
  nodes_NAME.tsv     one row per graph term: domain, depths, genes through it, genes ending at it,
                     proteins (labels) through it, and with --bakta the genomes those come from
  edges_NAME.tsv     one row per graph edge, child to parent, with its relation
  graph_NAME.graphml the same nodes and edges, for Cytoscape or Gephi
Standard library only.
"""
import argparse
import collections
import os
import re
import xml.etree.ElementTree as ET

from go_dag import Dag, RELATION_SETS, obo_version, parse_obo
from go_extract import GO_ID, bakta_rows, read_genomes
from go_rollup import deepest, read_gene_list, read_pairs, write_table


def unique_gene(label):
    """gene=dnaA_2 -> gene=dnaA; protein=... labels are left as they are."""
    return re.sub(r"_\d+$", "", label) if label.startswith("gene=") else label


UNIREF90 = re.compile(r"UniRef90_\w+")
QUICKGO_CHART = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}/chart"


def bakta_gene(gene_id, row):
    """Bakta gene name; else the UniRef90 cluster; else the feature's own ID."""
    if row.get("Gene"):
        return row["Gene"]
    cluster = UNIREF90.search(row.get("DbXrefs", ""))
    return cluster.group() if cluster else gene_id


def read_bakta_genes(path, genomes):
    """({feature ID: GO IDs}, {feature ID: gene}) for the Bakta features that list GO terms."""
    pairs, names = {}, {}
    for gene_id, row in bakta_rows(path, genomes):
        found = set(GO_ID.findall(row.get("DbXrefs", "")))
        if found:
            pairs[gene_id], names[gene_id] = found, bakta_gene(gene_id, row)
    return pairs, names


def merge(pairs, labels, dag, key=unique_gene):
    """Return {gene: (merged labels, live terms)} for the labels that have GO terms; key(label) names the gene."""
    genes = collections.defaultdict(lambda: (set(), set()))
    for label in sorted(labels & set(pairs)):
        live = {t for t in map(dag.resolve, pairs[label]) if t}
        if live:
            merged, found = genes[key(label)]
            merged.add(label)
            found |= live
    return dict(genes)


def subtree(listed, dag):
    return set(listed).union(*(dag.ancestors(t) for t in listed))


def graph(genes, dag, terms, relations, pairs=None):
    """Return (node rows, edge rows) for the union of the genes' subtrees.

    Given the {label: GO IDs} pairs behind the genes, nodes also count the labels (proteins) whose
    own subtree passes through the term, and, for genome|locus labels, the genomes they come from.
    """
    through, ending = collections.Counter(), collections.Counter()
    proteins, genomes, cache = collections.Counter(), collections.defaultdict(set), {}
    for merged, listed in genes.values():
        through.update(subtree(listed, dag))
        ending.update(deepest(listed, dag))
        for label in merged if pairs else ():
            key = frozenset(pairs[label])
            if key not in cache:  # copies of a gene usually list the same IDs
                cache[key] = subtree({t for t in map(dag.resolve, key) if t}, dag)
            proteins.update(cache[key])
            if "|" in label:
                for t in cache[key]:
                    genomes[t].add(label.split("|")[0])
    nodes = [{"go_id": t, "name": terms[t]["name"], "domain": terms[t]["domain"],
              "longest_depth": dag.longest_depth(t), "shortest_depth": dag.shortest_depth(t),
              "genes_through": through[t], "genes_deepest": ending[t],
              **({"proteins_through": proteins[t]} if pairs else {}),
              **({"genomes_through": len(genomes[t])} if genomes else {})} for t in sorted(through)]
    edges = sorted({(child, parent, relation) for child in through for relation, parent in terms[child]["edges"]
                    if relation in relations and parent in dag.parents[child]})
    return nodes, [{"child": c, "parent": p, "relation": r} for c, p, r in edges]


def write_graphml(path, nodes, edges):
    root = ET.Element("graphml", xmlns="http://graphml.graphdrawing.org/xmlns")
    keys = [(k, "node", "int" if isinstance(v, int) else "string") for k, v in nodes[0].items() if k != "go_id"]
    keys.append(("relation", "edge", "string"))
    for name, kind, kind_type in keys:
        ET.SubElement(root, "key", {"id": name, "for": kind, "attr.name": name, "attr.type": kind_type})
    body = ET.SubElement(root, "graph", id="G", edgedefault="directed")
    for node in nodes:
        element = ET.SubElement(body, "node", id=node["go_id"])
        for k, v in node.items():
            if k != "go_id":
                ET.SubElement(element, "data", key=k).text = str(v)
    for edge in edges:
        element = ET.SubElement(body, "edge", source=edge["child"], target=edge["parent"])
        ET.SubElement(element, "data", key="relation").text = edge["relation"]
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def write_group(out, name, genes, dag, terms, relations, per_genome=False, pairs=None):
    """Write the five tables and GraphML of one group; return (nodes, edges)."""
    gene_rows, subtree_rows = [], []
    for gene, (merged, listed) in sorted(genes.items()):
        tips, tree = deepest(listed, dag), subtree(listed, dag)
        row = {"gene": gene, "labels_merged": len(merged)}
        if per_genome:
            row["genomes"] = len({label.split("|")[0] for label in merged})
        row.update({"labels": ";".join(sorted(merged)), "listed_terms": len(listed),
                    "deepest_terms": len(tips), "subtree_terms": len(tree),
                    "quickgo_chart": QUICKGO_CHART.format(",".join(sorted(tips)))})
        gene_rows.append(row)
        subtree_rows += [{"gene": gene, "go_id": t, "deepest": "yes" if t in tips else "no"} for t in sorted(tree)]
    nodes, edges = graph(genes, dag, terms, relations, pairs)
    write_table(os.path.join(out, f"genes_{name}.tsv"), gene_rows, "\t")
    write_table(os.path.join(out, f"subtrees_{name}.tsv"), subtree_rows, "\t")
    write_table(os.path.join(out, f"nodes_{name}.tsv"), nodes, "\t")
    write_table(os.path.join(out, f"edges_{name}.tsv"), edges, "\t")
    write_graphml(os.path.join(out, f"graph_{name}.graphml"), nodes, edges)
    return nodes, edges


def read_selection(path):
    """The first word of every line; # starts a comment."""
    with open(path, encoding="utf-8") as handle:
        return {words[0] for words in (line.split("#")[0].split() for line in handle) if words}


def selected(pairs, names, wanted):
    """(labels picked by wanted, entries of wanted that picked nothing).

    An entry picks a protein by genome|locus tag, a gene in every genome by its gene key (Bakta
    gene name such as ompC, or UniRef90_...), or a gene in one genome by genome|gene key.
    """
    picked, used = set(), set()
    for label in pairs:
        for key in (label, names[label], f"{label.split('|')[0]}|{names[label]}"):
            if key in wanted:
                picked.add(label)
                used.add(key)
    return picked, sorted(wanted - used)


def run_bakta(args, dag, terms, relations):
    genomes = read_genomes(args.genomes) if args.genomes else None
    pairs, names = read_bakta_genes(args.pairs, genomes)
    print(f"{len(pairs)} Bakta features with GO in {args.pairs}")
    if args.select:
        labels, unmatched = selected(pairs, names, read_selection(args.select))
        genes = merge(pairs, labels, dag, key=names.get)
        nodes, edges = write_group(args.out, args.name, genes, dag, terms, relations, per_genome=True, pairs=pairs)
        found = len({label.split("|")[0] for label in labels})
        print(f"{args.name}: {len(labels)} proteins from {found} genomes, merged into {len(genes)} genes; "
              f"graph has {len(nodes)} terms and {len(edges)} edges")
        if unmatched:
            print(f"{len(unmatched)} entries of {args.select} match no protein with GO: {', '.join(unmatched[:20])}"
                  + (" ..." if len(unmatched) > 20 else ""))
        return
    genes = merge(pairs, set(pairs), dag, key=names.get)
    nodes, edges = write_group(args.out, args.name, genes, dag, terms, relations, per_genome=True, pairs=pairs)
    print(f"{args.name}: merged into {len(genes)} genes; graph has {len(nodes)} terms and {len(edges)} edges")

    by_genome = collections.defaultdict(set)
    for label in pairs:
        by_genome[label.split("|")[0]].add(label)
    os.makedirs(os.path.join(args.out, args.name + "s"), exist_ok=True)
    rows = []
    for genome, labels in sorted(by_genome.items()):
        own = merge(pairs, labels, dag, key=names.get)
        nodes, edges = graph(own, dag, terms, relations, pairs)
        write_graphml(os.path.join(args.out, args.name + "s", f"{genome}.graphml"), nodes, edges)
        rows.append({"genome": genome, "features_with_go": len(labels), "genes": len(own),
                     "terms": len(nodes), "edges": len(edges)})
    write_table(os.path.join(args.out, args.name + "s.tsv"), rows, "\t")
    missing = sorted(genomes - set(by_genome)) if genomes else []
    print(f"{len(rows)} genome graphs in {os.path.join(args.out, args.name + 's')}"
          + (f"; no Bakta GO for {', '.join(missing)}" if missing else ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("obo")
    parser.add_argument("pairs", help="notebook or gene<TAB>GO table, as for go_rollup.py; with --bakta, "
                                      "a folder of Bakta runs")
    parser.add_argument("out")
    parser.add_argument("--group", action="append", metavar="NAME=CSV",
                        help="labels of one group, from a CSV with a 'gene' column")
    parser.add_argument("--bakta", action="store_true", help="pairs is a folder of Bakta runs")
    parser.add_argument("--name", default="host", help="with --bakta: group name (writes NAME, NAMEs/, NAMEs.tsv)")
    parser.add_argument("--genomes", help="with --bakta: read only the genomes in the first column of this "
                                          "table (header row skipped)")
    parser.add_argument("--select", metavar="FILE",
                        help="with --bakta: draw only the proteins this file names, one per line: genome|locus tag, "
                             "a gene key (ompC, UniRef90_...) for every genome, or genome|gene key. Writes one graph "
                             "NAME and no per-genome graphs")
    parser.add_argument("--relations", choices=list(RELATION_SETS), default="is_a+part_of")
    args = parser.parse_args()
    if args.select and not args.bakta:
        parser.error("--select needs --bakta")
    os.makedirs(args.out, exist_ok=True)

    terms, alt = parse_obo(args.obo)
    relations = RELATION_SETS[args.relations]
    dag = Dag(terms, alt, relations)
    print(f"GO {obo_version(args.obo)}, edges {args.relations}")
    if args.bakta:
        return run_bakta(args, dag, terms, relations)
    if not args.group:
        parser.error("--group is required unless --bakta is given")
    pairs = read_pairs(args.pairs)
    print(f"{len(pairs)} labels with GO in {args.pairs}")
    for spec in args.group:
        name, path = spec.split("=", 1)
        labels = read_gene_list(path)
        genes = merge(pairs, labels, dag)
        nodes, edges = write_group(args.out, name, genes, dag, terms, relations, pairs=pairs)
        used = sum(len(m) for m, _ in genes.values())
        print(f"{name}: {len(labels)} labels, {used} with live GO, merged into {len(genes)} genes; "
              f"graph has {len(nodes)} terms and {len(edges)} edges")


if __name__ == "__main__":
    main()
