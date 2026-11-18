"""List every distinct protein in a folder of Bakta runs, with how many genomes have it and its GO.

Proteins merge across genomes the way go_graph.py --bakta merges them: by Bakta gene name, else by
UniRef90 cluster, else the CDS stands alone (genome|locus tag). One row per protein, most widespread
first. With the go_source.tsv that uniref90_go_fill.py writes next to its tables, the row also
counts where each copy's GO came from.

  python scripts/go_protein_list.py data/go_data/gaborieau/host/bakta_filled proteins.tsv

The protein column is what go_graph.py --select takes.
"""
import argparse
import collections
import csv
import os

from go_extract import GO_ID, bakta_rows, read_genomes
from go_graph import UNIREF90, bakta_gene

SOURCES = ("bakta", "uniref90", "uniref90_unisave", "none")


def proteins(path, genomes=None, source=None):
    """{protein: row dict}, built from every CDS in the Bakta tables."""
    found = {}
    for cds, row in bakta_rows(path, genomes):
        if row.get("Type") != "cds":
            continue
        key = bakta_gene(cds, row)
        p = found.setdefault(key, {"cds": 0, "genomes": set(), "with_go": 0, "go": set(),
                                   "product": collections.Counter(), "uniref90": collections.Counter(),
                                   "source": collections.Counter()})
        go = GO_ID.findall(row.get("DbXrefs", ""))
        p["cds"] += 1
        p["genomes"].add(cds.split("|")[0])
        p["with_go"] += bool(go)
        p["go"].update(go)
        p["product"][row.get("Product", "")] += 1
        cluster = UNIREF90.search(row.get("DbXrefs", ""))
        if cluster:
            p["uniref90"][cluster.group()] += 1
        if source is not None:
            p["source"][source.get(cds, "")] += 1
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bakta", help="folder of <genome>/<genome>.tsv Bakta tables")
    parser.add_argument("out")
    parser.add_argument("--genomes", help="only the genomes in the first column of this table (header row skipped)")
    args = parser.parse_args()
    source = None
    source_file = os.path.join(args.bakta, "go_source.tsv")
    if os.path.exists(source_file):
        with open(source_file) as handle:
            source = {f"{r['genome']}|{r['locus_tag']}": r["go_source"] for r in csv.DictReader(handle, delimiter="\t")}
    found = proteins(args.bakta, read_genomes(args.genomes) if args.genomes else None, source)

    with open(args.out, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t", lineterminator="\n")
        out.writerow(["protein", "product", "uniref90", "cds", "genomes", "cds_with_go"]
                     + [f"go_{s}" for s in SOURCES if source] + ["go_terms", "go_ids"])
        for key, p in sorted(found.items(), key=lambda kv: (-len(kv[1]["genomes"]), -kv[1]["cds"], kv[0])):
            out.writerow([key, p["product"].most_common(1)[0][0],
                          p["uniref90"].most_common(1)[0][0] if p["uniref90"] else "",
                          p["cds"], len(p["genomes"]), p["with_go"]]
                         + [p["source"][s] for s in SOURCES if source]
                         + [len(p["go"]), ",".join(sorted(p["go"]))])
    print(f"{len(found)} proteins from {sum(p['cds'] for p in found.values())} CDS -> {args.out}")


if __name__ == "__main__":
    main()
