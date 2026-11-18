"""Fill missing GO terms in Bakta tables from UniRef90, writing one table per genome.

Bakta GO comes first. A CDS with no GO from Bakta but with a UniRef90 cluster gets the GO terms
of the cluster's representative UniProtKB entry (looked up on rest.uniprot.org and cached).
Representatives that are UniParc entries (UPI...) have no GO. UniRef50 is not used.

The input is either one Bakta .tsv covering many genomes, split by its Sequence Id column (one
run over all phages, each phage one sequence), or a folder of Bakta runs, <genome>/<genome>.tsv
(one run per host). Writes to OUT:
  <genome>/<genome>.tsv   the Bakta table with GO added (the layout go_extract.py and
                          go_graph.py --bakta read)
  go_source.tsv           per CDS: genome, locus tag, go_source (bakta, uniref90,
                          uniref90_unisave or none), uniref90_rep
  uniref90_go_cache.tsv   UniProt answers, reused by later runs
  unisave_go_cache.tsv    with --unisave: GO of each representative's last version in UniSave

Many representatives in the Bakta database (v6.0, 2025) have since been deleted from UniProtKB,
so the current lookup finds no GO for them. --unisave takes the GO their last version had, from
UniSave; that GO is as old as the deletion, so go_source marks it uniref90_unisave.

  python scripts/uniref90_go_fill.py data/phage/bakta/phages.tsv results/gaborieau/phage/by_phage
  python scripts/uniref90_go_fill.py data/bakta_hosts results/gaborieau/host/bakta_filled --genomes interaction_matrix.csv

Standard library only; needs internet access for accessions not yet in the cache.
"""
import argparse
import concurrent.futures
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from go_extract import GO_ID, bakta_tables, read_genomes

UR90 = re.compile(r"UniRef90_(\w+)")
ACC = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")  # UniProtKB accession
URL = "https://rest.uniprot.org/uniprotkb/accessions?fields=accession,go_id&format=tsv&accessions="
BATCH = 400


def save(cache, known):
    with open(cache, "w") as handle:
        handle.writelines(f"{a}\t{','.join(g)}\n" for a, g in sorted(known.items()))


def fetch_go(accessions, cache):
    """{accession: [GO IDs]}, reading and extending the cache TSV."""
    known = {}
    if os.path.exists(cache):
        with open(cache) as handle:
            known = {a: g.split(",") if g else [] for a, g in csv.reader(handle, delimiter="\t")}
    todo = sorted(set(accessions) - set(known))
    known.update({a: [] for a in todo if not ACC.fullmatch(a)})  # UniParc (UPI...) representatives: no GO
    todo = [a for a in todo if ACC.fullmatch(a)]
    for i in range(0, len(todo), BATCH):
        batch = todo[i:i + BATCH]
        with urllib.request.urlopen(URL + urllib.parse.quote(",".join(batch)), timeout=120) as reply:
            rows = list(csv.DictReader(reply.read().decode().splitlines(), delimiter="\t"))
        found = {r["Entry"]: GO_ID.findall(r["Gene Ontology IDs"]) for r in rows}
        known.update({a: found.get(a, []) for a in batch})  # missing = obsolete or UniParc-only: no GO
        print(f"UniProt: {min(i + BATCH, len(todo))}/{len(todo)} accessions", file=sys.stderr)
        save(cache, known)
    save(cache, known)
    return known


def unisave_go(accession):
    """GO IDs in the last version of an entry, from UniSave, which keeps deleted entries too."""
    base = f"https://rest.uniprot.org/unisave/{accession}"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(base + "?format=json", timeout=60) as reply:
                versions = json.load(reply)["results"]
            if not versions:
                return []
            with urllib.request.urlopen(f"{base}?format=txt&versions={versions[0]['entryVersion']}", timeout=60) as reply:
                return sorted(set(re.findall(r"^DR   GO; (GO:\d{7})", reply.read().decode(), flags=re.M)))
        except (urllib.error.URLError, OSError, ValueError):
            time.sleep(5 * 2 ** attempt)
    return None  # failed: left out of the cache, so the next run asks again


def fetch_unisave(accessions, cache, workers=8):
    """{accession: [GO IDs]} from UniSave for the given UniProtKB accessions, cached in a TSV."""
    known = {}
    if os.path.exists(cache):
        with open(cache) as handle:
            known = {a: g.split(",") if g else [] for a, g in csv.reader(handle, delimiter="\t")}
    todo = sorted(set(accessions) - set(known))
    with concurrent.futures.ThreadPoolExecutor(workers) as pool:
        for i, (accession, found) in enumerate(zip(todo, pool.map(unisave_go, todo)), 1):
            if found is not None:
                known[accession] = found
            if i % 500 == 0 or i == len(todo):
                print(f"UniSave: {i}/{len(todo)} accessions", file=sys.stderr)
                save(cache, known)
    save(cache, known)
    return known


def read_table(path):
    """(header lines up to #Sequence Id, rows split on tabs)."""
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("#Sequence Id"))
    return lines[:start + 1], [line.split("\t") for line in lines[start + 1:] if line]


def genomes_of(path, genomes):
    """Yield (genome, header, rows): split by sequence for one file, one per run for a folder."""
    if os.path.isfile(path):
        header, rows = read_table(path)
        by_genome = {}
        for row in rows:
            by_genome.setdefault(row[0], []).append(row)
        for genome, own in by_genome.items():
            if genomes is None or genome in genomes:
                yield genome, header, own
    else:
        for genome, table in bakta_tables(path, genomes):
            yield (genome, *read_table(table))


def missing(row):
    """The UniRef90 representative of a CDS without GO, else None."""
    if row[1] != "cds" or GO_ID.search(row[8]):
        return None
    match = UR90.search(row[8])
    return match.group(1) if match else None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="one Bakta .tsv over many genomes, or a folder of <genome>/<genome>.tsv runs")
    parser.add_argument("out")
    parser.add_argument("--genomes", help="only the genomes in the first column of this table (header row skipped)")
    parser.add_argument("--unisave", action="store_true",
                        help="for representatives still without GO (mostly entries UniProt has deleted since the "
                             "Bakta database was built), take the GO of their last version from UniSave; "
                             "slow (one entry at a time, 8 in parallel), cached in unisave_go_cache.tsv")
    args = parser.parse_args()
    genomes = read_genomes(args.genomes) if args.genomes else None
    os.makedirs(args.out, exist_ok=True)

    need = {rep for _, _, rows in genomes_of(args.input, genomes) for rep in map(missing, rows) if rep}
    print(f"{len(need)} UniRef90 representatives to look up", file=sys.stderr)
    go = fetch_go(need, os.path.join(args.out, "uniref90_go_cache.tsv"))
    old = {}
    if args.unisave:
        old = fetch_unisave([a for a in need if ACC.fullmatch(a) and not go[a]],
                            os.path.join(args.out, "unisave_go_cache.tsv"))

    counts, written = {"bakta": 0, "uniref90": 0, "uniref90_unisave": 0, "none": 0}, 0
    with open(os.path.join(args.out, "go_source.tsv"), "w") as source:
        source.write("genome\tlocus_tag\tgo_source\tuniref90_rep\n")
        for genome, header, rows in genomes_of(args.input, genomes):
            for row in rows:
                if row[1] != "cds":
                    continue
                rep = UR90.search(row[8])
                if GO_ID.search(row[8]):
                    how = "bakta"
                elif rep and go.get(rep.group(1)):
                    row[8] = ", ".join(filter(None, [row[8]] + go[rep.group(1)]))
                    how = "uniref90"
                elif rep and old.get(rep.group(1)):
                    row[8] = ", ".join(filter(None, [row[8]] + old[rep.group(1)]))
                    how = "uniref90_unisave"
                else:
                    how = "none"
                counts[how] += 1
                source.write(f"{genome}\t{row[5]}\t{how}\t{rep.group(1) if rep else ''}\n")
            os.makedirs(os.path.join(args.out, genome), exist_ok=True)
            with open(os.path.join(args.out, genome, genome + ".tsv"), "w") as handle:
                handle.write("\n".join(header + ["\t".join(r) for r in rows]) + "\n")
            written += 1
    print(f"{written} genomes, {sum(counts.values())} CDS: {counts}")


if __name__ == "__main__":
    main()
