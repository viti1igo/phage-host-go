"""Split one Bakta run over many phages into one folder per phage, filling GO from UniRef90.

Bakta GO comes first. A CDS with no GO from Bakta gets the GO terms of its UniRef90 cluster's
representative UniProtKB entry (looked up on rest.uniprot.org, cached in CACHE). UniRef50 is
not used. Writes OUT/<phage>/<phage>.tsv (the layout go_extract.py and go_graph.py --bakta read)
and OUT/go_source.tsv saying where each CDS's GO came from.

  python scripts/phage_go_fill.py data/go_data/gaborieau/phage/bakta/all/phages.tsv data/go_data/gaborieau/phage/bakta/by_phage
"""
import csv
import os
import re
import sys
import urllib.parse
import urllib.request

GO_ID = re.compile(r"GO:\d{7}")
UR90 = re.compile(r"UniRef90_(\w+)")
ACC = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")  # UniProtKB accession
URL = "https://rest.uniprot.org/uniprotkb/accessions?fields=accession,go_id&format=tsv&accessions="


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
    for i in range(0, len(todo), 400):
        batch = todo[i:i + 400]
        with urllib.request.urlopen(URL + urllib.parse.quote(",".join(batch)), timeout=120) as reply:
            rows = list(csv.DictReader(reply.read().decode().splitlines(), delimiter="\t"))
        found = {r["Entry"]: GO_ID.findall(r["Gene Ontology IDs"]) for r in rows}
        known.update({a: found.get(a, []) for a in batch})  # missing = obsolete or UniParc-only: no GO
        print(f"UniProt: {min(i + 400, len(todo))}/{len(todo)} accessions", file=sys.stderr)
        save(cache, known)
    save(cache, known)
    return known


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    tsv, out = sys.argv[1:]
    with open(tsv) as handle:
        lines = handle.read().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("#Sequence Id"))
    header, body = lines[:start + 1], [l.split("\t") for l in lines[start + 1:] if l]
    need = [m.group(1) for f in body if f[1] == "cds" and not GO_ID.search(f[8])
            for m in [UR90.search(f[8])] if m]
    os.makedirs(out, exist_ok=True)
    go = fetch_go(need, os.path.join(out, "uniref90_go_cache.tsv"))
    by_phage, source = {}, []
    for f in body:
        if f[1] == "cds":
            m = UR90.search(f[8])
            if GO_ID.search(f[8]):
                how = "bakta"
            elif m and go.get(m.group(1)):
                f[8] = ", ".join(filter(None, [f[8]] + go[m.group(1)]))
                how = "uniref90"
            else:
                how = "none"
            source.append((f[0], f[5], how, m.group(1) if m else ""))
        by_phage.setdefault(f[0], []).append("\t".join(f))
    for phage, rows in by_phage.items():
        os.makedirs(os.path.join(out, phage), exist_ok=True)
        with open(os.path.join(out, phage, phage + ".tsv"), "w") as handle:
            handle.write("\n".join(header + rows) + "\n")
    with open(os.path.join(out, "go_source.tsv"), "w") as handle:
        handle.write("phage\tlocus_tag\tgo_source\tuniref90_rep\n")
        handle.writelines("\t".join(s) + "\n" for s in source)
    counts = {k: sum(s[2] == k for s in source) for k in ("bakta", "uniref90", "none")}
    print(f"{len(by_phage)} phages, {len(source)} CDS: {counts}")


if __name__ == "__main__":
    main()
