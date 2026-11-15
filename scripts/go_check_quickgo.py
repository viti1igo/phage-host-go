"""Check go_dag.py's ancestors against QuickGO on a reproducible random sample.

  python scripts/go_check_quickgo.py go-basic.obo quickgo_check.csv
  python scripts/go_check_quickgo.py go-basic.obo out.csv --sample 300 --seed 0
  python scripts/go_check_quickgo.py go-basic.obo out.csv --terms terms.txt

--terms checks every GO ID found in a file instead of a random sample, e.g. the deepest terms
of a gene set. Requests go 20 terms at a time with a pause between them (--pause seconds).

The same OBO file and seed always pick the same terms, so two runs are comparable.
The sample is split evenly across F, P and C. Standard library only.

QuickGO serves its own GO release, so a difference can be a release difference rather
than a bug: the script prints both releases and lists every differing term.
`go-basic.obo` also leaves out edges that cross GO domains, which QuickGO keeps, so
some `part_of` differences are expected by design.
"""
import argparse
import csv
import json
import random
import re
import time
import urllib.error
import urllib.request

from go_dag import Dag, RELATION_SETS, obo_version, parse_obo

API = "https://www.ebi.ac.uk/QuickGO/services/ontology/go"
QUICKGO_RELATIONS = {"is_a": "is_a", "is_a+part_of": "is_a,part_of"}
BATCH = 20


def sample_terms(dag, terms, count, seed):
    """Count terms, split evenly across the three domains, chosen by seed."""
    rng = random.Random(seed)
    chosen = []
    for domain in "FPC":
        pool = sorted(t for t in dag.live - dag.roots if terms[t]["domain"] == domain)
        chosen += rng.sample(pool, min(count // 3, len(pool)))
    return sorted(chosen)


def fetch(path):
    with urllib.request.urlopen(API + path, timeout=60) as response:
        return json.load(response)


def read_terms(path, dag):
    """The live GO IDs found anywhere in a file, sorted; alternate IDs become current ones."""
    with open(path, encoding="utf-8") as handle:
        found = {dag.resolve(t) for t in re.findall(r"GO:\d{7}", handle.read())}
    return sorted(t for t in found if t)


def quickgo_ancestors(ids, relations, pause=0.0):
    """{term: set of ancestors} from QuickGO, without the term itself; terms it lacks are left out."""
    found = {}
    for start in range(0, len(ids), BATCH):
        batch = ids[start:start + BATCH]
        data = fetch(f"/terms/{','.join(batch)}/ancestors?relations={relations}")
        for result in data["results"]:
            found[result["id"]] = set(result.get("ancestors", [])) - {result["id"]}
        time.sleep(pause)
    return found


def rows(chosen, terms, dag, theirs):
    """One row per term: our ancestors against QuickGO's."""
    out = []
    for term in chosen:
        ours = dag.ancestors(term)
        if term not in theirs:
            out.append({"go_id": term, "domain": terms[term]["domain"], "name": terms[term]["name"],
                        "ours": len(ours), "quickgo": "", "match": "missing", "only_ours": "", "only_quickgo": ""})
            continue
        their = theirs[term]
        out.append({"go_id": term, "domain": terms[term]["domain"], "name": terms[term]["name"],
                    "ours": len(ours), "quickgo": len(their), "match": "yes" if ours == their else "no",
                    "only_ours": " ".join(sorted(ours - their)),
                    "only_quickgo": " ".join(sorted(their - ours))})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("obo")
    parser.add_argument("out")
    parser.add_argument("--sample", type=int, default=300)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--relations", choices=list(QUICKGO_RELATIONS), default="is_a+part_of")
    parser.add_argument("--terms", help="check every GO ID in this file instead of a random sample")
    parser.add_argument("--pause", type=float, default=1.0, help="seconds between requests")
    args = parser.parse_args()

    terms, alt = parse_obo(args.obo)
    dag = Dag(terms, alt, RELATION_SETS[args.relations])
    chosen = read_terms(args.terms, dag) if args.terms else sample_terms(dag, terms, args.sample, args.seed)
    if not chosen:
        raise SystemExit("no live GO IDs to check")
    try:
        theirs = quickgo_ancestors(chosen, QUICKGO_RELATIONS[args.relations], args.pause)
        about = fetch("/about")
    except urllib.error.URLError as error:
        raise SystemExit(f"QuickGO request failed, so nothing was checked: {error}")

    table = rows(chosen, terms, dag, theirs)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)

    differ = [r for r in table if r["match"] == "no"]
    missing = [r["go_id"] for r in table if r["match"] == "missing"]
    print(f"our GO release: {obo_version(args.obo)}")
    print(f"QuickGO release: {about.get('go', {}).get('version', 'unknown')}")
    source = f"terms from {args.terms}" if args.terms else f"seed {args.seed}"
    print(f"{len(table)} terms checked with {args.relations}, {source}: "
          f"{len(table) - len(differ) - len(missing)} match, {len(differ)} differ, {len(missing)} not in QuickGO")
    if missing:
        print("  not in QuickGO:", " ".join(missing))
    for row in differ:
        print(f"  {row['go_id']} {row['name']} ({row['domain']}): "
              f"only ours [{row['only_ours']}] only QuickGO [{row['only_quickgo']}]")
    print("written to", args.out)


if __name__ == "__main__":
    main()
