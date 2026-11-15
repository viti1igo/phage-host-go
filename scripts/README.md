# Script reference

The pipeline and the commands that made the results are in the [main README](../README.md). This page gives each script's inputs, options and outputs. Every script uses only the Python standard library, except that FASTA input to `go_extract.py` needs [DIAMOND](https://github.com/bbuchfink/diamond) (tested with 2.2.8). Run commands from the repo root.

## `go_extract.py`

Turns proteins, or another tool's GO output, into cleaned GO terms and summaries.

| Input | Command |
|---|---|
| Bakta (one `.tsv`, or a folder of `<genome>/<genome>.tsv` runs) | `python scripts/go_extract.py run data/bakta_hosts out --format bakta --genomes interaction_matrix.csv --obo data/go-basic.obo`. `--genomes` keeps only the genomes in the first column of a table with a header row; gene IDs become `genome\|locus tag` |
| Protein FASTA | `python scripts/go_extract.py run proteins.faa out --db data/sprot_db --sprot data/uniprot_sprot.dat.gz --obo data/go-basic.obo` |
| eggNOG-mapper | `python scripts/go_extract.py run sample.emapper.annotations out --format emapper --obo data/go-basic.obo` |
| InterProScan TSV | `python scripts/go_extract.py run sample.tsv out --format interproscan --obo data/go-basic.obo` |
| Plain table | `python scripts/go_extract.py run pairs.tsv out --format table --obo data/go-basic.obo` |

A plain table has one gene per line: the gene ID, a tab, then GO IDs separated by commas or semicolons. A gene can also appear on several lines.

**FASTA input.** Build the DIAMOND database once, from the same Swiss-Prot file that is later used to look up GO, so the two can't come from different releases (about 40 seconds):

```bash
curl -L -o data/uniprot_sprot.dat.gz https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
python scripts/go_extract.py build-db data/uniprot_sprot.dat.gz data/sprot_db
```

DIAMOND finds up to 5 Swiss-Prot hits per protein. The script keeps the hit with the highest bit score among those that pass every threshold, and copies that entry's GO terms. Swiss-Prot has few phage proteins, so only about 5% of phage proteins get GO this way.

| Option | Default | Meaning |
|---|---|---|
| `--min-identity` | 30 | FASTA: minimum % identical amino acids in the aligned region |
| `--min-coverage` | 50 | FASTA: minimum % of the query protein covered by the alignment |
| `--max-evalue` | 1e-5 | FASTA: maximum E-value |
| `--threads` | 4 | DIAMOND threads |
| `--slim` | `goslim_metagenomics` | GO slim to summarise with; any subset named in the OBO file works |
| `--level` | 3 | Level for the secondary summary, counted from the domain root |

**Outputs**

| File | Contents |
|---|---|
| `pairs.tsv` | One line per gene and GO ID, as found, before any clean-up |
| `hits.tsv`, `diamond_hits.tsv` | FASTA input only: the hit kept for each protein, and all DIAMOND hits before the thresholds |
| `provenance.json` | Command, date, input and its SHA-256, GO release, database and its SHA-256, thresholds, DIAMOND version, counts, and the pairs dropped because their GO ID is retired or unknown in this GO release |
| `gene_terms.tsv` | Every live term per gene, its depths, whether it is a deepest term and whether it is a GO leaf |
| `rollup_goslim_metagenomics.tsv` | Each deepest term with the most specific GO slim terms above it, or `other` |
| `slim_summary.csv` | Per GO domain: genes, deepest terms, % covered by the slim, distinct slim terms, slim terms per gene |
| `rollup_level_3.tsv` | Each deepest term with its ancestors at the chosen level |
| `levels.csv` | For every level: distinct terms, terms per gene, and how many deepest terms are already above it |
| `measurements.csv` | Whether genes list whole branches or independent terms, for three sets of GO edges |

## `go_rollup.py`

The summaries above, from any gene-to-GO table. `go_extract.py run` calls it.

```bash
python scripts/go_rollup.py data/go-basic.obo pairs.tsv out --slim goslim_metagenomics --level 3
```

Options: `--depth longest|shortest` (how levels are counted), `--relations is_a|is_a+part_of|is_a+part_of+regulates`, and `--genes NAME=CSV` (repeatable) to summarise a group of genes listed in a CSV with a `gene` column.

## `phage_go_fill.py`

```bash
python scripts/phage_go_fill.py phages.tsv out_folder
```

Reads one Bakta `.tsv` covering many phages and writes `out_folder/<phage>/<phage>.tsv`, the layout `go_extract.py` and `go_graph.py --bakta` read. A CDS that Bakta gave GO keeps it. A CDS without GO but with a UniRef90 cluster gets the GO terms of the cluster's representative UniProtKB entry, looked up on `rest.uniprot.org` 400 at a time. Representatives that are UniParc entries (`UPI…`) have no GO. Also writes `go_source.tsv` (per CDS: `bakta`, `uniref90` or `none`) and `uniref90_go_cache.tsv`, which later runs reuse instead of asking UniProt again.

## `go_graph.py`

```bash
python scripts/go_graph.py data/go-basic.obo data/bakta_hosts out --bakta --genomes interaction_matrix.csv              # writes NAME=host
python scripts/go_graph.py data/go-basic.obo results/gaborieau/phage/by_phage out --bakta --name phage --genomes matrix_phages.csv
```

With `--bakta`, a gene is the Bakta gene name, else its UniRef90 cluster, else the feature's own ID; copies across genomes merge. For each group NAME it writes `genes_NAME.tsv`, `subtrees_NAME.tsv`, `nodes_NAME.tsv`, `edges_NAME.tsv` and `graph_NAME.graphml`, plus one graph per genome in `NAMEs/` and a summary `NAMEs.tsv`. `--relations` picks the GO edges (default `is_a+part_of`).

Without `--bakta`, the input is a gene-to-GO table and `--group NAME=CSV` (repeatable) names the genes of each graph; labels that differ only in a `_N` suffix (`gene=dnaA_2`) merge into one gene.

The `quickgo_chart` column of `genes_NAME.tsv` links to QuickGO's chart of the gene's deepest terms. QuickGO uses its own GO release and also draws `regulates` edges, so it can show more than the subtree here.

## `go_graph_export.py`

```bash
python scripts/go_graph_export.py graph_host.graphml [more.graphml ...] out_folder
```

Writes `NAME.gexf` (Gephi) and `NAME.cyjs` (Cytoscape) for every input. Layout: one column per GO domain, one row per longest depth (roots on top), each row ordered by its parents' positions so edges cross less. Colour by domain, size by genes through the term.

## `go_check_quickgo.py`

```bash
python scripts/go_check_quickgo.py data/go-basic.obo out.csv --terms deepest_terms.txt     # every GO ID in a file
python scripts/go_check_quickgo.py data/go-basic.obo out.csv --sample 300 --seed 0          # a repeatable random sample
```

Compares our ancestors of each term with QuickGO's, 20 terms per request with `--pause` seconds (default 1) between requests, and prints both GO releases. `go-basic.obo` leaves out edges between GO domains that QuickGO keeps, so some `part_of` differences are expected. Needs internet access to `www.ebi.ac.uk`.

## `go_dag.py`

The library behind every script:

```python
from go_dag import Dag, RELATION_SETS, parse_obo
terms, alt = parse_obo("data/go-basic.obo")
dag = Dag(terms, alt, RELATION_SETS["is_a+part_of"])
dag.ancestors(dag.resolve("GO:0044763"))
```

`Dag` gives parents, ancestors, longest and shortest depth, the terms at a level, GO slim terms, and resolves alternate IDs; obsolete terms are left out.

## Limits

- GO from FASTA is copied from the most similar Swiss-Prot protein, so it can be wrong where proteins are similar but differ in function.
- Results depend on the GO and UniProt releases, so always keep `provenance.json` with the results.
- GO slims don't fit every organism. `goslim_metagenomics` covers only about half of phage cellular component terms.
