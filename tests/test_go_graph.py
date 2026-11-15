"""Run from anywhere: python tests/test_go_graph.py"""
import collections
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from go_dag import Dag, RELATION_SETS, parse_obo
from go_graph import graph, merge, read_bakta_genes, subtree, unique_gene

OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000002
name: a
namespace: biological_process
is_a: GO:0000001

[Term]
id: GO:0000003
name: b
namespace: biological_process
is_a: GO:0000002
relationship: part_of GO:0000004

[Term]
id: GO:0000004
name: c
namespace: biological_process
is_a: GO:0000001

[Term]
id: GO:0000005
name: old
namespace: biological_process
is_obsolete: true
"""


class GraphTest(unittest.TestCase):
    def setUp(self):
        with tempfile.NamedTemporaryFile("w", suffix=".obo", delete=False) as handle:
            handle.write(OBO)
        self.terms, alt = parse_obo(handle.name)
        os.unlink(handle.name)
        self.relations = RELATION_SETS["is_a+part_of"]
        self.dag = Dag(self.terms, alt, self.relations)

    def test_merge_and_graph(self):
        self.assertEqual(unique_gene("gene=dnaA_12"), "gene=dnaA")
        self.assertEqual(unique_gene("protein=putative_2"), "protein=putative_2")
        pairs = {"gene=x_1": {"GO:0000003"}, "gene=x_2": {"GO:0000002", "GO:0000005"},
                 "gene=y": {"GO:0000004"}, "gene=z": {"GO:0000005"}, "gene=other": {"GO:0000002"}}
        genes = merge(pairs, {"gene=x_1", "gene=x_2", "gene=y", "gene=z"}, self.dag)
        self.assertEqual(set(genes), {"gene=x", "gene=y"})  # z has only an obsolete term
        self.assertEqual(genes["gene=x"], ({"gene=x_1", "gene=x_2"}, {"GO:0000002", "GO:0000003"}))

        nodes, edges = graph(genes, self.dag, self.terms, self.relations)
        for edge in edges:  # every edge is in the OBO file
            self.assertIn((edge["relation"], edge["parent"]), self.terms[edge["child"]]["edges"])
        self.assertEqual(len(edges), 4)
        through = collections.Counter(t for _, listed in genes.values() for t in subtree(listed, self.dag))
        self.assertEqual({n["go_id"]: n["genes_through"] for n in nodes}, dict(through))
        self.assertEqual({n["go_id"]: n["genes_deepest"] for n in nodes if n["genes_deepest"]},
                         {"GO:0000003": 1, "GO:0000004": 1})

    def test_bakta_genes_merge_across_genomes(self):
        header = "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
        rows = {"S1": ["c\tcds\t1\t9\t+\tA_1\tdnaA\tx\tGO:0000003, UniRef:UniRef90_P1",
                       "c\tcds\t1\t9\t+\tA_2\t\tx\tGO:0000004, UniRef:UniRef90_Q9",
                       "c\tcds\t1\t9\t+\tA_3\t\tx\tGO:0000002",
                       "c\tcds\t1\t9\t+\tA_4\tnoGO\tx\tUniRef:UniRef90_P2"],
                "S2": ["c\tcds\t1\t9\t+\tB_1\tdnaA\tx\tGO:0000002"],
                "S3": ["c\tcds\t1\t9\t+\tC_1\tdnaA\tx\tGO:0000004"]}
        with tempfile.TemporaryDirectory() as folder:
            for genome, lines in rows.items():
                os.makedirs(os.path.join(folder, genome))
                with open(os.path.join(folder, genome, genome + ".tsv"), "w") as handle:
                    handle.write(header + "\n".join(lines) + "\n")
            pairs, names = read_bakta_genes(folder, {"S1", "S2"})  # S3 is not asked for
        self.assertEqual(names, {"S1|A_1": "dnaA", "S1|A_2": "UniRef90_Q9", "S1|A_3": "S1|A_3", "S2|B_1": "dnaA"})
        genes = merge(pairs, set(pairs), self.dag, key=names.get)
        self.assertEqual(genes["dnaA"], ({"S1|A_1", "S2|B_1"}, {"GO:0000002", "GO:0000003"}))
        self.assertEqual(set(genes), {"dnaA", "UniRef90_Q9", "S1|A_3"})

    def test_export_layout_and_files(self):
        import go_graph_export as ex
        from go_graph import write_graphml
        genes = {"g": (set(), {"GO:0000003"})}
        nodes, edges = graph(genes, self.dag, self.terms, self.relations)
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "h.graphml")
            write_graphml(path, nodes, edges)
            read_nodes, read_edges = ex.read_graphml(path)
            pos = ex.layout(read_nodes, read_edges)
            self.assertEqual(pos["GO:0000001"][1], 0)  # the root is on the top row
            self.assertTrue(all(pos[c][1] > pos[p][1] for c, p, _ in read_edges))  # children below parents
            ex.write_gexf(os.path.join(folder, "h.gexf"), read_nodes, read_edges, pos)
            ex.write_cyjs(os.path.join(folder, "h.cyjs"), "h", read_nodes, read_edges, pos)
            import json, xml.etree.ElementTree as ET
            self.assertEqual(len(ET.parse(os.path.join(folder, "h.gexf")).getroot().iter().__next__()), 1)
            self.assertEqual(len(json.load(open(os.path.join(folder, "h.cyjs")))["elements"]["edges"]), len(edges))


if __name__ == "__main__":
    unittest.main()
