"""Run from anywhere: python tests/test_go_dag.py

Expected values were traced by hand from GO release 2022-03-22. The tests that need it look for
code/GO/go-basic.obo (DeepPBI-KG ships that release there) or data/go-basic-2022-03-22.obo:

  curl -L -o data/go-basic-2022-03-22.obo http://release.geneontology.org/2022-03-22/ontology/go-basic.obo

Without either file those tests are skipped; the rest run on inline data.
"""
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
OBO_2022 = next((p for p in (os.path.join(ROOT, "code", "GO", "go-basic.obo"),
                             os.path.join(ROOT, "data", "go-basic-2022-03-22.obo")) if os.path.exists(p)), None)
needs_obo_2022 = unittest.skipUnless(OBO_2022, "GO release 2022-03-22 not found; see the top of this file")
from go_dag import Dag, RELATION_SETS, parse_obo
from go_rollup import deepest, level_rows, measure, read_pairs, resolve_genes, slim_row
from go_check_quickgo import read_terms, rows, sample_terms

BP, MF, CC = "GO:0008150", "GO:0003674", "GO:0005575"
NUCLEUS, CYTOPLASM = "GO:0005634", "GO:0005737"
CELLULAR_ANATOMICAL_ENTITY, INTRACELLULAR_STRUCTURE = "GO:0110165", "GO:0005622"
DNA_REPAIR, REGULATION_OF_DNA_REPAIR = "GO:0006281", "GO:0006282"


@needs_obo_2022
class GoDagTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.terms, alt = parse_obo(OBO_2022)
        cls.dags = {name: Dag(cls.terms, alt, rel) for name, rel in RELATION_SETS.items()}
        cls.is_a = cls.dags["is_a"]
        cls.part_of = cls.dags["is_a+part_of"]
        cls.regulates = cls.dags["is_a+part_of+regulates"]

    def test_term_fields(self):
        self.assertEqual(self.terms[NUCLEUS]["name"], "nucleus")
        self.assertEqual(self.terms[NUCLEUS]["domain"], "C")
        self.assertEqual(self.terms[DNA_REPAIR]["domain"], "P")
        self.assertEqual(self.terms["GO:0016787"]["domain"], "F")

    def test_roots_are_the_three_domains(self):
        for dag in (self.part_of, self.regulates):
            self.assertEqual(dag.roots, {BP, MF, CC})

    def test_parents_follow_chosen_relations(self):
        self.assertEqual(self.is_a.parents[CYTOPLASM], {CELLULAR_ANATOMICAL_ENTITY})
        self.assertEqual(self.part_of.parents[CYTOPLASM], {CELLULAR_ANATOMICAL_ENTITY, INTRACELLULAR_STRUCTURE})
        self.assertEqual(self.part_of.parents[REGULATION_OF_DNA_REPAIR], {"GO:0051052", "GO:2001020"})
        self.assertEqual(self.regulates.parents[REGULATION_OF_DNA_REPAIR], {"GO:0051052", "GO:2001020", DNA_REPAIR})

    def test_ancestors_through_two_paths(self):
        # nucleus -> intracellular membrane-bounded organelle, whose two parents meet again at organelle.
        is_a = {"GO:0043231", "GO:0043227", "GO:0043229", "GO:0043226", CELLULAR_ANATOMICAL_ENTITY, CC}
        self.assertEqual(self.is_a.ancestors(NUCLEUS), is_a)
        # part_of adds one term, from intracellular organelle.
        self.assertEqual(self.part_of.ancestors(NUCLEUS), is_a | {INTRACELLULAR_STRUCTURE})

    def test_regulates_only_reaches_with_regulates(self):
        self.assertNotIn(DNA_REPAIR, self.part_of.ancestors(REGULATION_OF_DNA_REPAIR))
        self.assertIn(DNA_REPAIR, self.regulates.ancestors(REGULATION_OF_DNA_REPAIR))

    def test_depth(self):
        self.assertEqual(self.part_of.longest_depth(CC), 0)
        self.assertEqual(self.part_of.shortest_depth(CC), 0)
        # catalytic activity -> molecular_function; hydrolase -> catalytic activity.
        self.assertEqual(self.is_a.longest_depth("GO:0016787"), 2)
        self.assertEqual(self.is_a.shortest_depth("GO:0016787"), 2)
        # Cytoplasm is one path under is_a, and two different lengths once part_of is added.
        self.assertEqual(self.is_a.longest_depth(CYTOPLASM), 2)
        self.assertEqual(self.is_a.shortest_depth(CYTOPLASM), 2)
        self.assertEqual(self.part_of.longest_depth(CYTOPLASM), 3)
        self.assertEqual(self.part_of.shortest_depth(CYTOPLASM), 2)
        self.assertEqual(self.is_a.longest_depth(NUCLEUS), 5)

    def test_resolve(self):
        self.assertEqual(self.is_a.resolve(NUCLEUS), NUCLEUS)
        self.assertEqual(self.is_a.resolve("GO:0044763"), "GO:0009987")  # alt_id of cellular process
        self.assertEqual(self.is_a.resolve("GO:0000004"), BP)  # alt_id of a root
        self.assertIsNone(self.is_a.resolve("GO:0000005"))  # obsolete
        self.assertIsNone(self.is_a.resolve("GO:9999999"))  # unknown

    def test_obsolete_terms_are_left_out(self):
        self.assertIn("GO:0000005", self.terms)
        self.assertNotIn("GO:0000005", self.is_a.live)

    def test_leaf(self):
        self.assertTrue(self.regulates.is_leaf("GO:0090344"))
        self.assertFalse(self.is_a.is_leaf(NUCLEUS))

    def test_at_level_longest(self):
        # Nucleus is at depth 5, and branching gives two terms at level 3 and two at level 2.
        self.assertEqual(self.part_of.at_level(NUCLEUS, 5), {NUCLEUS})
        self.assertEqual(self.part_of.at_level(NUCLEUS, 9), {NUCLEUS})
        self.assertEqual(self.part_of.at_level(NUCLEUS, 4), {"GO:0043231"})
        self.assertEqual(self.part_of.at_level(NUCLEUS, 3), {"GO:0043227", "GO:0043229"})
        self.assertEqual(self.part_of.at_level(NUCLEUS, 2), {"GO:0043226", INTRACELLULAR_STRUCTURE})
        self.assertEqual(self.part_of.at_level(NUCLEUS, 1), {CELLULAR_ANATOMICAL_ENTITY})
        self.assertEqual(self.part_of.at_level(NUCLEUS, 0), {CC})

    def test_at_level_depends_on_depth_type(self):
        # Cytoplasm: longest depth 3 through part_of, shortest depth 2 through is_a.
        self.assertEqual(self.part_of.at_level(CYTOPLASM, 2), {INTRACELLULAR_STRUCTURE})
        self.assertEqual(self.part_of.at_level(CYTOPLASM, 2, depth="shortest"), {CYTOPLASM})
        self.assertEqual(self.part_of.at_level(CYTOPLASM, 1, depth="shortest"), {CELLULAR_ANATOMICAL_ENTITY})

    def test_at_level_always_finds_a_term(self):
        for dag in (self.is_a, self.part_of, self.regulates):
            for depth, measure in (("longest", dag.longest_depth), ("shortest", dag.shortest_depth)):
                for term in dag.live:
                    for level in (0, 1, 3, 5):
                        found = dag.at_level(term, level, depth)
                        self.assertTrue(found, (term, level, depth))
                        if measure(term) > level:
                            self.assertTrue(found <= dag.ancestors(term))
                            self.assertEqual({measure(a) for a in found}, {level})

    def test_whole_graph_is_consistent(self):
        for dag in self.dags.values():
            for term in dag.live:
                self.assertNotIn(term, dag.ancestors(term))
                self.assertLessEqual(dag.shortest_depth(term), dag.longest_depth(term))
                for parent in dag.parents[term]:
                    self.assertLess(dag.longest_depth(parent), dag.longest_depth(term))


@needs_obo_2022
class RollupTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        terms, alt = parse_obo(OBO_2022)
        cls.terms, cls.dag = terms, Dag(terms, alt, RELATION_SETS["is_a+part_of"])

    def test_read_pairs_table(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "pairs.tsv")
            with open(path, "w") as handle:
                handle.write("g1\tGO:0005634;GO:0043231\ng2\tGO:0005737,GO:0003674\ng3\t\ng1\tGO:0044763\n")
            self.assertEqual(read_pairs(path), {"g1": {"GO:0005634", "GO:0043231", "GO:0044763"},
                                                "g2": {"GO:0005737", "GO:0003674"}})

    def test_resolve_genes_splits_domains_and_counts_drops(self):
        by_domain, dropped = resolve_genes({"g": {"GO:0005634", "GO:0044763", "GO:0000005", "GO:9999999"}},
                                           self.terms, self.dag)
        self.assertEqual(by_domain["C"], {"g": {"GO:0005634"}})
        self.assertEqual(by_domain["P"], {"g": {"GO:0009987"}})  # alternate ID updated
        self.assertEqual(dropped, 2)

    def test_deepest_terms(self):
        self.assertEqual(deepest({NUCLEUS, "GO:0043231", CC, CYTOPLASM}, self.dag), {NUCLEUS, CYTOPLASM})

    def test_level_rows(self):
        # g1 reaches level 3 at two terms; g2 stops above it.
        rows = level_rows({"g1": {NUCLEUS, "GO:0043231"}, "g2": {CELLULAR_ANATOMICAL_ENTITY}}, self.dag, "longest", [3])
        self.assertEqual(rows[0]["genes"], 2)
        self.assertEqual(rows[0]["distinct_terms"], 3)
        self.assertEqual(rows[0]["median_terms_per_gene"], 1.5)
        self.assertEqual(rows[0]["pct_deepest_terms_above_level"], 50)


class SlimTest(unittest.TestCase):
    """root -> a -> b -> c, plus root -> d -> c, and x directly under root. Slim: a, b, d, root."""

    @classmethod
    def setUpClass(cls):
        records = [("GO:1", "", ""), ("GO:2", "GO:1", ""), ("GO:3", "GO:2", ""), ("GO:4", "GO:3", "GO:5"),
                   ("GO:5", "GO:1", ""), ("GO:6", "GO:1", "")]
        slim = {"GO:1", "GO:2", "GO:3", "GO:5"}
        text = ""
        for term, parent, other in records:
            text += f"[Term]\nid: {term}\nname: {term}\nnamespace: biological_process\n"
            text += "subset: test_slim\n" if term in slim else ""
            text += "".join(f"is_a: {p} ! parent\n" for p in (parent, other) if p)
            text += "\n"
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "mini.obo")
            with open(path, "w") as handle:
                handle.write(text)
            cls.terms, alt = parse_obo(path)
        cls.dag = Dag(cls.terms, alt, RELATION_SETS["is_a"])
        cls.slim = cls.dag.slim_terms(cls.terms, "test_slim")

    def test_slim_terms_leave_out_the_root(self):
        self.assertEqual(self.slim, {"GO:2", "GO:3", "GO:5"})
        with self.assertRaises(ValueError):
            self.dag.slim_terms(self.terms, "no_such_slim")

    def test_to_slim_keeps_only_the_most_specific(self):
        self.assertEqual(self.dag.to_slim("GO:3", self.slim), {"GO:3"})  # itself, not its slim parent GO:2
        self.assertEqual(self.dag.to_slim("GO:4", self.slim), {"GO:3", "GO:5"})  # two branches
        self.assertEqual(self.dag.to_slim("GO:6", self.slim), set())  # only the root is above it

    def test_slim_row(self):
        row = slim_row({"g1": {"GO:4", "GO:3"}, "g2": {"GO:6"}, "g3": {"GO:5"}}, self.dag, self.slim)
        self.assertEqual(row["genes"], 3)
        self.assertEqual(row["deepest_terms"], 3)  # GO:3 is above GO:4 in g1
        self.assertAlmostEqual(row["pct_deepest_terms_covered"], 200 / 3)
        self.assertEqual(row["distinct_slim_terms"], 2)
        self.assertEqual(row["median_slim_terms_per_gene"], 1)
        self.assertEqual(row["genes_only_other"], 1)


@needs_obo_2022
class QuickGoCheckTest(unittest.TestCase):
    """Offline: the sampling and comparison logic, without calling QuickGO."""

    @classmethod
    def setUpClass(cls):
        cls.terms, alt = parse_obo(OBO_2022)
        cls.dag = Dag(cls.terms, alt, RELATION_SETS["is_a+part_of"])

    def test_sample_is_repeatable_and_even_across_domains(self):
        first = sample_terms(self.dag, self.terms, 30, seed=0)
        self.assertEqual(first, sample_terms(self.dag, self.terms, 30, seed=0))
        self.assertNotEqual(first, sample_terms(self.dag, self.terms, 30, seed=1))
        for domain in "FPC":
            self.assertEqual(sum(self.terms[t]["domain"] == domain for t in first), 10)
        self.assertFalse(set(first) & self.dag.roots)

    def test_rows_report_each_side_of_a_difference(self):
        true_set = set(self.dag.ancestors(NUCLEUS))
        same, differ = rows([NUCLEUS], self.terms, self.dag, {NUCLEUS: true_set}), None
        self.assertEqual(same[0]["match"], "yes")
        changed = (true_set - {CC}) | {"GO:9999999"}
        differ = rows([NUCLEUS], self.terms, self.dag, {NUCLEUS: changed})[0]
        self.assertEqual(differ["match"], "no")
        self.assertEqual(differ["only_ours"], CC)
        self.assertEqual(differ["only_quickgo"], "GO:9999999")
        self.assertEqual(rows([NUCLEUS], self.terms, self.dag, {})[0]["match"], "missing")

    def test_read_terms_from_any_file(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as handle:
            handle.write(f"gene\tgo\nx\t{NUCLEUS},{NUCLEUS}\ny\tGO:9999999 {CC}\n")
        try:
            self.assertEqual(read_terms(handle.name, self.dag), sorted({NUCLEUS, CC}))
        finally:
            os.unlink(handle.name)


class MeasureTest(unittest.TestCase):
    """R <- A <- B <- D (is_a), and F directly under R."""

    @classmethod
    def setUpClass(cls):
        obo = "format-version: 1.2\n" + "".join(
            f"\n[Term]\nid: GO:{t}\nnamespace: biological_process\n" + (f"is_a: GO:{p}\n" if p else "")
            for t, p in [("R", ""), ("A", "R"), ("B", "A"), ("D", "B"), ("F", "R")])
        with tempfile.NamedTemporaryFile("w", suffix=".obo", delete=False) as handle:
            handle.write(obo)
        try:
            terms, alt = parse_obo(handle.name)
        finally:
            os.unlink(handle.name)
        cls.dag = Dag(terms, alt, RELATION_SETS["is_a"])

    def test_term_with_its_ancestor_is_a_branch(self):
        m = measure([{"GO:D", "GO:B"}], self.dag)
        self.assertEqual((m["proteins"], m["pct_with_term_and_ancestor"]), (1, 100))
        self.assertEqual(m["pct_annotations_with_listed_ancestor"], 50)
        self.assertEqual(m["pct_annotations_that_are_listed_ancestors"], 50)
        self.assertEqual(m["mean_pct_branch_listed"], 50)  # Tip D: branch {B, A}, B listed.
        self.assertEqual((m["pct_annotations_on_leaf"], m["median_depth"]), (50, 2.5))

    def test_unrelated_terms_are_independent(self):
        m = measure([{"GO:D", "GO:F"}], self.dag)  # F's branch is only the root, so it is skipped.
        self.assertEqual((m["pct_only_independent_multi"], m["pct_tips_empty_branch"]), (100, 100))


if __name__ == "__main__":
    unittest.main()
