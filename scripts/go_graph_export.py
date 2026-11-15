"""Turn go_graph.py GraphML files into Gephi (.gexf) and Cytoscape (.cyjs) files with a layout.

GraphML holds no positions, so a 5,000-term graph opens as one blob. This lays each GO domain
out in its own column, one row per longest depth (roots on top), and orders each row by the
mean position of its parents so edges cross less. Nodes are coloured by domain and sized by
the genes that pass through them.

  python scripts/go_graph_export.py data/go_data/gaborieau/host/graph/hosts/DIJ07.graphml out_folder
  python scripts/go_graph_export.py data/go_data/gaborieau/host/graph/hosts/*.graphml out_folder

Writes NAME.gexf (Gephi: File > Open) and NAME.cyjs (Cytoscape: File > Import > Network from
File) for every input. Standard library only.
"""
import collections
import json
import math
import os
import sys
import xml.etree.ElementTree as ET
from xml.sax.saxutils import quoteattr

NS = "{http://graphml.graphdrawing.org/xmlns}"
COLOURS = {"P": (47, 109, 181), "F": (194, 90, 27), "C": (47, 138, 87)}
COLUMN = {"P": 0, "F": 1, "C": 2}
WIDTH, GAP, ROW = 4000, 600, 250


def read_graphml(path):
    root = ET.parse(path).getroot()
    nodes = {n.get("id"): {d.get("key"): d.text for d in n.findall(NS + "data")} for n in root.iter(NS + "node")}
    edges = [(e.get("source"), e.get("target"), e.find(NS + "data").text) for e in root.iter(NS + "edge")]
    return nodes, edges


def layout(nodes, edges):
    """{term: (x, y)}: a column per domain, a row per longest depth, rows ordered by parents."""
    parents = collections.defaultdict(list)
    for child, parent, _ in edges:
        parents[child].append(parent)
    rows = collections.defaultdict(list)
    for term, v in nodes.items():
        rows[v["domain"], int(v["longest_depth"])].append(term)
    order = {}
    for _ in range(2):  # a second pass lets rows settle on their parents' final order
        for key in sorted(rows, key=lambda k: k[1]):
            row = rows[key]
            row.sort(key=lambda t: (sum(order.get(p, 0.5) for p in parents[t]) / len(parents[t])
                                    if parents[t] else 0.5, -int(nodes[t]["genes_through"])))
            order.update({t: (i + 0.5) / len(row) for i, t in enumerate(row)})
    return {t: (COLUMN[v["domain"]] * (WIDTH + GAP) + order[t] * WIDTH, int(v["longest_depth"]) * ROW)
            for t, v in nodes.items()}


def size(v, most):
    return round(4 + 26 * math.sqrt(int(v["genes_through"]) / most), 2)


def write_gexf(path, nodes, edges, pos):
    most = max(int(v["genes_through"]) for v in nodes.values())
    keys = [k for k in next(iter(nodes.values())) if k != "name"]
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<gexf xmlns="http://www.gexf.net/1.2draft" xmlns:viz="http://www.gexf.net/1.2draft/viz" version="1.2">',
           '<graph defaultedgetype="directed"><attributes class="node">']
    out += [f'<attribute id="{k}" title="{k}" type="{"string" if k == "domain" else "integer"}"/>' for k in keys]
    out.append('</attributes><attributes class="edge"><attribute id="relation" title="relation" type="string"/>'
               '</attributes><nodes>')
    for term, v in nodes.items():
        r, g, b = COLOURS[v["domain"]]
        x, y = pos[term]
        values = "".join(f"<attvalue for={quoteattr(k)} value={quoteattr(v[k])}/>" for k in keys)
        out.append(f'<node id="{term}" label={quoteattr(term + " " + v["name"])}><attvalues>{values}</attvalues>'
                   f'<viz:color r="{r}" g="{g}" b="{b}"/><viz:position x="{x:.1f}" y="{-y:.1f}" z="0"/>'
                   f'<viz:size value="{size(v, most)}"/></node>')
    out.append("</nodes><edges>")
    out += [f'<edge id="{i}" source="{c}" target="{p}"><attvalues><attvalue for="relation" value="{r}"/>'
            f"</attvalues></edge>" for i, (c, p, r) in enumerate(edges)]
    out.append("</edges></graph></gexf>")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(out))


def write_cyjs(path, name, nodes, edges, pos):
    most = max(int(v["genes_through"]) for v in nodes.values())
    data = {"data": {"name": name}, "elements": {
        "nodes": [{"data": {"id": t, "name": t, "label": v["name"], "colour": "#%02x%02x%02x" % COLOURS[v["domain"]],
                            "size": size(v, most), **{k: (v[k] if k == "domain" else int(v[k]))
                                                      for k in v if k not in ("name", "domain")},
                            "domain": v["domain"]},
                   "position": {"x": round(pos[t][0], 1), "y": round(pos[t][1], 1)}} for t, v in nodes.items()],
        "edges": [{"data": {"id": f"e{i}", "source": c, "target": p, "relation": r}}
                  for i, (c, p, r) in enumerate(edges)]}}
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    *inputs, out = sys.argv[1:]
    os.makedirs(out, exist_ok=True)
    for path in inputs:
        name = os.path.splitext(os.path.basename(path))[0]
        nodes, edges = read_graphml(path)
        pos = layout(nodes, edges)
        write_gexf(os.path.join(out, name + ".gexf"), nodes, edges, pos)
        write_cyjs(os.path.join(out, name + ".cyjs"), name, nodes, edges, pos)
        print(f"{name}: {len(nodes)} terms, {len(edges)} edges")


if __name__ == "__main__":
    main()
