#!/usr/bin/env bash

set -euo pipefail

WORK_DIR="${WORK_DIR:-$HOME/hare_work}"
TREE_SOURCE="${TREE_SOURCE:-$WORK_DIR/subset_parallel/pipeline_bulletproof_final/astral_species_tree_with_all_samples_grouped_modern_ancient_separate.tre}"
VCF="${VCF:-$WORK_DIR/MyHare_with_all_samples.vcf.gz}"
DSUITE="${DSUITE:-$HOME/Dsuite/Build/Dsuite}"
SETS_DIR="${SETS_DIR:-$WORK_DIR/for_data}"
OUT_DIR="${OUT_DIR:-$WORK_DIR/dsuite_results/fbranch_new}"

AMERICANUS_TREE="$OUT_DIR/astral_new_americanus_outgroup.tre"
RABBIT_TREE="$OUT_DIR/astral_new_rabbit_outgroup.tre"

if [[ ! -s "$TREE_SOURCE" ]]; then
    echo "Error: ASTRAL tree not found: $TREE_SOURCE" >&2
    exit 1
fi
if [[ ! -s "$VCF" ]]; then
    echo "Error: VCF not found: $VCF" >&2
    exit 1
fi
if [[ ! -x "$DSUITE" ]]; then
    echo "Error: Dsuite executable not found: $DSUITE" >&2
    exit 1
fi
if [[ ! -s "$SETS_DIR/sets_hare.txt" || ! -s "$SETS_DIR/sets_krol.txt" ]]; then
    echo "Error: sets_hare.txt and sets_krol.txt are required in $SETS_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

make_tree() {
    local source_outgroup_tip=$1
    local sets_file=$2
    local output_tree=$3

    python3 - "$TREE_SOURCE" "$sets_file" "$source_outgroup_tip" "$output_tree" <<'PY'
import sys


class Node:
    def __init__(self, name=""):
        self.name = name
        self.children = []


def parse_tree(text):
    pos = 0

    def skip_space():
        nonlocal pos
        while pos < len(text) and text[pos].isspace():
            pos += 1

    def skip_branch_length():
        nonlocal pos
        skip_space()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            while pos < len(text) and text[pos] not in ",();":
                pos += 1

    def read_label():
        nonlocal pos
        skip_space()
        start = pos
        while pos < len(text) and text[pos] not in ",();:":
            pos += 1
        return text[start:pos].strip()

    def parse_subtree():
        nonlocal pos
        skip_space()
        if text[pos] == "(":
            pos += 1
            node = Node()
            while True:
                node.children.append(parse_subtree())
                skip_space()
                if text[pos] == ",":
                    pos += 1
                elif text[pos] == ")":
                    pos += 1
                    break
                else:
                    raise ValueError("Invalid Newick structure")
            node.name = read_label()
            skip_branch_length()
            return node
        node = Node(read_label())
        skip_branch_length()
        return node

    root = parse_subtree()
    return root


def prune(node, keep):
    if not node.children:
        return node if node.name in keep else None
    children = [pruned for child in node.children if (pruned := prune(child, keep))]
    if not children:
        return None
    if len(children) == 1:
        return children[0]
    node.name = ""
    node.children = children
    return node


def leaves(node):
    if not node.children:
        return [node]
    result = []
    for child in node.children:
        result.extend(leaves(child))
    return result


def connect(node, graph):
    graph.setdefault(node, [])
    for child in node.children:
        graph[node].append(child)
        graph.setdefault(child, []).append(node)
        connect(child, graph)


def reroot_at_outgroup(tree, outgroup):
    graph = {}
    connect(tree, graph)
    leaf = next((node for node in graph if node.name == outgroup and not node.children), None)
    if leaf is None:
        raise ValueError(f"Outgroup tip not found: {outgroup}")
    parent = graph[leaf][0]

    def orient(node, previous):
        result = Node(node.name if not node.children else "")
        result.children = [orient(neighbor, node) for neighbor in graph[node] if neighbor is not previous]
        return result

    rest = orient(parent, leaf)
    root = Node()
    root.children = [rest, Node(outgroup)]
    return root


def serialize(node):
    if node.children:
        return "(" + ",".join(serialize(child) for child in node.children) + ")"
    return node.name


source_file, sets_file, source_outgroup_tip, output_file = sys.argv[1:]
with open(source_file) as handle:
    tree = parse_tree(handle.read())

populations = set()
with open(sets_file) as handle:
    for line in handle:
        fields = line.split()
        if len(fields) >= 2 and fields[1] != "xxx":
            populations.add(fields[1])

if "Outgroup" not in populations:
    raise ValueError(f"No Outgroup entry in {sets_file}")

keep = (populations - {"Outgroup"}) | {source_outgroup_tip}
tree_leaves = {leaf.name for leaf in leaves(tree)}
missing = sorted(keep - tree_leaves)
if missing:
    raise ValueError("Tree is missing required tips: " + ", ".join(missing))

tree = prune(tree, keep)
for leaf in leaves(tree):
    if leaf.name == source_outgroup_tip:
        leaf.name = "Outgroup"

with open(output_file, "w") as handle:
    handle.write(serialize(tree) + ";\n")
PY
    echo "Created tree: $output_tree"
}

make_tree "Lepus_americanus" "$SETS_DIR/sets_hare.txt" "$AMERICANUS_TREE"
make_tree "Oryctolagus_cuniculus" "$SETS_DIR/sets_krol.txt" "$RABBIT_TREE"

if [[ "${SKIP_DSUITE:-0}" == "1" ]]; then
    echo "SKIP_DSUITE=1: tree preparation only"
    exit 0
fi

run_fbranch() {
    local label=$1
    local tree=$2
    local sets=$3
    local prefix="$OUT_DIR/$label"

    echo "Running Dtrios: $label"
    "$DSUITE" Dtrios -t "$tree" -o "$prefix" "$VCF" "$sets"
    echo "Running Fbranch: $label"
    "$DSUITE" Fbranch "$tree" "${prefix}_tree.txt" | tee "${prefix}_fbranch.txt"
}

run_fbranch "americanus_outgroup" "$AMERICANUS_TREE" "$SETS_DIR/sets_hare.txt"
run_fbranch "rabbit_outgroup" "$RABBIT_TREE" "$SETS_DIR/sets_krol.txt"

echo "Fbranch results saved in: $OUT_DIR"
