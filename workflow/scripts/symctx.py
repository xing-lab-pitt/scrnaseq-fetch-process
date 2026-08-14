#!/usr/bin/env python3
"""Print just the slice of Python source you actually need.

Reading a whole module to find one function spends context on every line you did
not want. This prints the structure, or one named symbol, so the rest never
enters the conversation. Stdlib only (ast) — no tree-sitter, no plugin.

    python symctx.py outline <file.py> [...]      symbols + signatures, bodies folded
    python symctx.py show <file.py> <symbol>      full source of one symbol
    python symctx.py find <dir> <name>            locate a symbol across a tree

`symbol` accepts dotted paths for methods: `MyClass.my_method`.

Designed to be called through context-mode's ctx_execute so the output is indexed
rather than dumped:

    ctx_execute(language="shell",
                code="python symctx.py show pipeline.py build_matrices",
                intent="how build_matrices assembles layers")
"""
from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

DEF = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)


def _sig(node) -> str:
    """One-line signature for a def/class, without its body."""
    if isinstance(node, ast.ClassDef):
        bases = ", ".join(ast.unparse(b) for b in node.bases)
        return f"class {node.name}({bases})" if bases else f"class {node.name}"
    a = node.args
    parts = []
    # positional-only, then normal, aligning defaults to the tail of the list
    pos = a.posonlyargs + a.args
    pad = len(pos) - len(a.defaults)
    for i, arg in enumerate(pos):
        s = arg.arg
        if arg.annotation:
            s += f": {ast.unparse(arg.annotation)}"
        if i >= pad:
            s += f"={ast.unparse(a.defaults[i - pad])}"
        parts.append(s)
        if a.posonlyargs and i == len(a.posonlyargs) - 1:
            parts.append("/")
    if a.vararg:
        parts.append(f"*{a.vararg.arg}")
    elif a.kwonlyargs:
        parts.append("*")
    for arg, d in zip(a.kwonlyargs, a.kw_defaults):
        s = arg.arg
        if arg.annotation:
            s += f": {ast.unparse(arg.annotation)}"
        if d is not None:
            s += f"={ast.unparse(d)}"
        parts.append(s)
    if a.kwarg:
        parts.append(f"**{a.kwarg.arg}")
    kw = "async def" if isinstance(node, ast.AsyncFunctionDef) else "def"
    ret = f" -> {ast.unparse(node.returns)}" if node.returns else ""
    return f"{kw} {node.name}({', '.join(parts)}){ret}"


def _parse(path: Path) -> ast.Module:
    try:
        return ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    except SyntaxError as e:
        sys.exit(f"{path}: syntax error at line {e.lineno}: {e.msg}")


def _walk(body, prefix=""):
    """Yield (dotted_name, node, depth) for defs, descending one level into classes."""
    for node in body:
        if isinstance(node, DEF):
            name = f"{prefix}{node.name}"
            yield name, node, prefix.count(".")
            if isinstance(node, ast.ClassDef):
                yield from _walk(node.body, prefix=f"{name}.")


def cmd_outline(paths) -> int:
    for p in paths:
        path = Path(p)
        if not path.is_file():
            print(f"!! not a file: {path}")
            continue
        tree = _parse(path)
        n_lines = len(path.read_text(encoding="utf-8", errors="replace").splitlines())
        print(f"=== {path}  ({n_lines} lines) ===")
        doc = ast.get_docstring(tree)
        if doc:
            print(f'  """{doc.strip().splitlines()[0]}"""')
        syms = list(_walk(tree.body))
        if not syms:
            print("  (no top-level defs or classes)")
        for name, node, depth in syms:
            end = getattr(node, "end_lineno", node.lineno)
            indent = "  " + "    " * depth
            print(f"{indent}{_sig(node):<70} L{node.lineno}-{end}")
        # module-level assignments are often the config knobs worth seeing
        consts = [t.id for n in tree.body if isinstance(n, ast.Assign)
                  for t in n.targets if isinstance(t, ast.Name) and t.id.isupper()]
        if consts:
            print(f"  module constants: {', '.join(consts)}")
    return 0


def cmd_show(path, symbol) -> int:
    path = Path(path)
    tree = _parse(path)
    src = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for name, node, _ in _walk(tree.body):
        if name == symbol or name.split(".")[-1] == symbol:
            start = min([node.lineno] + [d.lineno for d in node.decorator_list])
            end = getattr(node, "end_lineno", node.lineno)
            print(f"=== {path}:{start}-{end}  {name} ===")
            for i in range(start, end + 1):
                print(f"{i:5d}  {src[i - 1]}")
            return 0
    avail = ", ".join(n for n, _, _ in _walk(tree.body)) or "(none)"
    sys.exit(f"symbol {symbol!r} not found in {path}. Available: {avail}")


def cmd_find(root, name) -> int:
    hits = 0
    for path in sorted(Path(root).rglob("*.py")):
        if any(part in {".git", "__pycache__", ".venv", "node_modules"}
               for part in path.parts):
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
        except SyntaxError:
            continue
        for sym, node, _ in _walk(tree.body):
            if name.lower() in sym.lower():
                print(f"{path}:{node.lineno}  {_sig(node)}")
                hits += 1
    if not hits:
        print(f"no symbol matching {name!r} under {root}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("outline", help="symbols + signatures, bodies folded")
    p.add_argument("paths", nargs="+")

    p = sub.add_parser("show", help="full source of one symbol")
    p.add_argument("path")
    p.add_argument("symbol")

    p = sub.add_parser("find", help="locate a symbol across a tree")
    p.add_argument("root")
    p.add_argument("name")

    a = ap.parse_args()
    if a.cmd == "outline":
        return cmd_outline(a.paths)
    if a.cmd == "show":
        return cmd_show(a.path, a.symbol)
    return cmd_find(a.root, a.name)


if __name__ == "__main__":
    sys.exit(main())
