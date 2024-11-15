#!/usr/bin/env python

"""
MIT License

Copyright (c) 2020 Marco Gorelli

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Check that test names start with `test`, and that test classes start with
`Test`.
"""
from __future__ import annotations

import ast
import os
from pathlib import Path
import sys
from collections.abc import Iterator, Sequence
import itertools

PRAGMA = "# skip name check"


def _find_names(node: ast.Module) -> Iterator[str]:
    for _node in ast.walk(node):
        if isinstance(_node, ast.Name):
            yield _node.id
        elif isinstance(_node, ast.Attribute):
            yield _node.attr


def _is_fixture(node: ast.expr) -> bool:
    if isinstance(node, ast.Call):
        node = node.func
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "fixture"
        and isinstance(node.value, ast.Name)
        and node.value.id == "pytest"
    )


def is_misnamed_test_func(
    node: ast.expr | ast.stmt, names: Sequence[str], line: str
) -> bool:
    return (
        isinstance(node, ast.FunctionDef)
        and not node.name.startswith("test")
        and names.count(node.name) == 0
        and not any(
            _is_fixture(decorator) for decorator in node.decorator_list
            )
        and PRAGMA not in line
        and node.name
        not in ("teardown_method", "setup_method",
                "teardown_class", "setup_class",
                "setup_module", "teardown_module",
                "_test_dependency",  # array_api_compat.tests.test_no_dependencies
            )
    )


def is_misnamed_test_class(
    node: ast.expr | ast.stmt, names: Sequence[str], line: str
) -> bool:
    return (
        isinstance(node, ast.ClassDef)
        and not node.name.startswith("Test")
        and names.count(node.name) == 0
        and PRAGMA not in line
        # Some of the KDTreeTest use a decorator to setup tests so these are
        # actually fine
        and "KDTreeTest" not in [
            decorator.id for decorator in node.decorator_list
            ]
    )


def main(content: str, file: str) -> int:
    lines = content.splitlines()
    tree = ast.parse(content)
    names = list(_find_names(tree))
    ret = 0
    for node in tree.body:
        if is_misnamed_test_func(node, names, lines[node.lineno - 1]):
            print(
                f"{file}:{node.lineno}:{node.col_offset} "
                f"found test function '{node.name}' which does not start with"
                " 'test'"
            )
            ret = 1
        elif is_misnamed_test_class(node, names, lines[node.lineno - 1]):
            print(
                f"{file}:{node.lineno}:{node.col_offset} "
                f"found test class '{node.name}' which does not start with"
                " 'Test'"
            )
            ret = 1
        if (
            isinstance(node, ast.ClassDef)
            and names.count(node.name) == 0
            and PRAGMA not in lines[node.lineno - 1]
        ):
            for _node in node.body:
                if is_misnamed_test_func(_node, names,
                                         lines[_node.lineno - 1]):
                    # It could be that this function is used somewhere by the
                    # parent class. For example, there might be a base class
                    # with
                    #
                    # class Foo:
                    #     def foo(self):
                    #          assert 1+1==2
                    #     def test_foo(self):
                    #         self.foo()
                    #
                    # and then some subclass overwrites `foo`. So, we check
                    # that `self.foo` doesn't appear in any of the test
                    # classes. Note some false negatives might get through,
                    # but that's OK. This is good enough that has helped
                    # identify several examples of tests not being run.
                    should_continue = False
                    for _file in itertools.chain(
                        Path("scipy").rglob("**/tests/**/test*.py"),
                        ["scipy/_lib/_testutils.py"],
                    ):
                        with open(os.path.join(_file)) as fd:
                            _content = fd.read()
                        if f"self.{_node.name}" in _content:
                            should_continue = True
                            break
                    if should_continue:
                        continue

                    print(
                        f"{file}:{_node.lineno}:{_node.col_offset} "
                        f"found test function '{_node.name}' which does not "
                        "start with 'test'"
                    )
                    ret = 1
    return ret


if __name__ == "__main__":
    ret = 0
    path = Path("scipy").rglob("**/tests/**/test*.py")

    for file in path:
        filename = os.path.basename(file)
        with open(file, encoding="utf-8") as fd:
            content = fd.read()
        ret |= main(content, file)

    sys.exit(ret)
