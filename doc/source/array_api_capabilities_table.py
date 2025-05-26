from docutils import nodes
from docutils.parsers.rst import directives
from docutils.statemachine import StringList

from sphinx.application import Sphinx
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import nested_parse_with_titles
from sphinx.util.typing import ExtensionMetadata

from tabulate import tabulate

from scipy._lib._array_api_docs_tables import calculate_table_statistics
from scipy._lib._array_api_docs_tables import make_flat_capabilities_table
from scipy._lib.tests.test_public_api import PUBLIC_MODULES


_flat_capabilities_table = None
_backends = None

def _get_flat_table_and_backends() -> tuple[list[dict[str, str]], list[str]]:
    global _flat_capabilities_table
    global _backends
    if _flat_capabilities_table is not None:
        return _flat_capabilities_table, _backends
    included_modules = [
        module for module in PUBLIC_MODULES
        if module not in {"scipy.linalg.cython_blas", "scipy.linalg.cython_lapack"}
    ]
    _flat_capabilities_table = make_flat_capabilities_table(included_modules)
    if not _flat_capabilities_table:
        return _flat_capabilities_table, []
    
    _backends = [key for key in _flat_capabilities_table[0]
                 if key not in ["module", "function"]]
    
    return _flat_capabilities_table, _backends


# Shortened names for use in table.
backend_names_map = {
    "jax.numpy cpu": "jax cpu",
    "jax.numpy gpu": "jax gpu",
    "jax.numpy jit": "jax jit",
    "dask.array": "dask",
    "dask.array lazy": "dask lazy",
}


class ArrayAPISupportPerModule(SphinxDirective):
    has_content = False
    option_spec = {
        module_name.replace("scipy.", ""): directives.unchanged
        for module_name in PUBLIC_MODULES
    }

    def run(self):
        flat_table, backends = _get_flat_table_and_backends()
        table_stats = calculate_table_statistics(flat_table)
        modules = self.options.keys()

        rows = []
        headers = ["module"]
        headers += [
            backend_names_map.get(backend, backend) for backend in backends
        ]
        for module in modules:
            link = self.options.get(module)
            module_text = f":doc:`{module} <{link}>`" if link else module
            row = [module_text]
            info = table_stats[f"scipy.{module}"]
            total = info.pop("total")
            for backend, count in info.items():
                cell_text = f"{count}/{total}"
                row.append(cell_text)
            rows.append(row)
        rst_table = tabulate(rows, headers=headers, tablefmt="rst")
        string_list = StringList(rst_table.splitlines())
        node = nodes.section()
        node.document = self.state.document
        nested_parse_with_titles(self.state, string_list, node)
        return node.children


class ArrayAPISupportPerFunction(SphinxDirective):
    has_content = False
    option_spec = {"module": directives.unchanged}

    def run(self):
        flat_table, backends = _get_flat_table_and_backends()
        module = self.options["module"]
        relevant_rows = (
            row for row in flat_table if row["module"] == f"scipy.{module}"
        )
        print(len(relevant_rows))

        headers = ["function"]
        headers += [
            backend_names_map.get(backend, backend) for backend in backends
        ]
        new_rows = []
        for row in relevant_rows:
            func = row["function"]
            new_row = [func]
            for backend in backends:
                supported = row["backend"]
                cell_text = "✅" if supported else "⛔"
                row.append(cell_text)
            rows.append(row)
        rst_table = tabulate(rows, headers=headers, tablefmt="rst")
        string_list = StringList(rst_table.splitlines())
        node = nodes.section()
        node.document = self.state.document
        nested_parse_with_titles(self.state, string_list, node)
        return node.children


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_directive(
        "array-api-support-per-module", ArrayAPISupportPerModule
    )
    app.add_directive(
        "array-api-support-per-function", ArrayAPISupportPerFunction
    )
    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
