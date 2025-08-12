from docutils import nodes
from docutils.parsers.rst import directives
from docutils.parsers.rst.states import Body
from docutils.statemachine import StringList

from sphinx.application import Sphinx
from sphinx.environment import BuildEnvironment
from sphinx.util.docutils import SphinxDirective
from sphinx.util.nodes import nested_parse_with_titles
from sphinx.util.typing import ExtensionMetadata

from tabulate import tabulate

from scipy._lib._array_api_docs_tables import BackendSupportStatus
from scipy._lib._array_api_docs_tables import calculate_table_statistics
from scipy._lib._array_api_docs_tables import make_flat_capabilities_table
from scipy._lib._public_api import PUBLIC_MODULES


def _make_reST_table(
        rows: list[str],
        headers: list[str],
        state: Body,
) -> list[nodes.Node]:
    """Return list of docutils nodes representing a reST table."""
    rst_table = tabulate(rows, headers=headers, tablefmt="rst")
    string_list = StringList(rst_table.splitlines())
    node = nodes.section()
    node.document = state.document
    nested_parse_with_titles(state, string_list, node)
    return node.children


def _get_flat_table_and_backends(
        env: BuildEnvironment, backend_type: str
) -> tuple[list[dict[str, str]], list[str]]:
    if not hasattr(env, "_array_api_capabilities_table_cache"):
        env._array_api_capabilities_table_cache = {}
    cache = env._array_api_capabilities_table_cache
    if backend_type in cache:
        return cache[backend_type]["flat_table"], cache[backend_type]["backends"]

    included_modules = [
        module for module in PUBLIC_MODULES
        # These are extension modules which should never be included. The introspection
        # in make_flat_capabilities_table will fail for these.
        if module not in {"scipy.linalg.cython_blas", "scipy.linalg.cython_lapack"}
    ]

    flat_table = make_flat_capabilities_table(
        included_modules, backend_type
    )

    if not flat_table:
        backends = []
    else:
        backends = [key for key in flat_table[0]
                    if key not in ["module", "function"]]

    cache[backend_type] = {"flat_table": flat_table, "backends": backends}
    return flat_table, backends


class ArrayAPISupportPerModule(SphinxDirective):
    has_content = False
    option_spec = {
        module_name.replace("scipy.", ""): directives.unchanged
        for module_name in PUBLIC_MODULES
    }
    option_spec["backend_type"] = directives.unchanged

    def run(self):
        backend_type = self.options.pop("backend_type")
        flat_table, backends = _get_flat_table_and_backends(self.env, backend_type)
        table_stats = calculate_table_statistics(flat_table)
        modules = self.options.keys()

        rows = []
        headers = ["module"]
        headers += backends

        for module in modules:
            label = self.options.get(module)
            module_text = f":ref:`{module} <{label}>`" if label else module
            info, accurate_count = table_stats[f"scipy.{module}"]
            total = info.pop("total")
            row = [f"{module_text} ({total})"]
            for backend, count in info.items():
                cell_text = f"{count/total:.0%}{'*' if not accurate_count else ''}"
                row.append(cell_text)
            rows.append(row)
        return _make_reST_table(rows, headers, self.state)


class ArrayAPISupportPerFunction(SphinxDirective):
    has_content = False
    option_spec = {
        "module": directives.unchanged, "backend_type": directives.unchanged
    }

    def _get_generated_doc_link_for_function(self, module, func):
        if module == "signal" and func == "czt":
            pagename = "czt-function"
        else:
            pagename = f"scipy.{module}.{func}"
        return (
            f":doc:`{func} <../../../reference/generated/{pagename}>`"
        )

    def run(self):
        backend_type = self.options["backend_type"]
        flat_table, backends = _get_flat_table_and_backends(self.env, backend_type)
        module = self.options["module"]
        relevant_rows = (
            row for row in flat_table if row["module"] == f"scipy.{module}"
        )

        headers = ["function"]
        headers += backends

        new_rows = []
        S = BackendSupportStatus
        missing_xp_capabilities = False
        for row in relevant_rows:
            func = row["function"]
            new_row = [self._get_generated_doc_link_for_function(module, func)]
            for backend in backends:
                supported = row[backend]
                cell_text = ""
                if supported == S.OUT_OF_SCOPE:
                    cell_text = "N/A"
                elif supported == S.YES:
                    cell_text = "✔️"
                elif supported == S.NO:
                    cell_text = "✖"
                else:
                    missing_xp_capabilities = True
                new_row.append(cell_text)
            new_rows.append(new_row)
        table_nodes = _make_reST_table(new_rows, headers, self.state)

        legend = nodes.admonition()
        legend += nodes.title(text="Legend")
        legend += nodes.paragraph(text="✔️ = supported")
        legend += nodes.paragraph(text="✖ = unsupported")
        legend += nodes.paragraph(text="N/A = out-of-scope")
        if missing_xp_capabilities:
            # Only include blank in the legend if there are actually functions
            # missing the xp_capabilities decorator.
            legend += nodes.paragraph(text="blank = not currently documented")
        return [legend] + table_nodes

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
