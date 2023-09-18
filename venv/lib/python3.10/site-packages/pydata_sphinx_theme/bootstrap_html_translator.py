"""A custom Sphinx HTML Translator for Bootstrap layout
"""
from packaging.version import Version
from docutils import nodes

import sphinx
from sphinx.writers.html5 import HTML5Translator
from sphinx.util import logging
from sphinx.ext.autosummary import autosummary_table

logger = logging.getLogger(__name__)


class BootstrapHTML5Translator(HTML5Translator):
    """Custom HTML Translator for a Bootstrap-ified Sphinx layout
    This is a specialization of the HTML5 Translator of sphinx.
    Only a couple of functions have been overridden to produce valid HTML to be
    directly styled with Bootstrap, and fulfill acessibility best practices.
    """

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.settings.table_style = "table"

    def starttag(self, *args, **kwargs):
        """ensure an aria-level is set for any heading role"""
        if kwargs.get("ROLE") == "heading" and "ARIA-LEVEL" not in kwargs:
            kwargs["ARIA-LEVEL"] = "2"
        return super().starttag(*args, **kwargs)

    def visit_table(self, node):
        # type: (nodes.Element) -> None
        # copy of sphinx source to *not* add 'docutils' and 'align-default' classes
        # but add 'table' class

        # generate_targets_for_table is deprecated in 4.0
        if Version(sphinx.__version__) < Version("4.0"):
            self.generate_targets_for_table(node)

        if Version(sphinx.__version__) < Version("4.3"):
            self._table_row_index = 0
        else:
            self._table_row_indices.append(0)

        classes = [cls.strip(" \t\n") for cls in self.settings.table_style.split(",")]

        # we're looking at the 'real_table', which is wrapped by an autosummary
        if isinstance(node.parent, autosummary_table):
            classes += ["autosummary"]

        # classes.insert(0, "docutils")  # compat
        # if 'align' in node:
        #     classes.append('align-%s' % node['align'])
        tag = self.starttag(node, "table", CLASS=" ".join(classes))
        self.body.append(tag)
