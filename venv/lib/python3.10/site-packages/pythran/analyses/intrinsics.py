""" Intrinsics gathers all intrinsics referenced in a module. """

from pythran.passmanager import ModuleAnalysis
import pythran.intrinsic as intrinsic
from pythran.utils import attr_to_path


class Intrinsics(ModuleAnalysis):

    """ Gather all intrinsics used in the module
    """

    def __init__(self):
        """ Result is a set of intrinsic values. """
        self.result = set()
        super(Intrinsics, self).__init__()

    def visit_Attribute(self, node):
        obj, _ = attr_to_path(node)
        if isinstance(obj, intrinsic.Intrinsic):
            self.result.add(obj)
        self.generic_visit(node)
