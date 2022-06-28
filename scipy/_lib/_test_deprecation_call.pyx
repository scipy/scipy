from ._test_deprecation_def cimport foo, foo_deprecated

def call():
    return foo(), foo_deprecated()
