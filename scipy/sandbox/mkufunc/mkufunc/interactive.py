import driver

from pypy.translator.translator import TranslationContext


DEFAULTS = {
  'translation.backend': None,
  'translation.type_system': None,
  'translation.verbose': False,
}

class Translation(object):

    def __init__(self, entry_point, argtypes=None, **kwds):
        self.driver = driver.TranslationDriver(overrides=DEFAULTS)
        self.config = self.driver.config

        self.entry_point = entry_point
        self.context = TranslationContext(config=self.config)

        # hook into driver events
        driver_own_event = self.driver._event
        def _event(kind, goal, func):
            self.driver_event(kind, goal, func)
            driver_own_event(kind, goal, func)
        self.driver._event = _event
        self.driver_setup = False

        self.update_options(argtypes, kwds)
        # for t.view() to work just after construction
        graph = self.context.buildflowgraph(entry_point)
        self.context._prebuilt_graphs[entry_point] = graph

    def driver_event(self, kind, goal, func):
        if kind == 'pre':
            self.ensure_setup()

    def ensure_setup(self, argtypes=None, policy=None, standalone=False):
        if not self.driver_setup:
            if standalone:
                assert argtypes is None
            else:
                if argtypes is None:
                    argtypes = []
            self.driver.setup(self.entry_point, argtypes, policy,
                              empty_translator=self.context)
            self.ann_argtypes = argtypes
            self.ann_policy = policy
            self.driver_setup = True
        else:
            # check consistency
            if standalone:
                assert argtypes is None
                assert self.ann_argtypes is None
            elif argtypes is not None and argtypes != self.ann_argtypes:
                raise Exception("inconsistent argtype supplied")
            if policy is not None and policy != self.ann_policy:
                raise Exception("inconsistent annotation polish supplied")

    def update_options(self, argtypes, kwds):
        if argtypes or kwds.get('policy') or kwds.get('standalone'):
            self.ensure_setup(argtypes, kwds.get('policy'),
                                        kwds.get('standalone'))
        kwds.pop('policy', None)
        kwds.pop('standalone', None)
        self.config.translation.set(**kwds)

    def ensure_opt(self, name, value=None, fallback=None):
        if value is not None:
            self.update_options(None, {name: value})
            return value
        val = getattr(self.config.translation, name, None)
        if fallback is not None and val is None:
            self.update_options(None, {name: fallback})
            return fallback
        if val is not None:
            return val
        raise Exception(
              "the %r option should have been specified at this point" %name)

    def ensure_type_system(self, type_system=None):
        if self.config.translation.backend is not None:
            return self.ensure_opt('type_system')
        return self.ensure_opt('type_system', type_system, 'lltype')

    def ensure_backend(self, backend=None):
        backend = self.ensure_opt('backend', backend)
        self.ensure_type_system()
        return backend

    # backend independent

    def annotate(self, argtypes=None, **kwds):
        self.update_options(argtypes, kwds)
        return self.driver.annotate()

    # type system dependent

    def rtype(self, argtypes=None, **kwds):
        self.update_options(argtypes, kwds)
        ts = self.ensure_type_system()
        return getattr(self.driver, 'rtype_'+ts)()

    # backend depedent

    def source(self, argtypes=None, **kwds):
        self.update_options(argtypes, kwds)
        backend = self.ensure_backend()
        self.driver.source_c()

    def compile(self, argtypes=None, **kwds):
        self.update_options(argtypes, kwds)
        backend = self.ensure_backend()
        self.driver.compile_c()
        return self.driver.c_entryp
