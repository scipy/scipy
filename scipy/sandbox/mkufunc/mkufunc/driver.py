import sys, os

from pypy.translator.translator import TranslationContext, graphof
from pypy.translator.tool.taskengine import SimpleTaskEngine
from pypy.translator.goal import query
from pypy.translator.goal.timing import Timer
from pypy.annotation import model as annmodel
from pypy.annotation.listdef import s_list_of_strings
from pypy.annotation import policy as annpolicy
from py.compat import optparse
from pypy.tool.udir import udir

import py
from pypy.tool.ansi_print import ansi_log
log = py.log.Producer("translation")
py.log.setconsumer("translation", ansi_log)

DEFAULTS = {
  'translation.gc': 'ref',
  'translation.cc': None,
  'translation.profopt': None,

  'translation.thread': False, # influences GC policy

  'translation.stackless': False,
  'translation.debug': True,
  'translation.insist': False,
  'translation.backend': 'c',
  'translation.fork_before': None,
  'translation.backendopt.raisingop2direct_call' : False,
  'translation.backendopt.merge_if_blocks': True,
}


def taskdef(taskfunc, deps, title, new_state=None, expected_states=[],
            idemp=False, earlycheck=None):
    taskfunc.task_deps = deps
    taskfunc.task_title = title
    taskfunc.task_newstate = None
    taskfunc.task_expected_states = expected_states
    taskfunc.task_idempotent = idemp
    taskfunc.task_earlycheck = earlycheck
    return taskfunc

# TODO:
# sanity-checks using states

_BACKEND_TO_TYPESYSTEM = {
    'c': 'lltype',
    'llvm': 'lltype'
}

def backend_to_typesystem(backend):
    return _BACKEND_TO_TYPESYSTEM.get(backend, 'ootype')

# set of translation steps to profile
PROFILE = set([])

class Instrument(Exception):
    pass


class ProfInstrument(object):
    name = "profinstrument"
    def __init__(self, datafile, compiler):
        self.datafile = datafile
        self.compiler = compiler

    def first(self):
        self.compiler._build()

    def probe(self, exe, args):
        from py.compat import subprocess
        env = os.environ.copy()
        env['_INSTRUMENT_COUNTERS'] = str(self.datafile)
        subprocess.call("'%s' %s" % (exe, args), env=env, shell=True)

    def after(self):
        # xxx
        os._exit(0)


class TranslationDriver(SimpleTaskEngine):

    def __init__(self, setopts=None, default_goal=None,
                 disable=[],
                 exe_name=None, extmod_name=None,
                 config=None, overrides=None):
        self.timer = Timer()
        SimpleTaskEngine.__init__(self)

        self.log = log

        if config is None:
            from pypy.config.pypyoption import get_pypy_config
            config = get_pypy_config(DEFAULTS, translating=True)
        self.config = config
        if overrides is not None:
            self.config.override(overrides)

        if setopts is not None:
            self.config.set(**setopts)

        self.exe_name = exe_name
        self.extmod_name = extmod_name

        self.done = {}

        self.disable(disable)

        if default_goal:
            default_goal, = self.backend_select_goals([default_goal])
            if default_goal in self._maybe_skip():
                default_goal = None

        self.default_goal = default_goal
        self.extra_goals = []
        self.exposed = []

        # expose tasks
        def expose_task(task, backend_goal=None):
            if backend_goal is None:
                backend_goal = task
            def proc():
                return self.proceed(backend_goal)
            self.exposed.append(task)
            setattr(self, task, proc)

        backend, ts = self.get_backend_and_type_system()
        for task in self.tasks:
            explicit_task = task
            parts = task.split('_')
            if len(parts) == 1:
                if task in ('annotate'):
                    expose_task(task)
            else:
                task, postfix = parts
                if task in ('rtype', 'backendopt', 'llinterpret',
                            'prehannotatebackendopt', 'hintannotate',
                            'timeshift'):
                    if ts:
                        if ts == postfix:
                            expose_task(task, explicit_task)
                    else:
                        expose_task(explicit_task)
                elif task in ('source', 'compile', 'run'):
                    if backend:
                        if backend == postfix:
                            expose_task(task, explicit_task)
                    elif ts:
                        if ts == backend_to_typesystem(postfix):
                            expose_task(explicit_task)
                    else:
                        expose_task(explicit_task)

    def set_extra_goals(self, goals):
        self.extra_goals = goals

    def get_info(self): # XXX more?
        d = {'backend': self.config.translation.backend}
        return d

    def get_backend_and_type_system(self):
        type_system = self.config.translation.type_system
        backend = self.config.translation.backend
        return backend, type_system

    def backend_select_goals(self, goals):
        backend, ts = self.get_backend_and_type_system()
        postfixes = [''] + ['_'+p for p in (backend, ts) if p]
        l = []
        for goal in goals:
            for postfix in postfixes:
                cand = "%s%s" % (goal, postfix)
                if cand in self.tasks:
                    new_goal = cand
                    break
            else:
                raise Exception, "cannot infer complete goal from: %r" % goal
            l.append(new_goal)
        return l

    def disable(self, to_disable):
        self._disabled = to_disable

    def _maybe_skip(self):
        maybe_skip = []
        if self._disabled:
            for goal in  self.backend_select_goals(self._disabled):
                maybe_skip.extend(self._depending_on_closure(goal))
        return dict.fromkeys(maybe_skip).keys()


    def setup(self, entry_point, inputtypes, policy=None, extra={}, empty_translator=None):
        standalone = inputtypes is None
        self.standalone = standalone

        if standalone:
            inputtypes = [s_list_of_strings]
        self.inputtypes = inputtypes

        if policy is None:
            policy = annpolicy.AnnotatorPolicy()
        if standalone:
            policy.allow_someobjects = False
        self.policy = policy

        self.extra = extra

        if empty_translator:
            translator = empty_translator
        else:
            translator = TranslationContext(config=self.config)

        self.entry_point = entry_point
        self.translator = translator
        self.libdef = None

        self.translator.driver_instrument_result = self.instrument_result

    def setup_library(self, libdef, policy=None, extra={}, empty_translator=None):
        self.setup(None, None, policy, extra, empty_translator)
        self.libdef = libdef

    def instrument_result(self, args):
        backend, ts = self.get_backend_and_type_system()
        if backend != 'c' or sys.platform == 'win32':
            raise Exception("instrumentation requires the c backend"
                            " and unix for now")
        from pypy.tool.udir import udir

        datafile = udir.join('_instrument_counters')
        makeProfInstrument = lambda compiler: ProfInstrument(datafile, compiler)

        pid = os.fork()
        if pid == 0:
            # child compiling and running with instrumentation
            self.config.translation.instrument = True
            self.config.translation.instrumentctl = (makeProfInstrument,
                                                     args)
            raise Instrument
        else:
            pid, status = os.waitpid(pid, 0)
            if os.WIFEXITED(status):
                status = os.WEXITSTATUS(status)
                if status != 0:
                    raise Exception, "instrumentation child failed: %d" % status
            else:
                raise Exception, "instrumentation child aborted"
            import array, struct
            n = datafile.size()//struct.calcsize('L')
            datafile = datafile.open('rb')
            counters = array.array('L')
            counters.fromfile(datafile, n)
            datafile.close()
            return counters

    def info(self, msg):
        log.info(msg)

    def _profile(self, goal, func):
        from cProfile import Profile
        from pypy.tool.lsprofcalltree import KCacheGrind
        d = {'func':func}
        prof = Profile()
        prof.runctx("res = func()", globals(), d)
        KCacheGrind(prof).output(open(goal + ".out", "w"))
        return d['res']

    def _do(self, goal, func, *args, **kwds):
        title = func.task_title
        if goal in self.done:
            self.log.info("already done: %s" % title)
            return
        else:
            self.log.info("%s..." % title)
        self.timer.start_event(goal)
        try:
            instrument = False
            try:
                if goal in PROFILE:
                    res = self._profile(goal, func)
                else:
                    res = func()
            except Instrument:
                instrument = True
            if not func.task_idempotent:
                self.done[goal] = True
            if instrument:
                self.proceed('compile')
                assert False, 'we should not get here'
        finally:
            self.timer.end_event(goal)
        return res

    def task_annotate(self):
        # includes annotation and annotatation simplifications
        translator = self.translator
        policy = self.policy
        self.log.info('with policy: %s.%s' %
                      (policy.__class__.__module__, policy.__class__.__name__))

        annmodel.DEBUG = self.config.translation.debug
        annotator = translator.buildannotator(policy=policy)

        if self.entry_point:
            s = annotator.build_types(self.entry_point, self.inputtypes)

            self.sanity_check_annotation()
            if self.standalone and s.knowntype != int:
                raise Exception("stand-alone program entry point must return an "
                                "int (and not, e.g., None or always raise an "
                                "exception).")
            annotator.simplify()
            return s
        else:
            assert self.libdef is not None
            for func, inputtypes in self.libdef.functions:
                annotator.build_types(func, inputtypes)
            self.sanity_check_annotation()
            annotator.simplify()
    #
    task_annotate = taskdef(task_annotate, [], "Annotating&simplifying")


    def sanity_check_annotation(self):
        translator = self.translator
        irreg = query.qoutput(query.check_exceptblocks_qgen(translator))
        if irreg:
            self.log.info("Some exceptblocks seem insane")

        lost = query.qoutput(query.check_methods_qgen(translator))
        assert not lost, "lost methods, something gone wrong with the annotation of method defs"

        so = query.qoutput(query.polluted_qgen(translator))
        tot = len(translator.graphs)
        percent = int(tot and (100.0*so / tot) or 0)
        # if there are a few SomeObjects even if the policy doesn't allow
        # them, it means that they were put there in a controlled way
        # and then it's not a warning.
        if not translator.annotator.policy.allow_someobjects:
            pr = self.log.info
        elif percent == 0:
            pr = self.log.info
        else:
            pr = log.WARNING
        pr("-- someobjectness %2d%% (%d of %d functions polluted by SomeObjects)" % (percent, so, tot))



    def task_rtype_lltype(self):
        rtyper = self.translator.buildrtyper(type_system='lltype')
        insist = not self.config.translation.insist
        rtyper.specialize(dont_simplify_again=True,
                          crash_on_first_typeerror=insist)
    #
    task_rtype_lltype = taskdef(task_rtype_lltype, ['annotate'], "RTyping")
    RTYPE = 'rtype_lltype'

    def task_rtype_ootype(self):
        # Maybe type_system should simply be an option used in task_rtype
        insist = not self.config.translation.insist
        rtyper = self.translator.buildrtyper(type_system="ootype")
        rtyper.specialize(dont_simplify_again=True,
                          crash_on_first_typeerror=insist)
    #
    task_rtype_ootype = taskdef(task_rtype_ootype, ['annotate'], "ootyping")
    OOTYPE = 'rtype_ootype'

    def task_prehannotatebackendopt_lltype(self):
        from pypy.translator.backendopt.all import backend_optimizations
        backend_optimizations(self.translator,
                              inline_threshold=0,
                              merge_if_blocks=True,
                              constfold=True,
                              raisingop2direct_call=False,
                              remove_asserts=True)
    #
    task_prehannotatebackendopt_lltype = taskdef(
        task_prehannotatebackendopt_lltype,
        [RTYPE],
        "Backendopt before Hint-annotate")

    def task_hintannotate_lltype(self):
        from pypy.jit.hintannotator.annotator import HintAnnotator
        from pypy.jit.hintannotator.model import OriginFlags
        from pypy.jit.hintannotator.model import SomeLLAbstractConstant

        get_portal = self.extra['portal']
        PORTAL, POLICY = get_portal(self)
        t = self.translator
        self.portal_graph = graphof(t, PORTAL)

        hannotator = HintAnnotator(base_translator=t, policy=POLICY)
        self.hint_translator = hannotator.translator
        hs = hannotator.build_types(self.portal_graph,
                                    [SomeLLAbstractConstant(v.concretetype,
                                                            {OriginFlags(): True})
                                     for v in self.portal_graph.getargs()])
        count = hannotator.bookkeeper.nonstuboriggraphcount
        stubcount = hannotator.bookkeeper.stuboriggraphcount
        self.log.info("The hint-annotator saw %d graphs"
                      " (and made stubs for %d graphs)." % (count, stubcount))
        n = len(list(hannotator.translator.graphs[0].iterblocks()))
        self.log.info("portal has %d blocks" % n)
        self.hannotator = hannotator
    #
    task_hintannotate_lltype = taskdef(task_hintannotate_lltype,
                                       ['prehannotatebackendopt_lltype'],
                                       "Hint-annotate")

    def task_timeshift_lltype(self):
        from pypy.jit.timeshifter.hrtyper import HintRTyper
        from pypy.jit.codegen import detect_cpu
        cpu = detect_cpu.autodetect()
        if cpu == 'i386':
            from pypy.jit.codegen.i386.rgenop import RI386GenOp as RGenOp
            RGenOp.MC_SIZE = 32 * 1024 * 1024
        elif cpu == 'ppc':
            from pypy.jit.codegen.ppc.rgenop import RPPCGenOp as RGenOp
            RGenOp.MC_SIZE = 32 * 1024 * 1024
        else:
            raise Exception('Unsuported cpu %r'%cpu)

        del self.hint_translator
        ha = self.hannotator
        t = self.translator
        # make the timeshifted graphs
        hrtyper = HintRTyper(ha, t.rtyper, RGenOp)
        hrtyper.specialize(origportalgraph=self.portal_graph, view=False)
    #
    task_timeshift_lltype = taskdef(task_timeshift_lltype,
                             ["hintannotate_lltype"],
                             "Timeshift")

    def task_backendopt_lltype(self):
        from pypy.translator.backendopt.all import backend_optimizations
        backend_optimizations(self.translator)
    #
    task_backendopt_lltype = taskdef(task_backendopt_lltype,
                                     [RTYPE,
                                      '??timeshift_lltype'],
                                     "lltype back-end optimisations")
    BACKENDOPT = 'backendopt_lltype'

    def task_backendopt_ootype(self):
        from pypy.translator.backendopt.all import backend_optimizations
        backend_optimizations(self.translator)
    #
    task_backendopt_ootype = taskdef(task_backendopt_ootype,
                                        [OOTYPE], "ootype back-end optimisations")
    OOBACKENDOPT = 'backendopt_ootype'


    def task_stackcheckinsertion_lltype(self):
        from pypy.translator.transform import insert_ll_stackcheck
        count = insert_ll_stackcheck(self.translator)
        self.log.info("inserted %d stack checks." % (count,))

    task_stackcheckinsertion_lltype = taskdef(
        task_stackcheckinsertion_lltype,
        ['?'+BACKENDOPT, RTYPE, 'annotate'],
        "inserting stack checks")
    STACKCHECKINSERTION = 'stackcheckinsertion_lltype'

    def possibly_check_for_boehm(self):
        if self.config.translation.gc == "boehm":
            from pypy.translator.tool.cbuild import check_boehm_presence
            from pypy.translator.tool.cbuild import CompilationError
            try:
                check_boehm_presence(noerr=False)
            except CompilationError, e:
                i = 'Boehm GC not installed.  Try e.g. "translate.py --gc=hybrid"'
                raise CompilationError('%s\n--------------------\n%s' % (e, i))

    def task_database_c(self):
        translator = self.translator
        if translator.annotator is not None:
            translator.frozen = True

        standalone = self.standalone

        if standalone:
            from pypy.translator.c.genc import CStandaloneBuilder as CBuilder
        else:
            from pypy.translator.c.genc import CExtModuleBuilder as CBuilder
        cbuilder = CBuilder(self.translator, self.entry_point,
                            config=self.config)
        cbuilder.stackless = self.config.translation.stackless
        if not standalone:     # xxx more messy
            cbuilder.modulename = self.extmod_name
        database = cbuilder.build_database()
        self.log.info("database for generating C source was created")
        self.cbuilder = cbuilder
        self.database = database
    #
    task_database_c = taskdef(task_database_c,
                            [STACKCHECKINSERTION, '?'+BACKENDOPT, RTYPE, '?annotate'],
                            "Creating database for generating c source",
                            earlycheck = possibly_check_for_boehm)

    def task_source_c(self):  # xxx messy
        translator = self.translator
        cbuilder = self.cbuilder
        database = self.database
        c_source_filename = cbuilder.generate_source(database)
        self.log.info("written: %s" % (c_source_filename,))
        self.c_source_filename = str(c_source_filename)
    #
    task_source_c = taskdef(task_source_c, ['database_c'], "Generating c source")

    def task_compile_c(self): # xxx messy
        cbuilder = self.cbuilder
        cbuilder.compile()

        if self.standalone:
            self.c_entryp = cbuilder.executable_name
            self.create_exe()
        else:
            self.c_entryp = cbuilder.get_entry_point()
    #
    task_compile_c = taskdef(task_compile_c, ['source_c'], "Compiling c source")


    def task_run_c(self):
        self.backend_run('c')
    #
    task_run_c = taskdef(task_run_c, ['compile_c'],
                         "Running compiled c source",
                         idemp=True)

    def task_llinterpret_lltype(self):
        from pypy.rpython.llinterp import LLInterpreter
        py.log.setconsumer("llinterp operation", None)

        translator = self.translator
        interp = LLInterpreter(translator.rtyper)
        bk = translator.annotator.bookkeeper
        graph = bk.getdesc(self.entry_point).getuniquegraph()
        v = interp.eval_graph(graph,
                              self.extra.get('get_llinterp_args',
                                             lambda: [])())

        log.llinterpret.event("result -> %s" % v)
    #
    task_llinterpret_lltype = taskdef(task_llinterpret_lltype,
                                      [STACKCHECKINSERTION, '?'+BACKENDOPT, RTYPE],
                                      "LLInterpreting")

    def task_source_llvm(self):
        translator = self.translator
        if translator.annotator is None:
            raise ValueError, "llvm requires annotation."

        from pypy.translator.llvm import genllvm

        self.llvmgen = genllvm.GenLLVM(translator, self.standalone)

        llvm_filename = self.llvmgen.gen_source(self.entry_point)
        self.log.info("written: %s" % (llvm_filename,))
    #
    task_source_llvm = taskdef(task_source_llvm,
                               [STACKCHECKINSERTION, BACKENDOPT, RTYPE],
                               "Generating llvm source")

    def task_compile_llvm(self):
        gen = self.llvmgen
        if self.standalone:
            exe_name = (self.exe_name or 'testing') % self.get_info()
            self.c_entryp = gen.compile_standalone(exe_name)
            self.create_exe()
        else:
            self.c_module, self.c_entryp = gen.compile_module()
    #
    task_compile_llvm = taskdef(task_compile_llvm,
                                ['source_llvm'],
                                "Compiling llvm source")

    def task_run_llvm(self):
        self.backend_run('llvm')
    #
    task_run_llvm = taskdef(task_run_llvm, ['compile_llvm'],
                            "Running compiled llvm source",
                            idemp=True)

    def task_source_js(self):
        from pypy.translator.js.js import JS
        self.gen = JS(self.translator, functions=[self.entry_point],
                      stackless=self.config.translation.stackless)
        filename = self.gen.write_source()
        self.log.info("Wrote %s" % (filename,))
    task_source_js = taskdef(task_source_js,
                        [OOTYPE],
                        'Generating Javascript source')

    def task_compile_js(self):
        pass
    task_compile_js = taskdef(task_compile_js, ['source_js'],
                              'Skipping Javascript compilation')

    def task_run_js(self):
        pass
    task_run_js = taskdef(task_run_js, ['compile_js'],
                              'Please manually run the generated code')

    def task_source_cli(self):
        from pypy.translator.cli.gencli import GenCli
        from pypy.translator.cli.entrypoint import get_entrypoint

        if self.entry_point is not None: # executable mode
            entry_point_graph = self.translator.graphs[0]
            entry_point = get_entrypoint(entry_point_graph)
        else:
            # library mode
            assert self.libdef is not None
            bk = self.translator.annotator.bookkeeper
            entry_point = self.libdef.get_entrypoint(bk)

        self.gen = GenCli(udir, self.translator, entry_point, config=self.config)
        filename = self.gen.generate_source()
        self.log.info("Wrote %s" % (filename,))
    task_source_cli = taskdef(task_source_cli, ["?" + OOBACKENDOPT, OOTYPE],
                             'Generating CLI source')

    def task_compile_cli(self):
        from pypy.translator.oosupport.support import unpatch_os
        from pypy.translator.cli.test.runtest import CliFunctionWrapper
        filename = self.gen.build_exe()
        self.c_entryp = CliFunctionWrapper(filename)
        # restore original os values
        if hasattr(self, 'old_cli_defs'):
            unpatch_os(self.old_cli_defs)

        self.log.info("Compiled %s" % filename)
        if self.standalone and self.exe_name:
            self.copy_cli_exe()
    task_compile_cli = taskdef(task_compile_cli, ['source_cli'],
                              'Compiling CLI source')

    def task_run_cli(self):
        pass
    task_run_cli = taskdef(task_run_cli, ['compile_cli'],
                              'XXX')

    def task_source_jvm(self):
        from pypy.translator.jvm.genjvm import GenJvm
        from pypy.translator.jvm.node import EntryPoint

        entry_point_graph = self.translator.graphs[0]
        is_func = not self.standalone
        entry_point = EntryPoint(entry_point_graph, is_func, is_func)
        self.gen = GenJvm(udir, self.translator, entry_point)
        self.jvmsource = self.gen.generate_source()
        self.log.info("Wrote JVM code")
    task_source_jvm = taskdef(task_source_jvm, ["?" + OOBACKENDOPT, OOTYPE],
                             'Generating JVM source')

    def task_compile_jvm(self):
        from pypy.translator.oosupport.support import unpatch_os
        from pypy.translator.jvm.test.runtest import JvmGeneratedSourceWrapper
        self.jvmsource.compile()
        self.c_entryp = JvmGeneratedSourceWrapper(self.jvmsource)
        # restore original os values
        if hasattr(self, 'old_cli_defs'):
            unpatch_os(self.old_cli_defs)
        self.log.info("Compiled JVM source")
        if self.standalone and self.exe_name:
            self.copy_jvm_jar()
    task_compile_jvm = taskdef(task_compile_jvm, ['source_jvm'],
                              'Compiling JVM source')

    def task_run_jvm(self):
        pass
    task_run_jvm = taskdef(task_run_jvm, ['compile_jvm'],
                           'XXX')

    def proceed(self, goals):
        if not goals:
            if self.default_goal:
                goals = [self.default_goal]
            else:
                self.log.info("nothing to do")
                return
        elif isinstance(goals, str):
            goals = [goals]
        goals.extend(self.extra_goals)
        goals = self.backend_select_goals(goals)
        return self._execute(goals, task_skip = self._maybe_skip())

    def from_targetspec(targetspec_dic, config=None, args=None,
                        empty_translator=None,
                        disable=[],
                        default_goal=None):
        if args is None:
            args = []

        driver = TranslationDriver(config=config, default_goal=default_goal,
                                   disable=disable)
        # patch some attributes of the os module to make sure they
        # have the same value on every platform.
        backend, ts = driver.get_backend_and_type_system()
        if backend in ('cli', 'jvm'):
            from pypy.translator.oosupport.support import patch_os
            driver.old_cli_defs = patch_os()

        target = targetspec_dic['target']
        spec = target(driver, args)

        try:
            entry_point, inputtypes, policy = spec
        except ValueError:
            entry_point, inputtypes = spec
            policy = None

        driver.setup(entry_point, inputtypes,
                     policy=policy,
                     extra=targetspec_dic,
                     empty_translator=empty_translator)

        return driver

    from_targetspec = staticmethod(from_targetspec)

    def prereq_checkpt_rtype(self):
        assert 'pypy.rpython.rmodel' not in sys.modules, (
            "cannot fork because the rtyper has already been imported")
    prereq_checkpt_rtype_lltype = prereq_checkpt_rtype
    prereq_checkpt_rtype_ootype = prereq_checkpt_rtype
