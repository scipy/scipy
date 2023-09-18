'''doit / click integration through custom class interface

This file contains tasks definitions for doit (https://pydoit.org).
And also a CLI interface using click (https://click.palletsprojects.com).

The CLI is ideal for project contributors while,
doit interface is better suited for authring the development tasks.


# USAGE:

```
class CLI(CliGroup):
    context = CONTEXT
    run_doit_task = run_doit_task

@click.group(cls=CLI)
@click.pass_context
def cli(ctx, **kwargs):
    """Developer Tool"""
    CLI.update_context(ctx, kwargs)
```

## 1 - click API

Commands can added using default Click API. i.e.

```
@cli.command()
@click.argument('extra_argv', nargs=-1)
@click.pass_obj
def python(ctx_obj, extra_argv):
    """Start a Python shell with PYTHONPATH set"""
```

## 2 - class based Click command definition

`CliGroup` provides an alternative class based API to create Click commands.

Just use the `cls_cmd` decorator. And define a `run()` method

```
@cli.cls_cmd('test')
class Test():
    """Run tests"""

    @classmethod
    def run(cls):
        print('Running tests...')
```

- Command may make use a Click.Group context defining a `ctx` class attribute
- Command options are also define as class attributes

```
@cli.cls_cmd('test')
class Test():
    """Run tests"""
    ctx = CONTEXT

    verbose = Option(
        ['--verbose', '-v'], default=False, is_flag=True, help="verbosity")

    @classmethod
    def run(cls, **kwargs): # kwargs contains options from class and CONTEXT
        print('Running tests...')
```

## 3 - class based interface can be run as a doit task by subclassing from Task

- Extra doit task metadata can be defined as class attribute TASK_META.
- `run()` method will be used as python-action by task

```
@cli.cls_cmd('test')
class Test(Task):   # Task base class, doit will create a task
    """Run tests"""
    ctx = CONTEXT

    TASK_META = {
        'task_dep': ['build'],
    }

    @classmethod
    def run(cls, **kwargs):
        pass
```

## 4 - doit tasks with cmd-action "shell" or dynamic metadata

Define method `task_meta()` instead of `run()`:

```
@cli.cls_cmd('refguide-check')
class RefguideCheck(Task):
    @classmethod
    def task_meta(cls, **kwargs):
        return {
```
'''

import re
import sys

import click
from click.globals import get_current_context
from rich_click import RichCommand, RichGroup
from doit.task import Task as DoitTask



class UnifiedContext():
    """Higher level to allow some level of Click.context with doit"""
    def __init__(self, options: dict):
        self.options = options
        self.vals = {}

    def get(self, save=None):
        # try to get from Click
        ctx = get_current_context(silent=True)
        if ctx:
            return ctx.obj
        else:
            if save:
                for name in self.options.keys():
                    if name in save:
                        self.vals[name] = save[name]
            return self.vals


param_type_map = {
    'text': str,
    'boolean': bool,
    'integer': int,
}

def param_click2doit(name: str, val: click.Parameter):
    """Convert click param to dict used by doit.cmdparse"""
    pd = {
        'name': name,
        'type': param_type_map[val.type.name],  # FIXME: add all types
        'default': val.default,
        'help': getattr(val, 'help', ''),
        'metavar': val.metavar,
    }
    for opt in val.opts:
        if opt[:2] == '--':
            pd['long'] = opt[2:]
        elif opt[0] == '-':
            pd['short'] = opt[1:]
    return pd


# convert click.types.ParamType.name to doit param type
CAMEL_PATTERN = re.compile(r'(?!^)([A-Z]+)')
def to_camel(name):
    return CAMEL_PATTERN.sub(r'-\1', name).lower()


def run_as_py_action(cls):
    """used by doit loader to create task instances"""
    if cls is Task:
        return
    task_kwargs = getattr(cls, 'TASK_META', {})
    return DoitTask(
        # convert name to kebab-case
        name=to_camel(cls.__name__),
        doc=cls.__doc__,
        actions=[cls.run],
        params=cls._params,
        **task_kwargs,
    )


class MetaclassDoitTask(type):
    def __new__(meta_cls, name, bases, dct):
        # params/opts from UnifiedContext and Option attributes
        cls = super().__new__(meta_cls, name, bases, dct)
        params = []
        if ctx := getattr(cls, 'ctx', None):
            for ctx_opt in ctx.options.values():
                params.append(param_click2doit(ctx_opt.name, ctx_opt))
        for attr_name, val in cls.__dict__.items():
            if isinstance(val, click.Parameter):
                params.append(param_click2doit(attr_name, val))
        cls._params = params

        task_basename = to_camel(cls.__name__)
        if hasattr(cls, 'task_meta'):
            def creator(**kwargs):
                task_meta = cls.task_meta(**kwargs)
                if 'basename' not in task_meta:
                    task_meta['basename'] = task_basename
                if 'doc' not in task_meta:
                    task_meta['doc'] = cls.__doc__
                return task_meta
            creator._task_creator_params = cls._params
        else:
            def creator():
                return run_as_py_action(cls)
        creator.basename = task_basename
        cls.create_doit_tasks = creator
        return cls


class Task(metaclass=MetaclassDoitTask):
    """Base class to define doit task and/or click command"""

    @classmethod
    def opt_defaults(cls):
        """helper used by another click commands to call this command"""
        return {p['name']: p['default'] for p in cls._params}


class CliGroup(RichGroup):
    """
    Sample for `run_doit_task()`

        from doit.api import run_tasks
        from doit.cmd_base import ModuleTaskLoader

        def run_doit_task(tasks):
            ":param tasks: (dict) task_name -> {options}"
            loader = ModuleTaskLoader(globals())
            doit_config = {
                'verbosity': 2,
                'reporter': 'zero',
            }
            return run_tasks(loader, tasks, extra_config={'GLOBAL': doit_config})
    """
    command_class = RichCommand
    context: UnifiedContext = None


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.params.extend(self.context.options.values())


    @staticmethod
    def run_doit_task():
        raise NotImplementedError()

    @classmethod
    def update_context(cls, ctx: click.core.Context, kwargs: dict):
        """copy received paramaters vals into group context

        Typically to be included in the cli callback
        """
        ctx.ensure_object(dict)
        for opt_name in cls.context.options.keys():
            ctx.obj[opt_name] = kwargs.get(opt_name)


    def cls_cmd(self, name):
        """class decorator, convert to click.Command"""
        def register_click(cls):
            # get options/arguments for class definition
            params = []
            for attr_name, attr_val in cls.__dict__.items():
                if isinstance(attr_val, click.Parameter):
                    params.append(attr_val)

            if issubclass(cls, Task):
                # run as doit task
                def callback(**kwargs):
                    sys.exit(self.run_doit_task.__func__({name: kwargs}))
            else:
                # run as plain function
                def callback(**kwargs):
                    sys.exit(cls.run(**kwargs))

            click_cmd = self.command_class(
                name=name,
                callback=callback,
                help=cls.__doc__,
                params=params,
            )
            self.add_command(click_cmd)
            return cls
        return register_click
