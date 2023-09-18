# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import pkgutil

from . import commands
from . import plugins


class PluginManager(object):
    """
    A class to load and manage plugins.

    By default in asv, plugins are searched for in the `asv.plugins`
    namespace package and in the `asv.commands` package.

    Then, any modules specified in the ``plugins`` entry in the
    ``asv.conf.json`` file are loaded.
    """
    def __init__(self):
        self._plugins = []

    def load_plugins(self, package):
        prefix = package.__name__ + '.'
        for module_finder, name, ispkg in pkgutil.iter_modules(package.__path__, prefix):
            __import__(name)
            mod = sys.modules[name]
            self.init_plugin(mod)
            self._plugins.append(mod)

    def import_plugin(self, name):
        extended = False
        if name.startswith('.'):
            extended = True
            sys.path.insert(0, '.')
            name = name[1:]
        try:
            mod = __import__(name, {}, {}, [], level=0)
            self.init_plugin(mod)
            self._plugins.append(mod)
        finally:
            if extended:
                del sys.path[0]

    def init_plugin(self, mod):
        if hasattr(mod, 'setup'):
            mod.setup()

    def run_hook(self, hook_name, args, kwargs):
        for plugin in self._plugins:
            if hasattr(plugin, hook_name):
                getattr(plugin, hook_name)(*args, **kwargs)


plugin_manager = PluginManager()
plugin_manager.load_plugins(commands)
plugin_manager.load_plugins(plugins)

commands.__doc__ = commands._make_docstring()
