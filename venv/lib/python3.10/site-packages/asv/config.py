# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys

from .console import log
from . import util

# TODO: Some verification of the config values

class Config(object):
    """
    Manages the configuration for a benchmark project.
    """
    api_version = 1

    def __init__(self):
        self.project = "project"
        self.project_url = "#"
        self.repo = None
        self.repo_subdir = ""
        self.branches = [None]
        self.pythons = ["{0[0]}.{0[1]}".format(sys.version_info)]
        self.matrix = {}
        self.exclude = []
        self.include = []
        self.env_dir = "env"
        self.benchmark_dir = "benchmarks"
        self.results_dir = "results"
        self.html_dir = "html"
        self.show_commit_url = "#"
        self.hash_length = 8
        self.environment_type = None
        self.install_timeout = 600
        self.dvcs = None
        self.regressions_first_commits = {}
        self.regressions_thresholds = {}
        self.plugins = []
        self.conda_channels = []
        self.build_command = None
        self.install_command = None
        self.uninstall_command = None

    @classmethod
    def load(cls, path=None):
        """
        Load a configuration from a file.  If no file is provided,
        defaults to `asv.conf.json`.
        """
        if not path:
            path = "asv.conf.json"

        if not os.path.isfile(path):
            raise util.UserError("Config file {0} not found.".format(path))

        d = util.load_json(path, cls.api_version, js_comments=True)
        try:
            return cls.from_json(d)
        except ValueError:
            raise util.UserError(
                "No repo specified in {0} config file.".format(path))

    @classmethod
    def from_json(cls, d):
        if 'wheel_cache_size' in d:
            log.warning("`wheel_cache_size` has been renamed to `build_cache_size`. Update your `asv.conf.json` accordingly.")
            d.setdefault('build_cache_size', d['wheel_cache_size'])

        conf = cls()
        conf.__dict__.update(d)

        if not getattr(conf, "repo", None):
            raise util.UserError(
                "No repo specified in config file.")

        if not getattr(conf, "branches", [None]):
            # If 'branches' attribute is present, at least some must
            # be listed.
            raise util.UserError(
                "No branches specified in config file.")

        return conf
