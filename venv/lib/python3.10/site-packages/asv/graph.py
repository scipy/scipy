# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import traceback

import six
from six.moves import xrange

from . import util
from . import step_detect

from .util import is_na, mean_na, geom_mean_na


# This is the maximum number of points to include in summary graphs.
# It is based on the number of pixels in the summary graph display on
# a recent Retina MacBook Pro (3840 pixels across the screen, divided
# by 5 summaries across, divided by 2 for good measure and to account
# for width of the line).
RESAMPLED_POINTS = (3840 / 5 / 2)


class GraphSet(object):
    """Manage multiple `Graph`"""

    def __init__(self):
        self._graphs = {}
        self._groups = {}
        super(GraphSet, self).__init__()

    def get_graph(self, benchmark_name, params):
        graph = Graph(benchmark_name, params)
        if graph.path not in self._graphs:
            self._graphs[graph.path] = graph
            self._groups.setdefault(benchmark_name, []).append(graph)
        return self._graphs[graph.path]

    def get_graph_group(self, benchmark_name):
        return self._groups.get(benchmark_name, [])

    def get_params(self):
        """Return all params used in graphs and their corresponding values set"""
        params = {}
        for graph in six.itervalues(self._graphs):
            for key, value in six.iteritems(graph.params):
                params.setdefault(key, set())
                if value:
                    params[key].add(value)
        return params

    def detect_steps(self, pool=None, dots=None):
        for graph in six.itervalues(self._graphs):
            graph.detect_steps(pool)
            if dots is not None and pool is None:
                dots()

        # Wait for results to compute
        for graph in six.itervalues(self._graphs):
            graph.get_steps()
            if dots is not None and pool is not None:
                dots()

    def get_summary_graphs(self, dots=None):
        for graphs in six.itervalues(self._groups):
            yield make_summary_graph(graphs)
            if dots is not None:
                dots()

    def save(self, html_dir, dots=None):
        for graph in six.itervalues(self._graphs):
            graph.save(html_dir)
            if dots is not None:
                dots()

    def __iter__(self):
        return six.iteritems(self._graphs)

    def __len__(self):
        return len(self._graphs)


class Graph(object):
    """
    Manages a single "line" in the resulting plots for the front end.

    Unlike "results", which contain the timings for a single commit,
    these contain the timings for a single benchmark.
    """
    def __init__(self, benchmark_name, params):
        """
        Initially the graph contains no data.  It must be added using
        multiple calls to `add_data_point`.

        Parameters
        ----------
        benchmark_name : str
            A unique string to identify the benchmark, and display in
            the frontend.

        params : dict of str -> str
            A dictionary of parameters describing the benchmark.

        """
        self.benchmark_name = benchmark_name
        self.params = params
        self.data_points = {}
        self.data_weights = {}

        self.path = self.get_file_path(self.params, benchmark_name)
        self.n_series = None
        self.scalar_series = True
        self._steps = None

    @classmethod
    def get_file_path(cls, params, benchmark_name):
        """
        Get a file path understood by the JS frontend, corresponding
        on the given parameters and benchmark_name.

        The implementation must match asv.js:graph_to_path
        """
        parts = ['graphs']
        l = list(six.iteritems(params))
        l.sort()
        for key, val in l:
            if val is None:
                part = '{0}-null'.format(key)
            elif val:
                part = '{0}-{1}'.format(key, val)
            else:
                part = '{0}'.format(key)
            parts.append(util.sanitize_filename(part))
        parts.append(util.sanitize_filename(benchmark_name))
        return os.path.join(*parts)

    def add_data_point(self, revision, value, weight=None):
        """
        Add a data point to the graph.

        Parameters
        ----------
        revision : int
            An integer value representing the commit revision in the commit log

        value : float or list
            The value(s) to plot in the benchmark.

        weight : float or list or None
            Weights corresponding to the values.
            Missing estimates are indicated with None.

        """
        self.data_points.setdefault(revision, [])
        self.data_weights.setdefault(revision, [])
        if not is_na(value):
            if not hasattr(value, '__len__'):
                value = [value]
                weight = [weight]
            else:
                self.scalar_series = False

            if self.n_series is None:
                self.n_series = len(value)
            elif len(value) != self.n_series:
                raise ValueError("Mismatching number of data series in graph")

            if weight is None:
                weight = [None]*len(value)

            self.data_points[revision].append(value)
            self.data_weights[revision].append(weight)

    def get_data(self):
        """
        Get the sorted and reduced data and weights.
        """

        if self.n_series is None:
            # No non-null data points
            self.n_series = 1

        def mean_axis0(v):
            if not v:
                return [None]*self.n_series
            return [mean_na(x[j] for x in v)
                    for j in xrange(self.n_series)]

        # Average data over commit log
        val = []
        for k in six.iterkeys(self.data_points):
            v = mean_axis0(self.data_points[k])
            w = mean_axis0(self.data_weights[k])
            val.append((k, v, w))
            del v, w
        val.sort()

        # Discard missing data at edges
        i = 0
        for i in xrange(len(val)):
            if any(not is_na(v) for v in val[i][1]):
                break
        else:
            i = len(val)

        j = i
        for j in xrange(len(val) - 1, i, -1):
            if any(not is_na(v) for v in val[j][1]):
                break

        val = val[i:j+1]

        # Single-element series
        if self.scalar_series:
            val = [(k, v[0], w[0]) for k, v, w in val]

        return val

    def save(self, html_dir):
        """
        Save the graph to a .json file used by the frontend.

        Parameters
        ----------
        html_dir : str
            The root of the HTML tree.
        """
        filename = os.path.join(html_dir, self.path + ".json")

        # Drop weights
        val = [v[:2] for v in self.get_data()]

        util.write_json(filename, val, compact=True)

    def detect_steps(self, pool=None):
        """
        Run step detection algorithm on the graph data.

        Afterward, the results can be obtained via get_steps()

        Parameters
        ----------
        pool : multiprocessing.Pool, optional
            Pool to use for asynchronous jobs.
            If not given, run in serial.

        """
        if self._steps is not None:
            # Already computed
            return

        val = self.get_data()

        if not val:
            # Nothing to compute
            self._steps = [[]]*self.n_series
            return

        if self.scalar_series:
            items = [val]
        else:
            items = [[(v[0], v[1][j], v[2][j]) for v in val] for j in range(self.n_series)]

        if pool is None:
            self._steps = [_compute_graph_steps(item, reraise=False) for item in items]
        else:
            self._steps = [pool.apply_async(_compute_graph_steps, (item,)) for item in items]

    def get_steps(self):
        """
        Return results from step detection.

        Returns
        -------
        steps : list of (left, right, val, min, err)
            Result of fitting a piecewise function to the graph.
            Missing data points do not necessarily belong in any piece.
            The values are: `left` (inclusive) and `right` (exclusive) specify
            a revision interval, `val` the median value in the interval, `min`
            the minimum value in the interval, and `err` the mean deviation from
            the median.

        """
        if self._steps is None:
            self.detect_steps()

        for j, item in enumerate(self._steps):
            if not isinstance(item, list):
                self._steps[j] = item.get()

        if self.scalar_series:
            return self._steps[0]
        else:
            return self._steps


def _compute_graph_steps(data, reraise=True):
    try:
        x = [d[0] for d in data]
        y = [d[1] for d in data]
        w = [d[2] for d in data]

        steps = step_detect.detect_steps(y, w)
        new_steps = []

        for left, right, cur_val, cur_min, cur_err in steps:
            new_steps.append((x[left], x[right-1] + 1, cur_val, cur_min, cur_err))

        return new_steps
    except BaseException as exc:
        if reraise:
            raise util.ParallelFailure(str(exc), exc.__class__, traceback.format_exc())
        else:
            raise


def make_summary_graph(graphs):
    val, n_series = _combine_graph_data(graphs)

    # Given multiple input series
    #
    #     val = [(x_0, (y[0,0], y[0,1], ..., y[0,n])),
    #            (x_1, (y[1,0], y[1,1], ..., y[1,n])),
    #            ... ]
    #
    # calculate summary data series
    #
    #     z = geom_mean(y, axis=1)
    #
    # Missing data in y is filled by the previous non-null
    # values (or the first non-null value, for nulls at the
    # beginning), to avoid meaningless jumps in the result.
    # Data points missing from all series are not filled.

    # Find first non-null values
    first_values = [None]*n_series
    for k, v in val:
        for j, x in enumerate(v):
            if first_values[j] is None and not is_na(x):
                first_values[j] = x
        if not any(is_na(x) for x in first_values):
            break

    first_values = [fv if fv is not None else 1.0
                    for fv in first_values]

    # Compute geom mean of filled series
    last_values = [None]*n_series
    new_val = []
    for k, v in val:
        # Fill missing data, unless it's missing from all
        # parameter combinations
        cur_vals = []
        if any(not is_na(x) for x in v):
            for j, x in enumerate(v):
                if is_na(x):
                    if last_values[j] is not None:
                        x = last_values[j]
                    else:
                        x = first_values[j]
                else:
                    last_values[j] = x

                cur_vals.append(x)

        # Geometric mean of values
        v = geom_mean_na(cur_vals)
        new_val.append((k, v))

    val = new_val

    # Resample
    val = resample_data(val)

    # Return as a graph
    graph = Graph(graphs[0].benchmark_name, {'summary': ''})
    for x, y in val:
        graph.add_data_point(x, y)
    return graph


def _combine_graph_data(graphs):
    """
    Concatenate data series from multiple graphs into a single series

    Returns
    -------
    val
        List of the form [(x_0, [y_00, y_01, ...]), (x_1, ...)]
        where the y-values are obtained by concatenating the data
        series. When some of the graphs do not have data for a given x-value,
        the missing data is indicated by None values.
    n_series
        Number of data series in output. Equal to the sum of n_series
        of the input graphs.

    """

    all_data = {}
    n_series = sum(graph.n_series if graph.n_series else 1
                   for graph in graphs)

    template = [None]*n_series
    pos = 0

    for graph in graphs:
        series = graph.get_data()
        for k, v, dv in series:
            prev = all_data.get(k, template)
            if graph.scalar_series:
                v = [v]
            all_data[k] = prev[:pos] + v + prev[pos+graph.n_series:]
        pos += graph.n_series

    val = list(all_data.items())
    val.sort()

    return val, n_series


def resample_data(val):
    if len(val) < RESAMPLED_POINTS:
        return val

    min_revision = min(x[0] for x in val)
    max_revision = max(x[0] for x in val)
    step_size = int((max_revision - min_revision) / RESAMPLED_POINTS)

    if step_size == 0:
        step_size = max_revision - min_revision + 1

    new_val = []
    j = 0
    for i in xrange(min_revision + step_size, max_revision + step_size, step_size):
        chunk = []
        while j < len(val) and val[j][0] < i:
            chunk.append(val[j][1])
            j += 1
        if len(chunk):
            new_val.append((i, mean_na(chunk)))
    return new_val
