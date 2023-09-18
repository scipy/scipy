'''
Minimal directed graph replacement for networkx.DiGraph

This has the sole advantage of being a standalone file that doesn't bring any
dependency with it.
'''

class DiGraph(object):
    def __init__(self):
        # adjacency[i][j] = True means j is a successor of i
        self._adjacency = {}
        self._edges = {}

    def successors(self, node):
        return (n for n in self._adjacency[node])

    def predecessors(self, node):
        return (k for k, v in self._adjacency.items() if node in v)

    def add_node(self, node):
        self._adjacency.setdefault(node, set())

    def add_edge(self, src, dest, **props):
        self.add_node(dest)
        self._adjacency.setdefault(src, set()).add(dest)
        self._edges[(src, dest)] = props

    @property
    def edges(self):
        return self._edges

    def has_edge(self, src, dest):
        return dest in self._adjacency[src]

    def remove_edge(self, src, dest):
        self._adjacency[src].remove(dest)
        del self._edges[(src, dest)]

    def __len__(self):
        return len(self._adjacency)

    def __iter__(self):
        return iter(self._adjacency.keys())

    def __contains__(self, value):
        return value in self._adjacency

    def __getitem__(self, node):
        return self._adjacency[node]

class Unfeasible(RuntimeError):
    pass

def has_path(graph, src, dest):
    visited = set()
    worklist = [src]
    while worklist:
        current = worklist.pop()
        if current in visited:
            continue
        visited.add(current)
        if dest in graph.successors(current):
            return True
        worklist.extend(graph.successors(current))
    return False

# Copied verbatim from NetworkX 2.6.1
#
# NetworkX is distributed with the 3-clause BSD license.
#
# ::
#
#   Copyright (C) 2004-2021, NetworkX Developers
#   Aric Hagberg <hagberg@lanl.gov>
#   Dan Schult <dschult@colgate.edu>
#   Pieter Swart <swart@lanl.gov>
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the NetworkX Developers nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

def _all_simple_paths_graph(G, source, targets, cutoff):
    visited = dict.fromkeys([source])
    stack = [iter(G[source])]
    while stack:
        children = stack[-1]
        child = next(children, None)
        if child is None:
            stack.pop()
            visited.popitem()
        elif len(visited) < cutoff:
            if child in visited:
                continue
            if child in targets:
                yield list(visited) + [child]
            visited[child] = None
            if targets - set(visited.keys()):  # expand stack until find all targets
                stack.append(iter(G[child]))
            else:
                visited.popitem()  # maybe other ways to child
        else:  # len(visited) == cutoff:
            for target in (targets & (set(children) | {child})) - set(visited.keys()):
                yield list(visited) + [target]
            stack.pop()
            visited.popitem()

def all_simple_paths(graph, src, target):
    return _all_simple_paths_graph(graph, src, {target},  len(graph) - 1)
