#!/usr/bin/env python3
# coding: utf-8

"""
MFBN (Multilevel framework for bipartite networks)

Copyright (C) 2017 Alan Valejo <alanvalejo@gmail.com> All rights reserved

This program comes with ABSOLUTELY NO WARRANTY. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS
WITH YOU.

Owner or contributors are not liable for any direct, indirect, incidental, special, exemplary, or consequential
damages, (such as loss of data or profits, and others) arising in any way out of the use of this software,
even if advised of the possibility of such damage.

This program is free software and distributed in the hope that it will be useful: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with this program. If not,
see http://www.gnu.org/licenses/.

Giving credit to the author by citing the papers.
"""

import operator
import numpy
import random
import math
import collections

from random import sample
from igraph import Graph

__maintainer__ = 'Alan Valejo'
__email__ = 'alanvalejo@gmail.com'
__author__ = 'Alan Valejo'
__credits__ = ['Alan Valejo']
__homepage__ = 'https://www.alanvalejo.com.br'
__license__ = 'GNU.GPL.v3'
__docformat__ = 'markdown en'
__version__ = '0.1'
__date__ = '2019-08-08'


def load_ncol(filename):
    """
    Load ncol npartite graph and generate special attributes
    """

    data = numpy.loadtxt(filename, skiprows=0, dtype=str)
    dict_edges = dict()
    for row in data:
        if len(row) == 3:
            dict_edges[(int(row[0]), int(row[1]))] = float(row[2])
        else:
            dict_edges[(int(row[0]), int(row[1]))] = 1

    edges, weights = list(zip(*dict_edges.items()))
    return edges, weights


class MGraph(Graph):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def load(self, network_filename, vertices, filename_type='ncol', type_filename=None):
        """
        filename_type: ncol, arff
        """

        edges, weights = None, None
        if filename_type == 'ncol':
            edges, weights = load_ncol(network_filename)

        self.add_vertices(sum(vertices))
        self.add_edges(edges)
        self['adjlist'] = list(map(set, self.get_adjlist()))
        self['vertices'] = vertices
        self['layers'] = len(vertices)
        self['level'] = [0] * self['layers']
        self['similarity'] = None
        self.es['weight'] = weights
        types = []
        for layer in range(self['layers']):
            types += [layer] * vertices[layer]
        self.vs['type'] = types
        for v in self.vs():
            v['weight'] = 1
            v['source'] = [v.index]
            v['name'] = [v.index]
            v['predecessor'] = [v.index]
            v['successor'] = [None]

        self['vertices_by_type'] = []
        for layer in range(self['layers']):
            self['vertices_by_type'].append(self.vs.select(type=layer).indices)

        # Not allow direct graphs
        if self.is_directed():
            self.to_undirected(combine_edges=None)

    def contract(self, matching):
        """
        Create coarse graph from matching of groups
        """

        # Contract vertices: Referencing the original graph of the coarse graph
        types = []
        weights = []
        sources = []
        predecessors = []
        matching = numpy.array(matching)
        uniqid = 0
        clusters = numpy.unique(matching)
        for cluster_id in clusters:
            vertices = numpy.where(matching == cluster_id)[0]
            weight = 0
            if len(vertices) > 0:
                source = []
                predecessor = []
                for vertex in vertices:
                    self.vs[vertex]['successor'] = uniqid
                    weight += self.vs[vertex]['weight']
                    source.extend(self.vs[vertex]['source'])
                    predecessor.append(vertex)
                weights.append(weight)
                types.append(self.vs[vertices[0]]['type'])
                sources.append(source)
                predecessors.append(predecessor)
                uniqid += 1

        # Create coarsened version
        coarse = MGraph()
        coarse.add_vertices(uniqid)
        coarse.vs['type'] = types
        coarse.vs['weight'] = weights
        coarse.vs['name'] = range(coarse.vcount())
        coarse.vs['successor'] = [None] * coarse.vcount()
        coarse.vs['source'] = sources
        coarse.vs['predecessor'] = predecessors
        coarse['layers'] = self['layers']
        coarse['similarity'] = None
        coarse['vertices'] = []

        coarse['vertices_by_type'] = []
        for layer in range(self['layers']):
            coarse['vertices_by_type'].append(coarse.vs.select(type=layer).indices)
            coarse['vertices'].append(len(coarse['vertices_by_type'][layer]))

        # Contract edges
        dict_edges = dict()
        for edge in self.es():
            v_successor = self.vs[edge.tuple[0]]['successor']
            u_successor = self.vs[edge.tuple[1]]['successor']

            # Add edge in coarsened graph
            if v_successor < u_successor:
                dict_edges[(v_successor, u_successor)] = dict_edges.get((v_successor, u_successor), 0) + edge['weight']
            else:
                dict_edges[(u_successor, v_successor)] = dict_edges.get((u_successor, v_successor), 0) + edge['weight']

        if len(dict_edges) > 0:
            edges, weights = list(zip(*dict_edges.items()))
            coarse.add_edges(edges)
            coarse.es['weight'] = weights
            coarse['adjlist'] = list(map(set, coarse.get_adjlist()))

        return coarse

    def gmb(self, vertices=None, reduction_factor=0.5, reverse=True):
        """
        Matches are restricted between vertices that are not adjacent
        but are only allowed to match with neighbors of its neighbors,
        i.e. two-hopes neighborhood
        """

        matching = numpy.array([-1] * self.vcount())
        matching[vertices] = vertices

        # Search two-hopes neighborhood for each vertex in selected layer
        dict_edges = dict()
        visited = [0] * self.vcount()
        for vertex in vertices:
            neighborhood = self.neighborhood(vertices=vertex, order=2)
            twohops = neighborhood[(len(self['adjlist'][vertex]) + 1):]
            for twohop in twohops:
                if visited[twohop] == 1:
                    continue
                dict_edges[(vertex, twohop)] = self['similarity'](vertex, twohop)
            visited[vertex] = 1

        # Select promising matches or pair of vertices
        visited = [0] * self.vcount()
        edges = sorted(dict_edges.items(), key=operator.itemgetter(1), reverse=reverse)
        merge_count = int(reduction_factor * len(vertices))
        for edge, value in edges:
            vertex = edge[0]
            neighbor = edge[1]
            if (visited[vertex] != 1) and (visited[neighbor] != 1):
                matching[neighbor] = vertex
                matching[vertex] = vertex
                visited[neighbor] = 1
                visited[vertex] = 1
                merge_count -= 1
            if merge_count == 0:
                break

        return matching

    def rgmb(self, matching, vertices=None, reduction_factor=0.5, seed_priority='random', reverse=True):
        """
        Matches are restricted between vertices that are not adjacent
        but are only allowed to match with neighbors of its neighbors,
        i.e. two-hopes neighborhood. This version use a random seed.
        """

        matching = numpy.array([-1] * self.vcount())
        matching[vertices] = vertices

        # Select seed set expansion
        if seed_priority == 'strength':
            vertices_score = numpy.array(self.strength(vertices, weights='weight'))
            dictionary = dict(zip(vertices, vertices_score))
            vertices_id = sorted(dictionary, key=dictionary.__getitem__, reverse=reverse)
        if seed_priority == 'degree':
            vertices_score = numpy.array(self.degree(vertices))
            dictionary = dict(zip(vertices, vertices_score))
            vertices_id = sorted(dictionary, key=dictionary.__getitem__, reverse=reverse)
        if seed_priority == 'random':
            vertices_id = vertices
            vertices_id = random.sample(vertices_id, len(vertices_id))

        # Find the matching
        visited = [0] * self.vcount()
        index = 0
        merge_count = int(reduction_factor * len(vertices))
        while merge_count > 0 and index < len(vertices):
            # Randomly select a vertex v of V
            vertex = vertices_id[index]
            if visited[vertex] == 1:
                index += 1
                continue
            # Select the edge (v, u) of E which maximum score
            # Tow hopes restriction: It ensures that the match only occurs
            # between vertices of the same type
            neighborhood = self.neighborhood(vertices=vertex, order=2)
            twohops = neighborhood[(len(self['adjlist'][vertex]) + 1):]
            # twohops = set((twohop for onehop in self['adjlist'][vertex] for twohop in self['adjlist'][onehop])) -
            # set([vertex])
            _max = 0.0
            neighbor = vertex
            for twohop in twohops:
                if visited[twohop] == 1:
                    continue
                # Calling a function of a module from a string
                score = self['similarity'](vertex, twohop)
                if score > _max:
                    _max = score
                    neighbor = twohop
            matching[neighbor] = vertex
            matching[vertex] = vertex
            visited[neighbor] = 1
            visited[vertex] = 1
            merge_count -= 1
            index += 1

        return matching

    def rm(self, vcount, reduction_factor=0.5):
        """
        Random Matching: Select a maximal matching using a
        randomized algorithm
        """

        matching = numpy.array([-1] * vcount)
        merge_count = int(reduction_factor * self.vcount())
        self.get_random_edges(merge_count, matching)
        return matching

    def lem(self, vcount, reduction_factor=0.5):
        """
        Heavy Light Matching: Search for a minimal matching using the
        weights of the edges of the graph.
        """

        matching = numpy.array([-1] * vcount)
        merge_count = int(reduction_factor * self.vcount())
        self.get_sorted_edges(merge_count, matching, reverse=False)
        return matching

    def hem(self, vcount, reduction_factor=0.5):
        """
        Heavy Edge Matching: Search for a maximal matching using the
        weights of the edges of the graph.
        """

        matching = numpy.array([-1] * vcount)
        merge_count = int(reduction_factor * self.vcount())
        self.get_sorted_edges(merge_count, matching, reverse=True)
        return matching

    def get_random_edges(self, merge_count, matching):
        """
        Return a random independent edge set in a graph, i.e., is a set
        of edges without common vertices random selected.
        """

        visited = [0] * self.vcount()
        edges = sample(self.es(), self.ecount())
        for edge in edges:
            if (visited[edge.tuple[0]] == 0) and (visited[edge.tuple[1]] == 0):
                u = self.vs[edge.tuple[0]]['name']
                v = self.vs[edge.tuple[1]]['name']
                matching[u] = u
                matching[v] = u
                visited[edge.tuple[1]] = 1
                visited[edge.tuple[0]] = 1
                merge_count -= 1
            if merge_count == 0:
                break

    def get_sorted_edges(self, merge_count, matching, reverse=True):
        """
        Search for a maximal matching using the weights of the edges of
        the graph. The aim is to find a maximal matching of the graph that
        minimizes the cut.
        """

        visited = [0] * self.vcount()
        edges = sorted(self.es(), key=lambda edge: edge['weight'], reverse=reverse)
        for edge in edges:
            if (visited[edge.tuple[0]] == 0) and (visited[edge.tuple[1]] == 0):
                u = self.vs[edge.tuple[0]]['name']
                v = self.vs[edge.tuple[1]]['name']
                matching[u] = u
                matching[v] = u
                visited[edge.tuple[1]] = 1
                visited[edge.tuple[0]] = 1
                merge_count -= 1
            if merge_count == 0:
                break

    def weighted_one_mode_projection(self, vertices):
        """
        Application of a one-mode projection to a bipartite network generates
        two unipartite networks, one for each layer, so that vertices with
        common neighbors are connected by edges in their respective projection.
        """

        graph = MGraph()
        graph.add_vertices(vertices)
        graph.vs['name'] = self.vs[vertices]['name']
        name_to_id = dict(zip(vertices, range(graph.vcount())))

        dict_edges = dict()
        visited = [0] * self.vcount()
        for vertex in vertices:
            neighborhood = self.neighborhood(vertices=vertex, order=2)
            twohops = neighborhood[(len(self['adjlist'][vertex]) + 1):]
            for twohop in twohops:
                if visited[twohop] == 1:
                    continue
                dict_edges[(name_to_id[vertex], name_to_id[twohop])] = self['similarity'](vertex, twohop)
            visited[vertex] = 1

        if len(dict_edges) > 0:
            edges, weights = list(zip(*dict_edges.items()))
            graph.add_edges(edges)
            graph.es['weight'] = weights

        return graph

    def mlpb(self, vertices=None, seed_priority='strength', reduction_factor=0.5, itr=10, tolerance=0.05,
             upper_bound=0.2, n=None, global_min_vertices=None, reverse=True):

        """ Matching via weight-constrained label propagation and neighborhood. """

        matching = numpy.array([-1] * self.vcount())
        matching[vertices] = vertices

        if global_min_vertices:
            min_vertices = global_min_vertices
        else:
            min_vertices = int((1 - reduction_factor) * len(vertices))
        if min_vertices < 1:
            min_vertices = 1

        max_size = int(math.ceil(((1.0 + upper_bound) * n) / min_vertices))
        number_of_vertices = len(vertices)
        weight_of_sv = self.vs['weight']
        label_dict = dict(zip(vertices, vertices))
        twohops_dict = collections.defaultdict(float)
        similarity_dict = collections.defaultdict(float)

        # Select seed set expansion: case of strength or degree seed
        if seed_priority == 'strength':
            vertices_score = numpy.array(self.strength(vertices, weights='weight'))
            dictionary = dict(zip(vertices, vertices_score))
            vertices_id = sorted(dictionary, key=dictionary.__getitem__, reverse=reverse)
        if seed_priority == 'degree':
            vertices_score = numpy.array(self.degree(vertices))
            dictionary = dict(zip(vertices, vertices_score))
            vertices_id = sorted(dictionary, key=dictionary.__getitem__, reverse=reverse)

        tolerance = tolerance * len(vertices)
        swap = tolerance + 1

        while (tolerance < swap) and itr:
            swap = 0
            itr -= 1

            # Select seed set expansion: case of random seed
            if seed_priority == 'random':
                vertices_id = vertices
                vertices_id = random.sample(vertices_id, len(vertices_id))

            for vertex in vertices_id:

                if self.degree(vertex) == 0:
                    continue

                # Tow hopes restriction: It ensures that the match only occurs
                # between vertices of the same type
                if not twohops_dict.get(vertex, False):
                    neighborhood = self.neighborhood(vertices=vertex, order=2)
                    twohops_dict[vertex] = neighborhood[(len(self['adjlist'][vertex]) + 1):]

                # Update neighborhood edge density
                Q = collections.defaultdict(float)
                for neighbor in twohops_dict[vertex]:
                    if weight_of_sv[label_dict[neighbor]] + self.vs[vertex]['weight'] <= max_size:
                        if vertex < neighbor:
                            u, v = vertex, neighbor
                        else:
                            u, v = neighbor, vertex
                        if not similarity_dict.get((u, v), False):
                            similarity_dict[(u, v)] = self['similarity'](u, v)
                        if similarity_dict[(u, v)] > 0.0:
                            Q[label_dict[neighbor]] += similarity_dict[(u, v)]

                total_similarity = sum(Q.values())
                for li in Q.keys():
                    Q[li] -= (total_similarity - Q[li])

                if Q:
                    # Select the dominant label
                    dominant_label = max(Q.items(), key=operator.itemgetter(1))[0]
                    prev_label = label_dict[vertex]
                    # If a dominant label was fund, match them together
                    # and update data structures
                    if dominant_label != prev_label:
                        swap += 1
                        # Update vertex label
                        label_dict[vertex] = dominant_label
                        # Update the super-vertex weight
                        weight_of_sv[prev_label] -= self.vs[vertex]['weight']
                        weight_of_sv[dominant_label] += self.vs[vertex]['weight']
                        # Verify the size-constraint restriction
                        if weight_of_sv[prev_label] == 0:
                            number_of_vertices -= 1
                        if number_of_vertices <= min_vertices:
                            tolerance = swap
                            break

        for key, value in label_dict.items():
            matching[key] = value

        return matching
