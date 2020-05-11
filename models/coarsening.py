#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coarsening

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

__maintainer__ = 'Alan Valejo'
__email__ = 'alanvalejo@gmail.com'
__author__ = 'Alan Valejo'
__credits__ = ['Alan Valejo']
__homepage__ = 'https://www.alanvalejo.com.br'
__license__ = 'GNU.GPL.v3'
__docformat__ = 'markdown en'
__version__ = '0.1'
__date__ = '2019-08-08'

import sys
import numpy
import multiprocessing as mp

from models.timing import Timing
from models.similarity import Similarity


def modified_starmap_async(function, kwargs):
    return function(**kwargs)


class Coarsening:

    def __init__(self, source_graph, **kwargs):

        prop_defaults = {
            'reduction_factor': None, 'max_levels': None, 'matching': None,
            'similarity': None, 'itr': None, 'upper_bound': None, 'seed_priority': None,
            'global_min_vertices': None, 'tolerance': None, 'reverse': None, 'threads': 1
        }

        self.__dict__.update(prop_defaults)
        self.__dict__.update(kwargs)

        self.source_graph = source_graph
        self.hierarchy_graphs = []
        self.hierarchy_levels = []

        # Create default values for optional parameters
        if self.reduction_factor is None:
            self.reduction_factor = [0.5] * self.source_graph['layers']
        if self.max_levels is None:
            self.max_levels = [3] * self.source_graph['layers']
        if self.matching is None:
            self.matching = ['rgmb'] * self.source_graph['layers']
        if self.similarity is None:
            self.similarity = ['common_neighbors'] * self.source_graph['layers']
        if self.itr is None:
            self.itr = [10] * self.source_graph['layers']
        if self.upper_bound is None:
            self.upper_bound = [2.0] * self.source_graph['layers']
        if self.global_min_vertices is None:
            self.global_min_vertices = [None] * self.source_graph['layers']
        if self.tolerance is None:
            self.tolerance = [0.01] * self.source_graph['layers']

        # Validation of list values
        if len(self.reduction_factor) == 1:
            self.reduction_factor = [self.reduction_factor[0]] * self.source_graph['layers']
        if len(self.max_levels) == 1:
            self.max_levels = [self.max_levels[0]] * self.source_graph['layers']
        if len(self.matching) == 1:
            self.matching = [self.matching[0]] * self.source_graph['layers']
        if len(self.similarity) == 1:
            self.similarity = [self.similarity[0]] * self.source_graph['layers']
        if len(self.itr) == 1:
            self.itr = [self.itr[0]] * self.source_graph['layers']
        if len(self.upper_bound) == 1:
            self.upper_bound = [self.upper_bound[0]] * self.source_graph['layers']
        if len(self.global_min_vertices) == 1:
            self.global_min_vertices = [self.global_min_vertices[0]] * self.source_graph['layers']
        if len(self.tolerance) == 1:
            self.tolerance = [self.tolerance[0]] * self.source_graph['layers']
        if len(self.seed_priority) == 1:
            self.seed_priority = [self.seed_priority[0]] * self.source_graph['layers']
        if len(self.reverse) == 1:
            self.reverse = [self.reverse[0]] * self.source_graph['layers']

        # Parameters dimension validation
        if self.source_graph['layers'] != len(self.reduction_factor):
            print('Number of layers and reduction_factor do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.max_levels):
            print('Number of layers and max_levels do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.matching):
            print('Number of layers and matching do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.similarity):
            print('Number of layers and similarity do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.itr):
            print('Number of layers and itr do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.upper_bound):
            print('Number of layers and upper_bound do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.global_min_vertices):
            print('Number of layers and global_min_vertices do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.tolerance):
            print('Number of layers and tolerance do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.seed_priority):
            print('Number of layers and seed_priority do not match.')
            sys.exit(1)
        if self.source_graph['layers'] != len(self.reverse):
            print('Number of layers and reverse do not match.')
            sys.exit(1)
        if self.threads > mp.cpu_count():
            print('Number of defined threads (' + str(self.threads) + ') cannot be greater than the real number'
                                                                      'of cors (' + str(mp.cpu_count()) + ').')
            sys.exit(1)

        # Matching method validation
        valid_matching = ['rgmb', 'gmb', 'mlpb', 'hem', 'lem', 'rm']
        for index, matching in enumerate(self.matching):
            matching = matching.lower()
            if matching not in valid_matching:
                print('Matching ' + matching + ' method is invalid.')
                sys.exit(1)
            self.matching[index] = matching

        # Seed priority validation
        valid_seed_priority = ['strength', 'degree', 'random']
        for index, seed_priority in enumerate(self.seed_priority):
            seed_priority = seed_priority.lower()
            if seed_priority not in valid_seed_priority:
                print('Seed priotiry ' + seed_priority + ' is invalid.')
                sys.exit(1)
            self.seed_priority[index] = seed_priority

        # Reverse validation
        for index, reverse in enumerate(self.reverse):
            if reverse.lower() in ('yes', 'true', 't', 'y', '1'):
                self.reverse[index] = True
            elif reverse.lower() in ('no', 'false', 'f', 'n', '0'):
                self.reverse[index] = False
            else:
                print('Boolean value expected in -rv.')
                sys.exit(1)

        # Similarity measure validation
        valid_similarity = [
            'common_neighbors', 'weighted_common_neighbors',
            'salton', 'preferential_attachment', 'jaccard', 'weighted_jaccard',
            'adamic_adar', 'resource_allocation', 'sorensen', 'hub_promoted',
            'hub_depressed', 'leicht_holme_newman'
        ]
        for index, similarity in enumerate(self.similarity):
            similarity = similarity.lower()
            if similarity not in valid_similarity:
                print('Similarity ' + similarity + ' misure is unvalid.')
                sys.exit(1)
            self.similarity[index] = similarity

        for layer in range(self.source_graph['layers']):
            if self.matching[layer] in ['rgmb', 'gmb', 'hem', 'lem', 'rm']:
                if self.global_min_vertices[layer] is not None:
                    self.global_min_vertices[layer] = None
                    text = 'Matching method ' + self.matching[layer]
                    text += ' (setted in layer '
                    text += str(layer) + ') does not accept -gmv parameter.'
                    print(text)
                if self.reduction_factor[layer] > 0.5:
                    self.reduction_factor[layer] = 0.5
                    text = 'Matching method ' + self.matching[layer]
                    text += ' (setted in layer '
                    text += str(layer) + ') does not accept -rf > 0.5.'
                    print(text)

    def run(self):

        graph = self.source_graph.copy()
        while True:

            level = graph['level']
            contract = False

            args = []
            for layer in range(graph['layers']):
                do_matching = True
                if self.global_min_vertices[layer] is None and level[layer] >= self.max_levels[layer]:
                    do_matching = False
                elif self.global_min_vertices[layer] and graph['vertices'][layer] <= self.global_min_vertices[layer]:
                    do_matching = False

                if do_matching:

                    contract = True
                    level[layer] += 1

                    graph['similarity'] = getattr(Similarity(graph, graph['adjlist']), self.similarity[layer])

                    kwargs = dict(reduction_factor=self.reduction_factor[layer])

                    if self.matching[layer] in ['mlpb', 'gmb', 'rgmb']:
                        kwargs['vertices'] = graph['vertices_by_type'][layer]
                        kwargs['reverse'] = self.reverse[layer]
                    if self.matching[layer] in ['mlpb', 'rgmb']:
                        kwargs['seed_priority'] = self.seed_priority[layer]
                    if self.matching[layer] in ['mlpb']:
                        kwargs['upper_bound'] = self.upper_bound[layer]
                        kwargs['n'] = self.source_graph['vertices'][layer]
                        kwargs['global_min_vertices'] = self.global_min_vertices[layer]
                        kwargs['tolerance'] = self.tolerance[layer]
                        kwargs['itr'] = self.itr[layer]

                    if self.matching[layer] in ['hem', 'lem', 'rm']:
                        one_mode_graph = graph.weighted_one_mode_projection(graph['vertices_by_type'][layer])
                        matching_function = getattr(one_mode_graph, self.matching[layer])
                    else:
                        matching_function = getattr(graph, self.matching[layer])

                    # Create a args for the engine multiprocessing.pool
                    args.append([(matching_function, kwargs)])

            if contract:
                # Create pools
                pool = mp.Pool(processes=self.threads)
                processes = []
                for arg in args:
                    processes.append(pool.starmap_async(modified_starmap_async, arg))

                # Merge chunked solutions
                matching = numpy.arange(graph.vcount())
                for process in processes:
                    result = process.get()[0]
                    vertices = numpy.where(result > -1)[0]
                    matching[vertices] = result[vertices]

                # Close processes
                pool.close()
                pool.join()

                # Contract current graph using the matching
                coarsened_graph = graph.contract(matching)
                coarsened_graph['level'] = level

                if coarsened_graph.vcount() == graph.vcount():
                    break

                self.hierarchy_graphs.append(coarsened_graph)
                self.hierarchy_levels.append(level[:])
                graph = coarsened_graph
            else:
                break
