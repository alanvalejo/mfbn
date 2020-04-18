#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coarsening
=====================================================

Copyright (C) 2017 Alan Valejo <alanvalejo@gmail.com> All rights reserved

In coarsening strategy a sequence (or hierarchy) of smaller networks is constructed from the original network, such
that $|V_0| > |V_1| > ... > |V_N|$. Such a hierarchy represents the network on multiple scales.

The original implementation (paper version) is deprecated. This software is a new version, more robust and fast.
I.e., there may be divergences between this version and the original algorithm. If you looking for the original
version used in the paper don't hesitate to contact the authors.

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

Giving credit to the author by citing the papers [1,2,3,4,5]

.. [1] Valejo, Alan and Faleiros, T. P. and Oliveira, Maria C. F. and Lopes, A. A., A coarsening method for bipartite
networks via weight-constrained label propagation, in ACM Computing Surveys, 2019,
doi: https://doi.org/10.1016/j.knosys.2020.105678

.. [2] Valejo, Alan and Oliveira, Maria C. F. and Filho, Geraldo P. R. and Lopes, A. A., Multilevel approach for
combinatorial optimization in bipartite network, in Knowledge-based systems, p. 45--61, vol. 151, 2018,
doi: https://doi.org/10.1016/j.knosys.2018.03.021

.. [3] Valejo, Alan and Ferreira, V. and Oliveira, Maria C. F. and Lopes, A. A., Community detection in bipartite
network: a modified coarsening approach, in International symposium on information management and big data (SIMBig),
track on social network and media analysis and mining (SNMAN). Part of the Communications in Computer and Information
Science book series (CCIS, volume 795), p. 123--136, 2017, doi: https://doi.org/10.1007/978-3-319-90596-9_9

.. [4] Valejo, Alan and Filho, Geraldo P. R. and Oliveira, Maria C. F. and Ferreira, V. and Lopes, A. A., One-mode
projection-based multilevel approach for community detection in bipartite networks, in International symposium on
information management and big data (SIMBig), track on social network and media analysis and mining (SNMAN),
p. 101-108, 2017 , doi: http://ceur-ws.org/Vol-2029/paper8.pdf

.. [5] Valejo, Alan and Ferreira, V. and Fabbri, R. and Oliveira, Maria C. F. and Lopes, A. A., A critical survey of the
multilevel method in complex networks, in ACM Computing Surveys, accepted paper, 2019

Warning
-------
The original implementations (i.e. paper versions [1,2,3,4]) are deprecated.
There may be divergences between this version and the original algorithm.
If you looking for the original version used in the paper don't hesitate to contact the authors.

This software is a new version, more robust and fast.
It is a beta version and has some bugs and inconsistencies.
The final version of this tool will be released is coming soon.
For now, you can use this version without guarantee of the results.
"""

import sys
import numpy
import os
import inspect

import json

import models.args as args
import models.helper as helper
import models.helpermgraph as helpermgraph
from models.timing import Timing
from models.similarity import Similarity

import sharedmem
from multiprocessing import Process

__maintainer__ = 'Alan Valejo'
__email__ = 'alanvalejo@gmail.com'
__author__ = 'Alan Valejo'
__credits__ = ['Alan Valejo']
__homepage__ = 'https://www.alanvalejo.com.br'
__license__ = 'GNU.GPL.v3'
__docformat__ = 'markdown en'
__version__ = '0.1'
__date__ = '2019-08-08'

def main():
    """
    Main entry point for the application when run from the command line.
    """

    # Timing instanciation
    timing = Timing(['Snippet', 'Time [m]', 'Time [s]'])

    with timing.timeit_context_add('Pre-processing'):

        # Setup parse options command line
        current_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        parser = args.setup_parser(current_path + '/args/coarsening.json')
        options = parser.parse_args()
        args.update_json(options)
        args.check_output(options)

        # Log instanciation
        log = helper.initialize_logger(dir='log', output='log')

        if options.input and options.vertices is None:
            log.warning('Vertices are required when input is given.')
            sys.exit(1)

        # Create default values for optional parameters
        if options.reduction_factor is None:
            options.reduction_factor = [0.5] * len(options.vertices)
        if options.max_levels is None:
            options.max_levels = [3] * len(options.vertices)
        if options.matching is None:
            options.matching = ['rgmb'] * len(options.vertices)
        if options.similarity is None:
            options.similarity = ['weighted_common_neighbors'] * len(options.vertices)
        if options.itr is None:
            options.itr = [10] * len(options.vertices)
        if options.upper_bound is None:
            options.upper_bound = [2.0] * len(options.vertices)
        if options.global_min_vertices is None:
            options.global_min_vertices = [None] * len(options.vertices)
        if options.tolerance is None:
            options.tolerance = [0.01] * len(options.vertices)

        # Validation of list values
        if len(options.reduction_factor) == 1:
            options.reduction_factor = [options.reduction_factor[0]] * len(options.vertices)
        if len(options.max_levels) == 1:
            options.max_levels = [options.max_levels[0]] * len(options.vertices)
        if len(options.matching) == 1:
            options.matching = [options.matching[0]] * len(options.vertices)
        if len(options.similarity) == 1:
            options.similarity = [options.similarity[0]] * len(options.vertices)
        if len(options.itr) == 1:
            options.itr = [options.itr[0]] * len(options.vertices)
        if len(options.upper_bound) == 1:
            options.upper_bound = [options.upper_bound[0]] * len(options.vertices)
        if len(options.global_min_vertices) == 1:
            options.global_min_vertices = [options.global_min_vertices[0]] * len(options.vertices)
        if len(options.tolerance) == 1:
            options.tolerance = [options.tolerance[0]] * len(options.vertices)

        # Verification of the dimension of the parameters
        if len(options.vertices) != len(options.reduction_factor):
            log.warning('Sizes of input arguments -v and -r do not match.')
            sys.exit(1)
        if len(options.vertices) != len(options.max_levels):
            log.warning('Sizes of input arguments -v and -m do not match.')
            sys.exit(1)
        if len(options.vertices) != len(options.matching):
            log.warning('Sizes of input arguments -v and -c do not match.')
            sys.exit(1)
        if len(options.vertices) != len(options.similarity):
            log.warning('Sizes of input arguments -v and -s do not match.')
            sys.exit(1)
        if len(options.vertices) != len(options.itr):
            log.warning('Size of input arguments -v and -i do not match.')
            sys.exit(1)
        if len(options.vertices) != len(options.upper_bound):
            log.warning('Size of input arguments -v and -ub do not match.')
            sys.exit(1)
        if len(options.vertices) != len(options.global_min_vertices):
            log.warning('Size of input arguments -v and -gmv do not match.')
            sys.exit(1)
        if len(options.vertices) != len(options.tolerance):
            log.warning('Size of input arguments -v and -t do not match.')
            sys.exit(1)

        # Validation of matching method
        valid_matching = ['rgmb', 'gmb', 'mlpb', 'hem', 'lem', 'rm']
        for index, matching in enumerate(options.matching):
            matching = matching.lower()
            if matching not in valid_matching:
                log.warning('Matching ' + matching + ' method is unvalid.')
                sys.exit(1)
            options.matching[index] = matching

        # Validation of sedd priority
        valid_seed_priority = ['strength', 'degree', 'random']
        for index, seed_priority in enumerate(options.seed_priority):
            seed_priority = seed_priority.lower()
            if seed_priority not in valid_seed_priority:
                log.warning('Seed priotiry ' + seed_priority + ' is unvalid.')
                sys.exit(1)
            options.seed_priority[index] = seed_priority

        # Validation reverse
        for index, reverse in enumerate(options.reverse):
            if reverse.lower() in ('yes', 'true', 't', 'y', '1'):
                options.reverse[index] = True
            elif reverse.lower() in ('no', 'false', 'f', 'n', '0'):
                options.reverse[index] = False
            else:
                log.warning('Boolean value expected in -rv.')
                sys.exit(1)

        # isinstance(data[i][k], bool)

        # Validation of similarity measure
        valid_similarity = ['common_neighbors', 'weighted_common_neighbors',
        'salton', 'preferential_attachment', 'jaccard', 'weighted_jaccard',
        'adamic_adar', 'resource_allocation', 'sorensen', 'hub_promoted',
        'hub_depressed', 'leicht_holme_newman']
        for index, similarity in enumerate(options.similarity):
            similarity = similarity.lower()
            if similarity not in valid_similarity:
                log.warning('Similarity ' + similarity + ' misure is unvalid.')
                sys.exit(1)
            options.similarity[index] = similarity

        for layer in range(len(options.vertices)):
            if options.matching[layer] in ['rgmb', 'gmb', 'hem', 'lem', 'rm']:
                if options.global_min_vertices[layer] is not None:
                    options.global_min_vertices[layer] = None
                    text = 'Matching method ' + options.matching[layer]
                    text += ' (setted in layer '
                    text += str(layer) + ') does not accept -gmv parameter.'
                    log.warning(text)
                if options.reduction_factor[layer] > 0.5:
                    options.reduction_factor[layer] = 0.5
                    text = 'Matching method ' + options.matching[layer]
                    text += ' (setted in layer '
                    text += str(layer) + ') does not accept -rf > 0.5.'
                    log.warning(text)

    # Load bipartite graph
    with timing.timeit_context_add('Load'):
        graph = helpermgraph.load(options.input, options.vertices)
        graph['level'] = [0] * graph['layers']
        source_ecount = graph.ecount()

    # Coarsening
    with timing.timeit_context_add('Coarsening'):
        hierarchy_graphs = []
        hierarchy_levels = []
        running = True
        while running:

            running = False

            membership = sharedmem.full(graph.vcount(), range(graph.vcount()), dtype='int')
            levels = graph['level']
            contract = False

            processes = []
            for layer in range(len(graph['vertices'])):
                matching_layer = True
                if (options.global_min_vertices[layer] is None):
                    if levels[layer] >= options.max_levels[layer]:
                        matching_layer = False
                elif (graph['vertices'][layer] <= options.global_min_vertices[layer]):
                    matching_layer = False

                if matching_layer:
                    contract = True
                    running = True
                    levels[layer] += 1

                    graph['similarity'] = getattr(Similarity(graph, graph['adjlist']), options.similarity[layer])
                    start = sum(graph['vertices'][0:layer])
                    end = sum(graph['vertices'][0:layer + 1])
                    vertices = range(start, end)

                    param = dict(reduction_factor=options.reduction_factor[layer])

                    if options.matching[layer] in ['mlpb', 'gmb', 'rgmb']:
                        param['vertices'] = vertices
                        param['reverse'] = options.reverse[layer]
                    if options.matching[layer] in ['mlpb', 'rgmb']:
                        param['seed_priority'] = options.seed_priority[layer]
                    if options.matching[layer] in ['mlpb']:
                        param['upper_bound'] = options.upper_bound[layer]
                        param['n'] = options.vertices[layer]
                        param['global_min_vertices'] = options.global_min_vertices[layer]
                        param['tolerance'] = options.tolerance[layer]
                        param['itr'] = options.itr[layer]

                    if options.matching[layer] in ['hem', 'lem', 'rm']:
                        one_mode_graph = graph.weighted_one_mode_projection(vertices)
                        matching_method = getattr(one_mode_graph, options.matching[layer])
                    else:
                        matching_method = getattr(graph, options.matching[layer])

                    processes.append(Process(target=matching_method, args=[membership], kwargs=param))

            for p in processes:
                p.start()
            for p in processes:
                p.join()

            if contract:
                coarse = graph.contract(membership)
                coarse['level'] = levels

                if coarse.vcount() == graph.vcount():
                    break

                graph = coarse

                if options.save_hierarchy or not running:
                    hierarchy_graphs.append(graph)
                    hierarchy_levels.append(levels[:])

    if not hierarchy_graphs:
        hierarchy_graphs.append(graph)
        hierarchy_levels.append(levels[:])

    # Save
    with timing.timeit_context_add('Save'):

        output = options.output
        for index, obj in enumerate(zip(hierarchy_levels, hierarchy_graphs)):
            levels, graph = obj
            index += 1

            if options.save_conf or options.show_conf:
                d = {}
                d['source_input'] = options.input
                d['source_vertices'] = [options.vertices[0], options.vertices[1]]
                d['source_vcount'] = options.vertices[0] + options.vertices[1]
                d['source_ecount'] = source_ecount
                d['ecount'] = graph.ecount()
                d['vcount'] = graph.vcount()
                d['vertices'] = graph['vertices']
                d['reduction_factor'] = options.reduction_factor
                d['max_levels'] = options.max_levels
                d['achieved_levels'] = graph['level']
                d['similarity'] = options.similarity
                d['matching'] = options.matching
                d['level'] = levels
                d['upper_bound'] = options.upper_bound
                d['global_min_vertices'] = options.global_min_vertices
                d['itr'] = options.itr

            if options.save_conf:
                with open(output + '-' + str(index) + '.conf', 'w+') as f:
                    json.dump(d, f, indent=4)

            if options.show_conf:
                print(json.dumps(d, indent=4))

            if options.save_ncol:
                graph.write(output + '-' + str(index) + '.ncol', format='ncol')

            if options.save_source:
                with open(output + '-' + str(index) + '.source', 'w+') as f:
                    for v in graph.vs():
                        f.write(' '.join(map(str, v['source'])) + '\n')

            if options.save_membership:
                membership = [0] * (options.vertices[0] + options.vertices[1])
                for v in graph.vs():
                    for source in v['source']:
                        membership[source] = v.index
                numpy.savetxt(output + '-' + str(index) + '.membership', membership, fmt='%d')

            if options.save_predecessor:
                with open(output + '-' + str(index) + '.predecessor', 'w+') as f:
                    for v in graph.vs():
                        f.write(' '.join(map(str, v['predecessor'])) + '\n')

            if options.save_successor:
                numpy.savetxt(output + '-' + str(index) + '.successor', graph.vs['successor'], fmt='%d')

            if options.save_weight:
                numpy.savetxt(output + '-' + str(index) + '.weight', graph.vs['weight'], fmt='%d')

            if options.save_gml:
                del graph['adjlist']
                del graph['similarity']
                graph['layers'] = str(graph['layers'])
                graph['vertices'] = ','.join(map(str, graph['vertices']))
                graph['level'] = ','.join(map(str, graph['level']))
                graph.vs['name'] = map(str, range(0, graph.vcount()))
                graph.vs['type'] = map(str, graph.vs['type'])
                graph.vs['weight'] = map(str, graph.vs['weight'])
                graph.vs['successor'] = map(str, graph.vs['successor'])
                for v in graph.vs():
                    v['source'] = ','.join(map(str, v['source']))
                    v['predecessor'] = ','.join(map(str, v['predecessor']))
                graph.write(output + '-' + str(index) + '.gml', format='gml')

            if not options.save_hierarchy:
                break

    if options.show_timing:
        timing.print_tabular()
    if options.save_timing_csv:
        timing.save_csv(output + '-timing.csv')
    if options.save_timing_json:
        timing.save_json(output + '-timing.json')

if __name__ == "__main__":
    sys.exit(main())
