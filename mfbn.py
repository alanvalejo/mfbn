#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MFBN: Multilevel framework for bipartite networks

Copyright (C) 2017 Alan Valejo <alanvalejo@gmail.com> All rights reserved

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

from models.mgraph import MGraph
from models.coarsening import Coarsening
import models.args as args

from models.timing import Timing

__maintainer__ = 'Alan Valejo'
__email__ = 'alanvalejo@gmail.com'
__author__ = 'Alan Valejo'
__credits__ = ['Alan Valejo']
__homepage__ = 'https://www.alanvalejo.com.br'
__license__ = 'GNU.GPL.v3'
__docformat__ = 'markdown en'
__version__ = '0.1.0'
__date__ = '2020-04-25'


def main():
    """
    Main entry point for the application when run from the command line.
    """

    # Timing instance
    timing = Timing(['Snippet', 'Time [m]', 'Time [s]'])

    with timing.timeit_context_add('Pre-processing'):

        # Setup parse options command line
        current_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        parser = args.setup_parser(current_path + '/args/mfbn.json')
        options = parser.parse_args()
        args.update_json(options)
        args.check_output(options)

        if options.input and options.vertices is None:
            print('Vertices are required when input is given.')
            sys.exit(1)

    # Load bipartite graph
    with timing.timeit_context_add('Load graph'):

        source_graph = MGraph()
        source_graph.load(options.input, options.vertices)

    # Coarsening
    with timing.timeit_context_add('Coarsening'):

        kwargs = dict(
            reduction_factor=options.reduction_factor, max_levels=options.max_levels,
            matching=options.matching, similarity=options.similarity, itr=options.itr,
            upper_bound=options.upper_bound, gmv=options.gmv,
            tolerance=options.tolerance, reverse=options.reverse, seed_priority=options.seed_priority,
            threads=options.threads
        )

        coarsening = Coarsening(source_graph, **kwargs)
        coarsening.run()

    # Save
    with timing.timeit_context_add('Save'):

        output = options.output
        for index, obj in enumerate(zip(coarsening.hierarchy_levels, coarsening.hierarchy_graphs)):
            level, coarsened_graph = obj
            index += 1

            if options.save_conf or options.show_conf:
                d = {
                    'source_input': options.input
                    , 'source_vertices': source_graph['vertices']
                    , 'source_vcount': source_graph.vcount()
                    , 'source_ecount': source_graph.ecount()
                    , 'coarsened_ecount': coarsened_graph.ecount()
                    , 'coarsened_vcount': coarsened_graph.vcount()
                    , 'coarsened_vertices': coarsened_graph['vertices']
                    , 'achieved_levels': coarsened_graph['level']
                    , 'reduction_factor': options.reduction_factor
                    , 'max_levels': options.max_levels
                    , 'similarity': options.similarity
                    , 'matching': options.matching
                    , 'upper_bound': options.upper_bound
                    , 'gmv': options.gmv
                    , 'itr': options.itr
                    , 'level': level
                }

            if options.save_conf:
                with open(output + '-' + str(index) + '-info.json', 'w+') as f:
                    json.dump(d, f, indent=4)

            if options.show_conf:
                print(json.dumps(d, indent=4))

            if options.save_ncol:
                coarsened_graph.write(output + '-' + str(index) + '.ncol', format='ncol')

            if options.save_source:
                with open(output + '-' + str(index) + '.source', 'w+') as f:
                    for v in coarsened_graph.vs():
                        f.write(' '.join(map(str, v['source'])) + '\n')

            if options.save_membership:
                membership = [0] * (source_graph['vertices'][0] + source_graph['vertices'][1])
                for v in coarsened_graph.vs():
                    for source in v['source']:
                        membership[source] = v.index
                numpy.savetxt(output + '-' + str(index) + '.membership', membership, fmt='%d')

            if options.save_predecessor:
                with open(output + '-' + str(index) + '.predecessor', 'w+') as f:
                    for v in coarsened_graph.vs():
                        f.write(' '.join(map(str, v['predecessor'])) + '\n')

            if options.save_successor:
                numpy.savetxt(output + '-' + str(index) + '.successor', coarsened_graph.vs['successor'], fmt='%d')

            if options.save_weight:
                numpy.savetxt(output + '-' + str(index) + '.weight', coarsened_graph.vs['weight'], fmt='%d')

            if options.save_gml:
                del coarsened_graph['adjlist']
                del coarsened_graph['similarity']
                coarsened_graph['layers'] = str(coarsened_graph['layers'])
                coarsened_graph['vertices'] = ','.join(map(str, coarsened_graph['vertices']))
                coarsened_graph['level'] = ','.join(map(str, coarsened_graph['level']))
                coarsened_graph.vs['name'] = map(str, range(0, coarsened_graph.vcount()))
                coarsened_graph.vs['type'] = map(str, coarsened_graph.vs['type'])
                coarsened_graph.vs['weight'] = map(str, coarsened_graph.vs['weight'])
                coarsened_graph.vs['successor'] = map(str, coarsened_graph.vs['successor'])
                for v in coarsened_graph.vs():
                    v['source'] = ','.join(map(str, v['source']))
                    v['predecessor'] = ','.join(map(str, v['predecessor']))
                coarsened_graph.write(output + '-' + str(index) + '.gml', format='gml')

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
