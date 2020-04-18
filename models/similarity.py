#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Similarity measures

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

import math

__maintainer__ = 'Alan Valejo'
__email__ = 'alanvalejo@gmail.com'
__author__ = 'Alan Valejo'
__credits__ = ['Alan Valejo']
__homepage__ = 'https://www.alanvalejo.com.br'
__license__ = 'GNU.GPL.v3'
__docformat__ = 'markdown en'
__version__ = '0.1'
__date__ = '2019-08-08'

class Similarity(object):

	graph, adjlist = (None,) * 2

	def __init__(self, graph, adjlist):
		self.graph = graph
		self.adjlist = adjlist

	def weight(self, i, j):
		""" Calculates pairwise weight edge on a geiven graph. """

		return self.graph[i, j]

	def preferential_attachment(self, i, j):
		""" Calculates pairwise preferential attachment similarities on a given unweighted graph. """

		return float(self.graph.degree(i)) * float(self.graph.degree(j))

	def common_neighbors(self, i, j):
		""" Calculates pairwise common neighbors similarity on a given unweighted graph. """

		return len(self.adjlist[i].intersection(self.adjlist[j]))

	def get_common_neighbors(self, i, j):
		""" Calculates pairwise common neighbors similarity on a given unweighted graph. """

		return self.adjlist[i].intersection(self.adjlist[j])

	def weighted_common_neighbors(self, i, j):
		"""
		Calculates pairwise common neighbors similarity uses the edge-weight information
		on a given unweighted graph.
		"""

		_sum = 0.0
		cn = self.adjlist[i].intersection(self.adjlist[j])
		for z in cn:
			_sum += (self.graph[i, z] + self.graph[j, z]) / 2
		return _sum

	def jaccard(self, i, j):
		""" Calculates pairwise jaccard similarity on a given unweighted graph. """

		isect = len(self.adjlist[i].intersection(self.adjlist[j]))
		union = (len(self.adjlist[i]) + len(self.adjlist[j]) - isect)
		return 0 if union == 0 else isect / float(union)

	def weighted_jaccard(self, i, j):
		""" Calculates pairwise jaccard similarity on a given unweighted graph. """

		_sum_isect = 0.0
		isect = self.adjlist[i].intersection(self.adjlist[j])
		for z in isect:
			_sum_isect += (self.graph[i, z] + self.graph[j, z]) / 2

		_sum_union = 0.0
		union = self.adjlist[i].union(self.adjlist[j]) - isect
		for z in union:
			_sum_union += (self.graph[i, z] + self.graph[j, z]) / 2

		return 0 if _sum_union == 0.0 else _sum_isect / _sum_union

	def salton(self, i, j):
		""" Calculates pairwise solton similarity on a given unweighted graph. """

		product = float(self.graph.degree(i)) * float(self.graph.degree(j))
		if product == 0.0:
			return 0.0

		isect = len(self.adjlist[i].intersection(self.adjlist[j]))
		return isect / math.sqrt(product)

	def adamic_adar(self, i, j):
		""" Calculates pairwise adamic adar similarity on a given unweighted graph. """

		score = 0.0
		for isect in self.adjlist[i].intersection(self.adjlist[j]):
			degree = self.graph.degree(isect)
			if degree != 0:
				score += 1 / math.log(degree)
		return score

	def resource_allocation(self, i, j):
		""" Calculates pairwise resource allocation similarity on a given unweighted graph. """

		score = 0.0
		for isect in self.adjlist[i].intersection(self.adjlist[j]):
			degree = self.graph.degree(isect)
			if degree != 0:
				score += 1 / degree
		return score

	def sorensen(self, i, j):
		""" Calculates pairwise sorensen similarity on a given unweighted graph. """

		_sum = float(self.graph.degree(i)) * float(self.graph.degree(j))
		if _sum == 0.0:
			return 0.0

		isect = 2 * len(self.adjlist[i].intersection(self.adjlist[j]))
		return isect / _sum

	def hub_promoted(self, i, j):
		""" Calculates pairwise hub promoted similarity on a given unweighted graph. """

		minimum = min(float(self.graph.degree(i)), float(self.graph.degree(j)))
		if minimum == 0.0:
			return 0.0

		isect = len(self.adjlist[i].intersection(self.adjlist[j]))
		return isect / minimum

	def hub_depressed(self, i, j):
		""" Calculates pairwise hub depressed similarity on a given unweighted graph. """

		maximum = max(float(self.graph.degree(i)), float(self.graph.degree(j)))
		if maximum == 0.0:
			return 0.0

		isect = len(self.adjlist[i].intersection(self.adjlist[j]))
		return isect / maximum

	def leicht_holme_newman(self, i, j):
		""" Calculates pairwise leicht holmeNewman similarity on a given unweighted graph. """

		product = float(self.graph.degree(i)) * float(self.graph.degree(j))
		if product == 0.0:
			return 0.0

		isect = len(self.adjlist[i].intersection(self.adjlist[j]))
		return isect / product

	def within_common_neighbors(self, i, j):
		"""
		Calculates pairwise inter-cluster common neighbors (WCN) similarity on a
		given unweighted graph. WCN considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		isect = self.adjlist[i].intersection(self.adjlist[j])
		within_isect = 0.0
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				within_isect += 1.0
		return within_isect

	def within_jaccard(self, i, j):
		"""
		Calculates pairwise inter-cluster jaccard (WJac) similarity on a
		given unweighted graph. WJac considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		isect = self.adjlist[i].intersection(self.adjlist[j])
		within_isect = 0.0
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				within_isect += 1.0
		union = (len(self.adjlist[i]) + len(self.adjlist[j]) - within_isect)
		return 0 if union == 0 else within_isect / float(union)

	def within_salton(self, i, j):
		"""
		Calculates pairwise inter-cluster salton (WSal) similarity on a
		given unweighted graph. WSal considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		product = float(self.graph.degree(i)) * float(self.graph.degree(j))

		if product == 0.0:
			return 0.0

		isect = self.adjlist[i].intersection(self.adjlist[j])
		within_isect = 0.0
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				within_isect += 1.0
		return within_isect / math.sqrt(product)

	def within_adamic_adar(self, i, j):
		"""
		Calculates pairwise inter-cluster adamic adar (WAA) similarity on a
		given unweighted graph. WAA considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		score = 0.0
		for isect in self.adjlist[i].intersection(self.adjlist[j]):
			if self.vs[isect]['membership'] == self.vs[i]['membership']:
				degree = self.graph.degree(isect)
				if degree != 0:
					score += 1 / math.log(degree)
		return score

	def within_resource_allocation(self, i, j):
		"""
		Calculates pairwise inter-cluster resource allocation (WRL) similarity on a
		given unweighted graph. WRL considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		score = 0.0
		for isect in self.adjlist[i].intersection(self.adjlist[j]):
			if self.vs[isect]['membership'] == self.vs[i]['membership']:
				degree = self.graph.degree(isect)
				if degree != 0:
					score += 1 / degree
		return score

	def within_sorensen(self, i, j):
		"""
		Calculates pairwise inter-cluster sorensen (WSor) similarity on a
		given unweighted graph. WSor considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		_sum = float(self.graph.degree(i)) * float(self.graph.degree(j))
		if _sum == 0.0:
			return 0.0

		isect = self.adjlist[i].intersection(self.adjlist[j])
		within_isect = 0.0
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				within_isect += 2.0
		return within_isect / _sum

	def within_hub_promoted(self, i, j):
		"""
		Calculates pairwise inter-cluster hub promoted (WHP) similarity on a
		given unweighted graph. WHP considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		minimum = min(float(self.graph.degree(i)), float(self.graph.degree(j)))
		if minimum == 0.0:
			return 0.0

		isect = self.adjlist[i].intersection(self.adjlist[j])
		within_isect = 0.0
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				within_isect += 1.0
		return within_isect / minimum

	def within_hub_depressed(self, i, j):
		"""
		Calculates pairwise inter-cluster hub depressed (WHD) similarity on a
		given unweighted graph. WHD considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		maximum = max(float(self.graph.degree(i)), float(self.graph.degree(j)))
		if maximum == 0.0:
			return 0.0

		isect = self.adjlist[i].intersection(self.adjlist[j])
		within_isect = 0.0
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				within_isect += 1.0
		return within_isect / maximum

	def within_leicht_holme_newman(self, i, j):
		"""
		Calculates pairwise inter-cluster leicht holme newman (WLHN) similarity on a
		given unweighted graph. WLHN considering solely the subset of within-cluster
		common neighbors instead of the set of all common neighbors
		"""

		product = float(self.graph.degree(i)) * float(self.graph.degree(j))
		if product == 0.0:
			return 0.0

		isect = self.adjlist[i].intersection(self.adjlist[j])
		within_isect = 0.0
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				within_isect += 1.0
		return within_isect / product

	def wic(self, i, j):
		"""
		Calculates similarity between a pair of vertices using information from
		intra-cluster or within-cluster (W) and inter-cluster (IC) common
		neighbors of these vertices.
		"""

		isect = self.adjlist[i].intersection(self.adjlist[j])
		nWcn = 0.0 # Intra cluster or intra community
		nIcn = 0.0 # Inter clusters or inter comunities
		for vertex in isect:
			if self.vs[vertex]['membership'] == self.vs[i]['membership']:
				nWcn += 1.0
			else:
				nIcn += 1.0
		return nWcn if (nIcn == 0.0) else nWcn / nIcn
