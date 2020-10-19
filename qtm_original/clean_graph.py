# Author: Gitanjali Bhattacharjee
# The functions in this file take the original graph of the road network and remove edges to which traffic can never be
# assigned.

import mahmodel_road_only as mahmodel

import bd, pickle
import util
import networkx as nx

def correct_graph(G): # Correct properties of the original graph from Miller.
	count = 0 # if capacity of edge is 0, its t_a should be 0
	count0 = 0 # if capacity of edge is 0, its t_0 should be 0
	count1 = 0 # if length of edge (distance_0) is 0, it shouldn't be able to be assigned traffic
	for edge in G.edges():
		if G[edge[0]][edge[1]]['capacity'] == 0:
			if G[edge[0]][edge[1]]['t_a'] != float('inf'):
				G[edge[0]][edge[1]]['t_a'] = float('inf')
				count += 1
			if G[edge[0]][edge[1]]['t_0'] != float('inf'):
				G[edge[0]][edge[1]]['t_0'] = float('inf')
				count0 += 1
		if G[edge[0]][edge[1]]['distance_0'] == 0:
			G[edge[0]][edge[1]]['capacity'] = 0
			G[edge[0]][edge[1]]['t_0'] = float('inf')
			G[edge[0]][edge[1]]['t_a'] = float('inf')
			G[edge[0]][edge[1]]['distance'] = 0

			count1 += 1
			#print G[edge[0]][edge[1]]['distance'], G[edge[0]][edge[1]]['capacity'], G[edge[0]][edge[1]]['t_a']

	#print 'corrected ', count, ' and ', count1, ' edges in graph'

	return G

def prune_existing_graph(G): # remove all edges that cannot be assigned trips, either because they have 0 capacity or 0 distance.

	G_copy = nx.DiGraph(G)  # copy so we can remove edges

	edges_to_remove = [(u, v) for (u, v) in G_copy.edges() if
					   G_copy[u][v]['capacity'] <= 0 or G_copy[u][v]['distance_0'] <= 0]

	G_copy.remove_edges_from(edges_to_remove)

	demand = bd.build_demand('input/BATS2000_34SuperD_TripTableData.csv',
							 'input/superdistricts_centroids_dummies.csv')

	supernodes = demand.keys()

	count_supers = 0
	for n in supernodes:
		if G_copy.has_node(n) is False:
			pass
		else:
			count_supers += 1

	assert count_supers == len(supernodes), 'ERROR: Deleted a supernode in prune_graph().'
	# print 'Number of supernodes = ', count_supers, len(supernodes)

	return G_copy

def save_corrected_graph():

	G = nx.read_gpickle("input/graphMTC_CentroidsLength3int.gpickle")
	G = mahmodel.add_superdistrict_centroids(G)
	G = correct_graph(G)
	G = prune_existing_graph(G)
	G = nx.freeze(G)

	demand = bd.build_demand('input/BATS2000_34SuperD_TripTableData.csv',
							 'input/superdistricts_centroids_dummies.csv')

	supernodes = demand.keys()

	for n, nbrsdict in G.adjacency_iter():
		if n in supernodes:
			print n, len(nbrsdict.keys())

	if nx.is_frozen(G):
		G_copy = nx.DiGraph(G)
		nx.write_gpickle(G_copy, 'input/graphMTC_GB.gpickle')
	else:
		nx.write_gpickle(G, 'input/graphMTC_GB.gpickle')

def load_corrected_graph():

	G = nx.read_gpickle('input/graphMTC_GB.gpickle')
	G = nx.freeze(G)
	# print len(G.edges())

	return G

def main():

	save_corrected_graph()
