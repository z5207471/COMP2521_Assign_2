// Centrality Measures API implementation
// COMP2521 Assignment 2
//
// Overview:
// Provides functions which calculate centrality measures of a graph via three
// distinct methods. These are: closeness centrality, betweenness centrality and
// normalised betweenness centrality.
//
// Author:
// Lachlan Scott (z5207471@unsw.edu.au)
//
// Written:
// 12/07/21 - 06/08/21

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

#define INFINITY INT_MAX

// Stores the calculated values of sigma-st(v) and sigma-st
// used for betweenness centrality calculation
typedef struct BCSigmas {
	int nSP; // Number of shortest paths from s to t
	int nSPThrough; // Number of shortest paths from s to t through v
} BCSigmas;

static NodeValues newNodeValues(int nV);
static BCSigmas traverseSPS(ShortestPaths sps, int dest, int target);
static void recurTraverseSPS(ShortestPaths sps, int vertex, int target, int tFound, BCSigmas *sigmas);
static double wassermanFaust(int N, int n, int distSum);
static int numReachable(ShortestPaths sps, int src);
static double normaliseBC(double value, int n);

// Obtains the closeness centrality values for each node in a graph
NodeValues closenessCentrality(Graph g) {
	NodeValues nvs = newNodeValues(GraphNumVertices(g));

	// Get closenness centrality value for each node
	for (int i = 0; i < nvs.numNodes; i++) {
		ShortestPaths sps = dijkstra(g, i);
		int distSum = 0;

		// Get sum of all shortest path dists and number of reachable nodes
		for (int j = 0; j < nvs.numNodes; j++) {
			// If a shortest path exists to the node, add its length
			// to the distance sum and increment number of reachable nodes
			if (sps.dist[j] != INFINITY) {
				distSum += sps.dist[j];
			}
		}

		int nReachable = numReachable(sps, i);

		freeShortestPaths(sps);

		if (nReachable == 1) {
			// Closeness val for isolated node should be 0
			nvs.values[i] = 0.0;
		} else {
			// Closeness val for non-isolated nodes calculated using Wasserman and Faust
			nvs.values[i] = wassermanFaust(nvs.numNodes, nReachable, distSum);
		}
	}

	return nvs;
}

// Obtains the betweenness centrality values for each node in a graph
NodeValues betweennessCentrality(Graph g) {
	NodeValues nvs = newNodeValues(GraphNumVertices(g));

	// Get betweenness centrality value for each node i
	for (int i = 0; i < nvs.numNodes; i++) {
		double centralitySum = 0.0;

		// j is the source node in each iteration
		for (int j = 0; j < nvs.numNodes; j++) {
			ShortestPaths sps = dijkstra(g, j);

			// k is the destination node in each iteration
			for (int k = 0; k < nvs.numNodes; k++) {

				// Must obey condition i != j != k
				if (j != k && i != j && i != k) {
					BCSigmas sigmas = { 0 };

					if (sps.dist[k] != INFINITY) {
						// If at least one shortest path from j to k exists, 
						// count shortest paths and paths through i
						sigmas = traverseSPS(sps, k, i);
					}

					// Add ratio of (shortest paths through i to shortest paths) to
					// total centrality sum
					if (sigmas.nSPThrough > 0 && sigmas.nSP > 0) {
						centralitySum += (double)sigmas.nSPThrough / (double)sigmas.nSP;
					}
				}
			}
			freeShortestPaths(sps);
		}
		nvs.values[i] = centralitySum;
	}

	return nvs;
}

// Obtains the normalised betweenness centrality values for each node
// in a graph
NodeValues betweennessCentralityNormalised(Graph g) {
	NodeValues nvs = betweennessCentrality(g);

	// Normalise each calculated BC value
	for (int i = 0; i < nvs.numNodes; i++) {
		nvs.values[i] = normaliseBC(nvs.values[i], GraphNumVertices(g));
	}

	return nvs;
}

// Creates a new NodeValues struct
// Produces a runtime error if the values array fails to be allocated
static NodeValues newNodeValues(int nV) {
	NodeValues nvs = {0};
	nvs.numNodes = nV;
	nvs.values = (double *)malloc(nvs.numNodes * sizeof(double));

	if (nvs.values == NULL) {
		fprintf(stderr, "Couldn't allocate NodeValues values array!\n");
        exit(EXIT_FAILURE);
	}

	return nvs;
}

// Traverses through a ShortestPaths pred array, counting the
// number of separate shortest paths from the destination vertex to the
// source vertex as well as the number that pass through the target vertex
static BCSigmas traverseSPS(ShortestPaths sps, int dest, int target) {
	BCSigmas sigmas = { 1, 0 };

	int targetFound = false;
	recurTraverseSPS(sps, dest, target, targetFound, &sigmas);

	return sigmas;
}

// Recursively traverses through a ShortestPaths pred array, counting the
// number of separate shortest paths from the destination vertex to the
// source vertex as well as the number that pass through the target vertex
static void recurTraverseSPS(ShortestPaths sps, int vertex, int target, int tFound, BCSigmas *sigmas) {
	// Check whether the target vertex lies along the current path
	int throughTarget = tFound;
	if (throughTarget == false && vertex == target) {
		sigmas->nSPThrough++; // Target found, at least one path through target exists
		throughTarget = true;
	}

	// Iterate through, count and traverse each predecessor of the current vertex
	int predCount = 0;
	PredNode *curr = sps.pred[vertex];
	while (curr != NULL) {
		predCount++;
		recurTraverseSPS(sps, curr->v, target, throughTarget, sigmas);
		curr = curr->next;
	}

	// Check whether there are multiple paths branching from the current vertex
	if (predCount > 1) {
		// If there are multiple branches, each one after the first counts as a
		// new shortest path
		sigmas->nSP += predCount - 1;

		// Once target has been found along a path, every path divergence
		// afterward is counted as a separate path containing the target
		if (throughTarget == true) {
			sigmas->nSPThrough += predCount - 1;
		}
	}
}

void showNodeValues(NodeValues nvs) {

}

// Frees a given NodeValues struct
void freeNodeValues(NodeValues nvs) {
	free(nvs.values);
}

// Apply the Wasserman and Faust formula to calculate
// closeness centrality
// Assumes the corresponding graph contains at least 3 nodes
static double wassermanFaust(int N, int n, int distSum) {
	double dN = (double)N;
	double dn = (double)n;
	double dDistSum = (double)distSum;

	return (((dn - 1) / (dN - 1)) * ((dn - 1) / (dDistSum)));
}

// Returns the number of nodes in a graph reachable from a given source node.
// Assumes that the input ShortestPaths struct has been generated with respect
// to the given source node.
// Assumes the corresponding graph contains at least 3 nodes
static int numReachable(ShortestPaths sps, int src) {
	int output = 0;
	for (int i = 0; i < sps.numNodes; i++) {
		if (sps.dist[i] != INFINITY) {
			// If a path distance is not infinite, the path exists
			output++;
		}
	}
	return output;
}

// Normalises a betweenness centrality value
// Assumes the corresponding graph contains at least 3 nodes
static double normaliseBC(double value, int n) {
	return ((1 / (((double)n - 1) * ((double)n - 2))) * value);
}
