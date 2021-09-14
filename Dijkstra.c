// Dijkstra API implementation
// COMP2521 Assignment 2
//
// Overview:
// Applies Dijkstra's algorithm to determine the shortest path between
// each node in a graph and a given "source" node.
//
// Author:
// Lachlan Scott (z5207471@unsw.edu.au)
//
// Written:
// 12/07/21 - 06/08/21

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Dijkstra.h"
#include "Graph.h"
#include "PQ.h"

#define INFINITY INT_MAX

static PredNode *predNodeNew(int v);

// Applies Dijkstra's algorithm to determine the shortest path between
// each node in a graph and a given "source" node.
ShortestPaths dijkstra(Graph g, Vertex src) {
	ShortestPaths sps = {0};

	sps.numNodes = GraphNumVertices(g);
	sps.src = src;

	PQ vSet = PQNew();

	sps.dist = malloc(sps.numNodes * sizeof(int));

	if (sps.dist == NULL) {
		fprintf(stderr, "Couldn't allocate dist array!\n");
        exit(EXIT_FAILURE);
	}

	// Initialise src distance to 0, all others to infinity
	// Insert all vertices into priority queue
	for (int i = 0; i < sps.numNodes; i++) {
		sps.dist[i] = (i == src) ? 0 : INFINITY;
		PQInsert(vSet, i, sps.dist[i]);
	}

	sps.pred = (PredNode **)malloc(sps.numNodes * sizeof(PredNode *));

	if (sps.pred == NULL) {
		fprintf(stderr, "Couldn't allocate pred array!\n");
        exit(EXIT_FAILURE);
	}

	// Initialise all predecessors to NULL
	for (int i = 0; i < sps.numNodes; i++) {
		sps.pred[i] = NULL;
	}

	// Algorithm completes when all nodes have been visited
	while (!PQIsEmpty(vSet)) {
		Vertex curr = PQDequeue(vSet);

		// If shortest distance to current node is infinity, the node cannot be reached
		// from the source node, skip it
		if (sps.dist[curr] == INFINITY) {
			continue;
		}

		AdjList currChild = GraphOutIncident(g, curr);

		// Iterate through the outlinks of the current node
		// Perform edge relaxation based on new edge information
		while (currChild != NULL) {
			int childV = currChild->v;
			int distToChild = sps.dist[curr] + currChild->weight;

			// Perform edge relaxation
			if (distToChild < sps.dist[childV]) {
				// Case where shorter path exists
				sps.dist[childV] = distToChild;
				PredNode *pn = predNodeNew(curr);

				// Free previous predecessor list
				PredNode *currPred = sps.pred[childV];
				while (currPred != NULL) {
					PredNode *prevPred = currPred;
					currPred = currPred->next;
					free(prevPred);
				}

				// Set node's new predecessor to current node
				// Update its position in the priority queue
				sps.pred[childV] = pn;
				PQUpdate(vSet, childV, sps.dist[childV]);
			} else if (distToChild == sps.dist[childV]) {
				// Case where equally short path exists
				PredNode *currPred = sps.pred[childV];

				// Find end of predecessor linked list
				while (currPred->next != NULL) {
					currPred = currPred->next;
				}

				PredNode *pn = predNodeNew(curr);
				currPred->next = pn;
			}
			currChild = currChild->next;
		}
	}
	PQFree(vSet);
	return sps;
}

// Allocate a new PredNode struct
// Produces a runtime error if the PredNode fails to be allocated
static PredNode *predNodeNew(int v) {
	PredNode *pn = (PredNode *)malloc(sizeof(PredNode));

	if (pn == NULL) {
		fprintf(stderr, "Couldn't allocate new PredNode!\n");
        exit(EXIT_FAILURE);
	}

	pn->v = v;
	pn->next = NULL;
	return pn;
}

void showShortestPaths(ShortestPaths sps) {

}

// Frees a ShortestPaths struct
void freeShortestPaths(ShortestPaths sps) {
	// Iterate through the predecessor array and free all PredNodes
	for (int i = 0; i < sps.numNodes; i++) {
		PredNode *curr = sps.pred[i];

		while (curr != NULL) {
			PredNode *prev = curr;
			curr = curr->next;
			free(prev);
		}
	}

	free(sps.pred);
	free(sps.dist);
}
