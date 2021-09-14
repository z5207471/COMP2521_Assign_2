// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2
//
// Overview:
// Performs hierarchical agglomerative clustering using the Lance-Williams 
// algorithm. Applies either single linkage or complete linkage method based
// on the "method" input parameter. Uses a pre-defined graph ADT.
//
// Author:
// Lachlan Scott (z5207471@unsw.edu.au)
//
// Written:
// 12/07/21 - 06/08/21

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX

// A set of two integers used to specify indices in an array of
// a pair of target elements (e.g. graph vertics, clusters)
typedef struct targets {
    int a;
    int b;
} targets;

// A set of integers representing the graph vertices contained within
// a Dendrogram struct
typedef struct vertexSet {
    int *vertices;
    int nV;
} vertexSet;

static double dMax(double a, double b);
static double dMin(double a, double b);
static double clusterDist(Graph g, vertexSet vSet1, vertexSet vSet2, int method);
static double vertexDist(Graph g, targets t);
static void mergeDends(Dendrogram *dends, int size, targets tDends);
static Dendrogram *allocDends(int size);
static void freeVertexSet(vertexSet *v);
static vertexSet getVertices(Dendrogram d, int maxVertices);
static void recurGetVertices(Dendrogram d, int *nV, int *vertices);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
    // Handle invalid method value input
    if (method != SINGLE_LINKAGE && method != COMPLETE_LINKAGE) {
        fprintf(stderr, "Invalid method int value (%d) given to LanceWilliamsHAC function.\n", method);
        exit(EXIT_FAILURE);
    }

    int nV = GraphNumVertices(g);

    // Initial number of dendrograms is number of graph vertices
    Dendrogram *dends = allocDends(nV);
    for (int i = 0; i < nV; i++) {
        dends[i]->vertex = i;
    }

    // Iterate over graph nV times to produce HAC dendrogram of given graph
    for (int numDends = nV; numDends > 1; numDends--) {
        targets dendsToMerge = { 0, 1 };
        double minClusterDist = INFINITY;

        // Iterate through each dendrogram, determine dist. to each other dendrogram
        for (int i = 0; i < numDends; i++) {
            vertexSet vSet1 = getVertices(dends[i], nV);

            for (int j = i + 1; j < numDends; j++) {
                vertexSet vSet2 = getVertices(dends[j], nV);

                double cDist = clusterDist(g, vSet1, vSet2, method);

                freeVertexSet(&vSet2);

                // Update dendrograms to be merged and minimum cluster distance
                if (cDist < minClusterDist) {
                    dendsToMerge.a = i;
                    dendsToMerge.b = j;
                    minClusterDist = cDist;
                }
            }
            freeVertexSet(&vSet1);
        }

        mergeDends(dends, numDends, dendsToMerge);
    }

    Dendrogram output = dends[0];
    free(dends);

	return output;
}

/**
 * Frees all memory associated with the given Dendrogram structure.
 */
void freeDendrogram(Dendrogram d) {
    if (d == NULL) {
        return;
    }

    freeDendrogram(d->left);
    freeDendrogram(d->right);
    free(d);
}

// Returns the maximum of two double values
static double dMax(double a, double b) {
    return (a > b) ? a : b;
}

// Returns the minimum of two double values
static double dMin(double a, double b) {
    return (a < b) ? a : b;
}

// Determines the distance between two clusters using either the single
// or complete linkage method
// Returns INFINITY if two clusters are not connected
static double clusterDist(Graph g, vertexSet vSet1, vertexSet vSet2, int method) {
    double cDist = INFINITY;

    // Iterate through vertices of each cluster to find cluster distance
    for (int i = 0; i < vSet1.nV; i++) {
        for (int j = 0; j < vSet2.nV; j++) {
            // Set current vertex pair to be measured
            targets tars = { vSet1.vertices[i], vSet2.vertices[j] };

            double vDist = vertexDist(g, tars);

            if (method == SINGLE_LINKAGE) {
                // Single linkage - find minimal distance between two clusters
                cDist = dMin(cDist, vDist);
            } else if (method == COMPLETE_LINKAGE) {
                // Complete linkage - find maximum distance between two clusters
                if (vDist != INFINITY) { // Only update dist. if vertices are connected
                    if (cDist == INFINITY) {
                        // Cluster dist has not been updated before
                        // Set to any value < infinity
                        cDist = vDist;
                    } else {
                        // Find maximum finite distance value
                        cDist = dMax(cDist, vDist);
                    }
                }
            }
        }
    }

    return cDist;
}

// Determines the distance between two vertices a and b in the graph
// Vertices are specified by the targets struct parameter t
// Returns INFINITY if the vertices are not connected
static double vertexDist(Graph g, targets t) {
    double distAB = INFINITY;

    AdjList aOutlinks = GraphOutIncident(g, t.a);

    // Determine whether an edge from a to b exists
    // If so, get its distance value
    while (aOutlinks != NULL) {
        if (aOutlinks->v == t.b) {
            distAB = (double)aOutlinks->weight;
            break;
        }
        aOutlinks = aOutlinks->next;
    }

    double distBA = INFINITY;

    AdjList bOutlinks = GraphOutIncident(g, t.b);

    // Determine whether an edge from b to a exists
    // If so, get its distance value
    while (bOutlinks != NULL) {
        if (bOutlinks->v == t.a) {
            distBA = (double)bOutlinks->weight;
            break;
        }
        bOutlinks = bOutlinks->next;
    }

    // Handle cases where 0 or 1 edges exist between a and b
    if (distAB == INFINITY && distBA == INFINITY) {
        return INFINITY;
    } else if (distAB == INFINITY) {
        return (1 / distBA);
    } else if (distBA == INFINITY) {
        return (1 / distAB);
    }

    // 2 edges exist, calculate distance using largest edge weight between a and b
    return (1 / dMax(distAB, distBA));
}

// Merges two Dendrograms, indicated by tDends, in an array of Dendrograms
// Shifts all values down in the array after merging to remove "gaps"
// Produces a runtime error if new Dendrogram fails to be allocated
static void mergeDends(Dendrogram *dends, int size, targets tDends) {
    assert(tDends.a > -1 && tDends.b > -1);
    
    // Combine target Dendrograms into a single Dendrogram, add it to dends
    Dendrogram combDend = (Dendrogram)malloc(sizeof(DNode));

    if (combDend == NULL) {
        fprintf(stderr, "Couldn't allocate Dendrogram!\n");
        exit(EXIT_FAILURE);
    }

    combDend->vertex = -1;
    combDend->left = dends[tDends.a];
    combDend->right = dends[tDends.b];
    dends[tDends.a] = combDend;
    dends[tDends.b] = NULL; // Unused Dendrogram is set to NULL

    // Shift all dendrograms after the NULL element down
    int foundEmpty = false;
    for (int i = 0; i < size; i++) {
        if (foundEmpty == true) {
            // Shift element down
            dends[i-1] = dends[i];
        } else if (dends[i] == NULL) {
            // Location of Dendrogram b found, shift all following Dendrograms down
            foundEmpty = true;
        }
    }
}

// Allocates the "dends" Dendrogram array used in the Lance-Williams formula
// Produces a runtime error if "dends"/a Dendrogram fails to be allocated
static Dendrogram *allocDends(int size) {
    Dendrogram *dends = (Dendrogram *)malloc(size * sizeof(Dendrogram));

    if (dends == NULL) {
        fprintf(stderr, "Couldn't allocate dends array!\n");
        exit(EXIT_FAILURE);
    }

    // Allocate all Dendrograms in the array
    for (int i = 0; i < size; i++) {
        dends[i] = (Dendrogram)malloc(sizeof(DNode));

        if (dends[i] == NULL) {
            fprintf(stderr, "Couldn't allocate Dendrogram!\n");
            exit(EXIT_FAILURE);
        }

        dends[i]->vertex = -1;
        dends[i]->left = NULL;
        dends[i]->right = NULL;
    }

    return dends;
}

// Frees a vertexSet struct
static void freeVertexSet(vertexSet *v) {
    free(v->vertices);
}

// Obtains the graph vertices contained within a Dendrogram
// Returns an empty vertexSet if no vertices are present
// Produces a runtime error if vertices array fails to be allocated
static vertexSet getVertices(Dendrogram d, int maxVertices) {
    vertexSet vSet = {0};
    vSet.vertices = (int *)malloc(maxVertices * sizeof(int));

    if (vSet.vertices == NULL) {
        fprintf(stderr, "Couldn't allocate vertices array!\n");
        exit(EXIT_FAILURE);
    }

    recurGetVertices(d, &vSet.nV, vSet.vertices);

    return vSet;
}

// Recursively traverse a Dendrogram to find all vertices it contains
static void recurGetVertices(Dendrogram d, int *nV, int *vertices) {
    // If vertex is not -1, current dendrogram is a leaf
    if (d->vertex > -1) {
        vertices[*nV] = d->vertex;
        (*nV)++;
        return;
    }

    recurGetVertices(d->left, nV, vertices);
    recurGetVertices(d->right, nV, vertices);
}
