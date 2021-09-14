// Lance-Williams Algorithm for Hierarchical Agglomerative Clustering
// COMP2521 Assignment 2

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "LanceWilliamsHAC.h"

#define INFINITY DBL_MAX
#define TRUE 1
#define FALSE 0

#define DEBUG FALSE

typedef struct targets {
    int a;
    int b;
} targets;

typedef struct vertexSet {
    int *vertices;
    int nV;
} vertexSet;

void freeDendrogram(Dendrogram d);
double dMax(double a, double b);
double dMin(double a, double b);
double clusterDist(Graph g, vertexSet vSet1, vertexSet vSet2, int method);
double vertexDist(Graph g, targets t);
double **AllocDistArray(int size);
void FreeDistArray(double **dist, int size);
void mergeDends(Dendrogram *dends, int size, targets tDends);
Dendrogram *AllocDends(int size);
void FreeDends(Dendrogram *dends, int size);
void freeVertexSet(vertexSet *v);
vertexSet getVertices(Dendrogram d, int maxVertices);
void RecurGetVertices(Dendrogram d, int *nV, int *vertices);

/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {
    if (DEBUG == TRUE) {
        printf("Starting\n");
        fflush(stdout);
    }
    // Handle invalid method value input
    if (method != SINGLE_LINKAGE && method != COMPLETE_LINKAGE) {
        fprintf(stderr, "Invalid method int value (%d) given to LanceWilliamsHAC function.\n", method);
        exit(EXIT_FAILURE);
    }

    int nV = GraphNumVertices(g);

    if (DEBUG == TRUE) {
        printf("Allocating Dendrogram array\n");
        fflush(stdout);
    }

    // Initial number of dendrograms is number of graph vertices
    Dendrogram *dends = AllocDends(nV);
    for (int i = 0; i < nV; i++) {
        if (DEBUG == TRUE) {
            printf("Allocating dend for vertex %d\n", i);
            fflush(stdout);
        }
        dends[i]->vertex = i;
    }

    if (DEBUG == TRUE) {
        printf("Finished allocating dends, starting HAC\n");
        fflush(stdout);
    }

    // Iterate n times to produce HAC dendrogram of given graph
    for (int numDends = nV; numDends > 1; numDends--) {
        int iterNum = nV - numDends + 1;
        if (DEBUG == TRUE) {
            printf("%d: Starting iteration %d\n", iterNum, iterNum);
            fflush(stdout);
        }
        targets dendsToMerge = { 0, 1 };
        double minClusterDist = INFINITY;

        for (int i = 0; i < numDends; i++) {
            if (DEBUG == TRUE) {
                printf("%d: Getting vertices of 1st dendrogram %d\n", iterNum, i);
                fflush(stdout);
            }
            vertexSet vSet1 = getVertices(dends[i], nV);

            if (DEBUG == TRUE) {
                printf("Vertices of 1st dendrogram are: ");
                for (int k = 0; k < vSet1.nV; k++) {
                    printf("%d ", vSet1.vertices[k]);
                }
                printf("\n");
                fflush(stdout);
            }

            for (int j = i + 1; j < numDends; j++) {
                if (DEBUG == TRUE) {
                    printf("%d: Getting vertices of 2nd dendrogram %d\n", iterNum, j);
                    fflush(stdout);
                }
                vertexSet vSet2 = getVertices(dends[j], nV);

                if (DEBUG == TRUE) {
                    printf("Vertices of 2nd dendrogram are: ");
                    for (int k = 0; k < vSet2.nV; k++) {
                        printf("%d ", vSet2.vertices[k]);
                    }
                    printf("\n");
                    fflush(stdout);
                }

                if (DEBUG == TRUE) {
                    printf("%d: \tGetting distance between dendrograms %d and %d\n", iterNum, i, j);
                    fflush(stdout);
                }
                double cDist = clusterDist(g, vSet1, vSet2, method);
                if (DEBUG == TRUE) {
                    printf("%d: \tDistance was %f\n", iterNum, cDist);
                    fflush(stdout);
                }
                freeVertexSet(&vSet2);
                if (cDist < minClusterDist) {
                    dendsToMerge.a = i;
                    dendsToMerge.b = j;
                    minClusterDist = cDist;
                    if (DEBUG == TRUE) {
                        printf("%d: New minimum distance of %f!\n", iterNum, minClusterDist);
                        fflush(stdout);
                    }
                }
            }
            freeVertexSet(&vSet1);
        }

        if (DEBUG == TRUE) {
            printf("%d: Starting to merge dends %d and %d\n", iterNum, dendsToMerge.a, dendsToMerge.b);
            fflush(stdout);
        }
        mergeDends(dends, numDends, dendsToMerge);
        if (DEBUG == TRUE) {
            printf("%d: Finished merging dends\n", iterNum);
            fflush(stdout);
        }
        if (DEBUG == TRUE) {
            printf("%d: Finished iteration %d\n", iterNum, iterNum);
            fflush(stdout);
        }
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

double dMax(double a, double b) {
    return (a > b) ? a : b;
}

double dMin(double a, double b) {
    return (a < b) ? a : b;
}

double clusterDist(Graph g, vertexSet vSet1, vertexSet vSet2, int method) {
    double cDist = INFINITY;

    for (int i = 0; i < vSet1.nV; i++) {
        for (int j = 0; j < vSet2.nV; j++) {
            targets tars = { vSet1.vertices[i], vSet2.vertices[j] };

            double vDist = vertexDist(g, tars);
            if (method == SINGLE_LINKAGE) {
                cDist = dMin(cDist, vDist);
            } else if (method == COMPLETE_LINKAGE) {
                if (vDist != INFINITY) {
                    if (cDist == INFINITY) {
                        cDist = vDist;
                    } else {
                        cDist = dMax(cDist, vDist);
                    }
                }
            }
        }
    }

    return cDist;
}

double vertexDist(Graph g, targets t) {
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

// Allocates the distance 2D double array of a given size to be used in an
// iteration of the Lance-Williams formula
// Produces runtime errors if the array or sub-arrays fail to be allocated
double **AllocDistArray(int size) {
    // Allocate dist array
    double **dist = (double **)malloc(size * sizeof(double *));

    if (dist == NULL) {
        fprintf(stderr, "Couldn't allocate dist array!\n");
        exit(EXIT_FAILURE);
    }

    // Allocate dist sub-arrays
    for (int i = 0; i < size; i++) {
        dist[i] = (double *)malloc(size * sizeof(double));

        if (dist[i] == NULL) {
            fprintf(stderr, "Couldn't allocate dist sub-array!\n");
            exit(EXIT_FAILURE);
        }
    }

    return dist;
}

// Frees the distance 2D double array used in an iteration of the
// Lance-Williams formula
void FreeDistArray(double **dist, int size) {
    // Free dist sub-arrays
    for (int i = 0; i < size; i++) {
        free(dist[i]);
    }

    // Free dist array
    free(dist);
}

void mergeDends(Dendrogram *dends, int size, targets tDends) {
    assert(tDends.a > -1 && tDends.b > -1);
    if (DEBUG == TRUE) {
        printf("Starting mergeDends\n");
        fflush(stdout);
    }
    
    // Combine target Dendrograms into a single Dendrogram, add it to dends
    Dendrogram combDend = (Dendrogram)malloc(sizeof(DNode));

    if (combDend == NULL) {
        fprintf(stderr, "Couldn't allocate Dendrogram!\n");
        exit(EXIT_FAILURE);
    }

    if (DEBUG == TRUE) {
        printf("Target a vertex: %d\n", (*dends[tDends.a]).vertex);
        fflush(stdout);
        printf("Target b vertex: %d\n", (*dends[tDends.b]).vertex);
        fflush(stdout);
        printf("Starting to create combDend\n");
        fflush(stdout);
    }

    combDend->vertex = -1;
    combDend->left = dends[tDends.a];
    combDend->right = dends[tDends.b];
    dends[tDends.a] = combDend;
    dends[tDends.b] = NULL;

    if (DEBUG == TRUE) {
        printf("Finished creating combDend\n");
        fflush(stdout);
    }

    // Shift all dendrograms after the empty element down
    int foundEmpty = FALSE;
    for (int i = 0; i < size; i++) {
        if (foundEmpty == TRUE) {
            dends[i-1] = dends[i];
        } else if (dends[i] == NULL) {
            foundEmpty = TRUE;
        }
    }
}

// Allocates the dends Dendrogram array used in an iteration of the
// Lance-Williams formula
// Produces a runtime error if dends fails to be allocated
Dendrogram *AllocDends(int size) {
    Dendrogram *dends = (Dendrogram *)malloc(size * sizeof(Dendrogram));

    if (dends == NULL) {
        fprintf(stderr, "Couldn't allocate dends array!\n");
        exit(EXIT_FAILURE);
    }

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

// Frees the dends Dendrogram array used in an iteration of the
// Lance-Williams formula
void FreeDends(Dendrogram *dends, int size) {
    for (int i = 0; i < size; i++) {
        freeDendrogram(dends[i]);
    }

    free(dends);
}

void freeVertexSet(vertexSet *v) {
    free(v->vertices);
}

vertexSet getVertices(Dendrogram d, int maxVertices) {
    vertexSet vSet = {0};
    vSet.vertices = (int *)malloc(maxVertices * sizeof(int));

    if (vSet.vertices == NULL) {
        fprintf(stderr, "Couldn't allocate vertices array!\n");
        exit(EXIT_FAILURE);
    }

    RecurGetVertices(d, &vSet.nV, vSet.vertices);

    return vSet;
}

void RecurGetVertices(Dendrogram d, int *nV, int *vertices) {
    // If vertex is not -1, current dendrogram is a leaf
    if (d->vertex > -1) {
        vertices[*nV] = d->vertex;
        (*nV)++;
        return;
    }

    RecurGetVertices(d->left, nV, vertices);
    RecurGetVertices(d->right, nV, vertices);
}