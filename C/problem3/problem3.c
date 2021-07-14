//
// Created by Federico Pigozzi on 14/07/21.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/problem.h"


/**
 *
 * Utility functions
 */

double euclidean_distance(double x1, double x2, double y1, double y2) {
    double x_diff = x1 - x2;
    double y_diff = y1 - y2;
    return sqrt(x_diff * x_diff + y_diff * y_diff);
}

int randomInt(unsigned int b){
    return rand() % (b + 1);
}

/**
 *
 * Problem instantiation and inspection
 */

struct problem *newProblem(char *filename, double tc, double cc) {
    FILE* fp;
    struct problem *new_problem = malloc(sizeof(struct problem));
    if (new_problem == NULL) {
        return new_problem;
    }
    char ch;

    // open file and return NULL if fail
    fp = fopen(filename, "r");
    if (fp == NULL) {
        return NULL;
    }

    // get number of buildings (rows) from file
    int n_buildings = 0;
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') {
            ++n_buildings;
        }
        //printf("%c", ch);
    }
    fclose(fp);

    // reopen file and fill array with buildings' coordinates
    FILE* new_fp;
    new_fp = fopen(filename, "r");
    if (new_fp == NULL) {
        return NULL;
    }
    int row = 0;
    int i = 0;
    char buffer[2];
    new_problem->coordinates = (struct coordinate*) malloc(n_buildings * sizeof(struct coordinate));
    if (new_problem->coordinates == NULL) {
        return NULL;
    }
    while ((ch = fgetc(new_fp)) != EOF) {
        if (ch == '\n') {
            ++row;
            new_problem->coordinates[row].y = atof(buffer);
            i = 0;
            continue;
        }
        else if (ch == ',') {
            new_problem->coordinates[row].x = atof(buffer);
            i = 0;
            continue;
        }
        buffer[i++] = ch;
    }

    // update problem struct fields and compute distances
    new_problem->n = n_buildings;
    new_problem->trench_cost = tc;
    new_problem->cable_cost = cc;
    new_problem->distances = (double **) malloc(sizeof(double*) * n_buildings);
    if (new_problem->distances == NULL) {
        return NULL;
    }
    for (int i = 0; i < n_buildings; ++i) {
        new_problem->distances[i] = (double *) malloc(sizeof(double) * n_buildings);
        if (new_problem->distances[i] == NULL) {
            return NULL;
        }
    }
    for (int i = 0; i < n_buildings; ++i) {
        for (int j = 0; j < n_buildings; ++j) {
            new_problem->distances[i][j] = euclidean_distance(new_problem->coordinates[i].x, new_problem->coordinates[j].x, new_problem->coordinates[i].y, new_problem->coordinates[j].y);
        }
    }

    // free resources and return
    fclose(new_fp);
    return new_problem;
}

int getNumObjectives(const struct problem *p) {
    return 1;
}

/**
 *
 * Memory management
 */

void freeProblem(struct problem *p) {
    for (int i = 0; i < p->n; ++i) {
        free(p->distances[i]);
    }
    free(p->distances);
    free(p->coordinates);
    free(p);
}

struct solution *allocSolution(struct problem *p) {
    struct solution *s = (struct solution*) malloc(sizeof(struct solution));
    if (s == NULL) {
        return NULL;
    }
    s->problem_instance = p;
    s->parents = (int*) malloc(sizeof(int) * p->n);
    if (s->parents == NULL) {
        return NULL;
    }
    s->score = -1.0;
    return s;
}

void freeSolution(struct solution *s) {
    free(s->parents);
    free(s);
}

struct move *allocMove(struct problem *p) {
    struct move *m = (struct move *) malloc(sizeof(struct move));
    if (m == NULL) {
        return m;
    }
    m->problem_instance = p;
    m->source_concerned = -1;
    //m->target_concerned = -1;
    m->new_parent = -1;
    m->new_score = 0.0;
    return m;
}

void freeMove(struct move *v) {
    free(v);
}

/**
 *
 * Reporting
 */

void printProblem(struct problem *p) {
    printf("Buildings: %d\n", p->n);
    for (int i = 0; i < p->n; ++i) {
        printf("B%d: (%.2f,%.2f); ", i, p->coordinates[i].x, p->coordinates[i].y);
    }
    printf("\nDistances:\n\t");
    for (int i=0; i < p->n; ++i) {
        printf("B%d\t", i);
    }
    printf("\n");
    for (int j = 0; j < p->n; ++j) {
        printf("B%d\t",j);
        for(int i = 0; i < p->n; ++i) {
            printf("%.2f\t", p->distances[j][i]);
        }
        printf("\n");
    }
}

void printSolution(struct solution *s) {
    printf("Solution - %d buildings:\n", s->problem_instance->n);
    printf("Build\tParent\tDistance\tPath to center\n");
    //for(int i = 0; i < n; ++i) {
    //    printf("B%d\t%i\t%f\t%f\n",i,s->parents[i],s->lengths_to_parent[i],s->paths_length_to_center[i]);
    //}
    printf("Score: %f\n",s->score);
    printf("Neighborhood :\n");
}


void printMove(struct move* v) {
    printf("Move - %d buildings:\n", v->problem_instance->n);
    printf("source node: %d, new parent: %d, score: %lf\n", v->source_concerned, /*v->target_concerned,*/ v->new_parent, v->new_score);
}

/**
 *
 * Operations on solutions
 */

struct solution * randomSolution(struct solution *s){
    const int N = s->problem_instance->n;
    const int n = N-2;
    int * c = malloc(n * sizeof(int));
    int * d = calloc(N, sizeof(int));
    for(int i = 0; i < N-2; ++i) {
        const int v = rand() % N;
        c[i] = v;
        d[v]++;
    }
    int i = 0;
    while (d[i] != 0) {
        i++;
    }
    int leaf = i;
    for (int j = 0; j < n; ++j) {
        const int v = c[j];
        s->parents[leaf] = v;
        if (--d[v] == 0 && v<i){
            leaf = v;
        }
        else {
            i++;
            while(d[i]!=0){
                i++;
            }
            leaf = i;
        }
    }
    s->parents[leaf] = N-1;
    free(d);
    free(c);
    s = resetRandomMoveWOR(s);
    return s;
}

struct solution *copySolution(struct solution *dest, const struct solution *src) {
    memcpy(dest->parents, src->parents, sizeof(int) * src->problem_instance->n);
    dest->score = src->score;
    return dest;
}

double computeShortestPath(struct solution *s, int node, double temp) {
    if (s->parents[node] == -1) {
        return temp;
    }
    temp += s->problem_instance->distances[node][s->parents[node]];
    return computeShortestPath(s, s->parents[node], temp);
}

double *getObjectiveVector(double *objv, struct solution *s) {
    double total_trench = 0.0;
    for (int i = 0; i < s->problem_instance->n; ++i) {
        total_trench +=  s->problem_instance->distances[i][s->parents[i]];
    }
    double total_path = 0.0;
    for (int i = 0; i < s->problem_instance->n; ++i) {
        total_path += computeShortestPath(s, i, 0.0);
    }
    double tot = s->problem_instance->trench_cost * total_trench + s->problem_instance->n * total_path;
    *objv = tot;
    s->score = tot;
    return objv;
}

struct solution *applyMove(struct solution *s, const struct move *v) {
    s->parents[v->source_concerned] = v->new_parent;
    s->score = v->new_score;
    return s;
}

int getNumberOfMovesPerNode(struct solution *s, int node) {
    int ans = 0;
    for (int i = 0; i < s->problem_instance->n; ++i) {
        if (s->parents[i] == node) {
            ans += 1;
        }
    }
    if (s->parents[node] != -1 && s->parents[s->parents[node]] != -1) {
        ans += 1;
    }
    return ans;
}

int getNeighbourhoodSize(struct solution *s) {
    int ans = 0;
    for (int node = 0; node < s->problem_instance->n; ++node) {
        ans += getNumberOfMovesPerNode(s, node);
    }
    return ans;
}

struct solution *resetRandomMoveWOR(struct solution *s) {

}

/**
 *
 * Operations on moves
 */

struct move* randomMove(struct move *v, const struct solution *s) {
    int chosen_move = randomInt(getNeighbourhoodSize(s));
    //printf("chosen move: %d\n", chosen_move);
    int chosen_node = -1;
    int chosen_new_parent = -1;
    int cum_sum = 0;
    int stop = 0;
    //printChildren(s);
    // iterate over all nodes
    for (int i = 0; i < s->problem_instance->n; ++i) {
        if (stop) {
            break;
        }
        // iterate over all possible moves from node
        for (int j = 0; j < getNumberOfMovesPerNode(s, i); ++j) {
            if (cum_sum >= chosen_move) {
                chosen_node = i;
                // if it must be the grandparent
                if (j == getNumberOfMovesPerNode(s, i) - 1 && s->parents[s->parents[i]] != -1) {
                    chosen_new_parent = s->parents[s->parents[i]];
                }
                // else it is one of the children
                else {
                    for (int k = 0; k < s->problem_instance->n; ++k) {
                        // check it is a child of the parent, and it is not the node under consideration
                        if (s->parents[k] == s->parents[i] && k != i) {
                            ++k;
                        }
                        // if it is what we want
                        if (k == j) {
                            chosen_new_parent = k;
                            break;
                        }
                    }
                }
                stop = 1;
                break;
            }
            ++cum_sum;
        }
    }
    //printf("chosen node: %d\n", chosen_node);
    //printf("chosen new parent: %d\n", chosen_new_parent);
    // update move
    v->source_concerned = chosen_node;
    //v->target_concerned = s->parents[chosen_node];
    v->new_parent = chosen_new_parent;
    return v;
}

struct move *copyMove(struct move *dest, const struct move *src) {
    dest->problem_instance = src->problem_instance;
    dest->source_concerned = src-> source_concerned;
    //dest->target_concerned = src->target_concerned;
    dest->new_parent = src-> new_parent;
    dest->new_score = src-> new_score;
    return dest;
}

double *getObjectiveIncrement(double *obji, struct move *v, struct solution *s) {
    struct solution *copy = allocSolution(s->problem_instance);
    copy = copySolution(copy, s);
    applyMove(copy, v);
    double increment = *getObjectiveVector(obji, copy);
    v->new_score = increment;
    *obji = increment;
    return obji;
}

int main(void) {
    return 0;
}

