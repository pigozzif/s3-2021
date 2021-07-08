#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/**
  * This file implements the functions defined in problem for the 
  * computer network problem.
  */


/**
  * Basic data structures
  */
struct coordinate {
    double x, y;
};

struct problem{
    int n;  // number of buildings
    double trench_cost;
    double cable_cost;
    double **distances;  // distance matrix among buildings
    struct coordinate *coordinates;
};

struct solution{
	struct problem *problem_instance;
	double *paths_length_to_center; // total path to the center from node n
	double *lengths_to_parent;// the distance from each node in the tree to the parent. The parent is the immediate node going up to the root of the tree.
	int * parents; // the array of indexes of the parent for each node.
	double score; // the objective function value for the solution.
};

struct move{
    struct problem * problem_instance;
    int node_concerned;
    int new_parent;
    double new_score;
};

/**
  * Problem instantiation and inspection
  */
double euclidean_distance(double x1, double x2, double y1, double y2) {
     double x_diff = x1 - x2;
     double y_diff = y1 - y2;
     return sqrt(x_diff * x_diff + y_diff * y_diff);
}

// TODO: distance matrix can be optimised
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
         for (int j = 0; j < n_buildings; ++j) {
              new_problem->distances[i][j] = euclidean_distance(new_problem->coordinates[i].x, new_problem->coordinates[i].y, new_problem->coordinates[j].x, new_problem->coordinates[j].y);
         }
    }

    // free resources and return
    fclose(new_fp);
    return new_problem;
}

int getNumObjectives(const struct problem *p) {
    return 1;
}

void freeProblem(struct problem *p) {
    for (int i = 0; i < p->n; ++i) {
        free(p->distances[i]);
    }
    free(p->distances);
    free(p->coordinates);
    free(p);
}

struct solution *allocSolution(struct problem *p) {
    struct solution *sol = (struct solution *) malloc(sizeof(struct solution));
    if (sol == NULL) {
        return sol;
    }
    sol->problem_instance = p;
    sol->paths_length_to_center = (double *) malloc(sizeof(double) * p->n);
    if (sol->paths_length_to_center == NULL) {
        return NULL;
    }
    sol->lengths_to_parent = (double *) malloc(sizeof(double) * p->n);
    if (sol->lengths_to_parent == NULL) {
        return NULL;
    }
    sol->parents = (int *) malloc(sizeof(int) * p->n);
    if (sol->parents == NULL) {
        return NULL;
    }
    for (int i = p->n - 1; i >= 0; --i) {
        sol->parents[i] = -1;
    }
    return sol;
}

void freeSolution(struct solution *s) {
    free(s->paths_length_to_center);
    free(s->lengths_to_parent);
    free(s->parents);
    free(s);
}

struct move *allocMove(struct problem *p) {
    struct move *m = (struct move *) malloc(sizeof(struct move));
    if (m == NULL) {
        return m;
    }
    m->problem_instance = p;
    m->node_concerned = -1;
    m->new_parent = -1;
    m->new_score = 0.0;
    return m;
}

void freeMove(struct move *v) {
    free(v);
}

/**
  * Reporting
  */
void printProblem(struct problem *p) 
{
    int i, j;
    if(p == NULL)
        return;
    printf("Buildings: %d\n", p->n);
    for(i=0; i<p->n; i++)
        printf("B%d: (%.2f,%.2f); ",i, p->coordinates[i].x, p->coordinates[i].y); 
    printf("\nDistances (written into problem struct):\n\t");
    for(i=0; i<p->n; i++)
        printf("B%d\t",i);
    printf("\n");
    for(j=0; j<p->n; j++)
    {
        printf("B%d\t",j);
        for(i=0; i<p->n; i++)
            printf("%.2f\t",p->distances[j][i]);
        printf("\n");
    }
    /* calculated from coordinates */
    printf("Distances (calculated):\n\t");
    for(i=0; i<p->n; i++)
        printf("B%d\t",i);
    printf("\n");
    for(j=0; j<p->n; j++)
    {
        printf("B%d\t",j);
        for(i=0; i<p->n; i++)
            printf("%.2f\t",euclidean_distance(
                p->coordinates[j].x, p->coordinates[i].x,
                p->coordinates[j].y, p->coordinates[i].y
                )
                );
        printf("\n");
    }
}

void printSolution(struct solution *s) { 
}

void printMove(struct move* v) {
}

/* coppying src solution (struct) elements to dest solution (struct) elements */
struct solution * copySolution(struct solution *dest, const struct solution *src)
{
    int i, n;
    if(dest == NULL || src == NULL)
        return NULL;
    dest->problem_instance = src->problem_instance;
    dest->score = src->score;
    n = src->problem_instance->n;
    for(i=0; i<n; i++)
    {
        dest->paths_length_to_center[i] = src->paths_length_to_center[i];
        dest->lengths_to_parent[i] = src->lengths_to_parent[i];
        dest->parents[i] = src->parents[i];
    }
    /* or maybe
    memcpy(dest->paths_length_to_center, src->paths_length_to_center, n*sizeof(double));
    memcpy(dest->lengths_to_parent, src->lengths_to_parent, n*sizeof(double));
    memcpy(dest->parents, src->parents, n*sizeof(int));
    */
   return dest;
}


// JUST FOR DEBUGGING
int main(void) {
    struct problem *p = newProblem("buildings.txt", 1.0, 1.0);
    printProblem(p);
    struct solution *s = allocSolution(p);
    printSolution(s);    
    struct move* m = allocMove(p);
    printMove(m);
    freeSolution(s);
    freeMove(m);
    freeProblem(p);
    return 0;
}
