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
        int **children;
        int *nbChildren;
};

struct move{
    struct problem * problem_instance;
    int node_concerned;
    int new_parent;
    double new_score;
};


/**
* This is a basic implementation of the neighborhood.
* it uses a list for maintaing the visited solutions in the
* neighborhood of a solution x.
* Note that this should be implemented as a set using
* some tree implementation, so we will have a O(log(n)) compelxity
* for the existance request.
* But for time reseason, we just use an array.
* TODO: implement a neighborhood using a set.
*
* Neighborhood definition:
*
*/
struct neighborhood {
	struct solution * x; // current solution
	unsigned int neighborhood_size;
	struct solution * neihborgs; // array of sampled neighbors.
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
    sol->nbChildren = (int*) malloc(sol->problem_instance->n * sizeof(int));
    for (int i = 0; i < sol->problem_instance->n; ++i) {
        sol->nbChildren[i] = 0;
    }
    sol->children = (int**) malloc(sol->problem_instance->n * sizeof(int*));
    for (int i = 0; i < sol->problem_instance->n; ++i) {
        sol->children[i] = (int*) malloc(sol->problem_instance->n * sizeof(int));
    }
    for (int i = 0; i < sol->problem_instance->n; ++i) {
        for (int j = 0; j < sol->problem_instance->n; ++j) {
            sol->children[i][j] = -1;
        }
    }
    return sol;
}

void freeSolution(struct solution *s) {
    free(s->paths_length_to_center);
    free(s->lengths_to_parent);
    free(s->parents);
    for (int i = 0; i < s->problem_instance->n; ++i) {
        free(s->children[i]);
    }
    free(s->children);
    free(s->nbChildren);
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
    printf("\nDistances:\n\t");
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
}

void printSolution(struct solution *s) 
{ 
    int i, n;
    if(s == NULL || s->problem_instance == NULL)
        return;
    n = s->problem_instance->n;
    printf("Solution - %d buildings:\n",n);
    printf("Build\tParent\tDistance\tPath to center\n");
    for(i=0; i<n; i++)
    {
        printf("B%d\t%i\t%f\t%f\n",i,s->parents[i],s->lengths_to_parent[i],s->paths_length_to_center[i]);
    }
    printf("Score: %f\n",s->score);
}

void printMove(struct move* v) 
{
    int n;
    if(v == NULL || v->problem_instance == NULL)
        return;
    n = v->problem_instance->n;
    printf("Move - %d buildings:\n",n);
    printf("Node: %d, new parent: %d, score: %lf\n", v->node_concerned, v->new_parent, v->new_score);
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



// TO TEST
/*
struct solution * randomSolution(struct solution *s){
    const int N = s->problem_instance->n; 
    const int n = N-2;
    int * c = malloc(n*sizeof(int));
    int * d = calloc(N, sizeof(int));
    for(int i=0;i<N-2;i++){
        const int v = rand()%N;
        c[i] = v;
        d[v]++;
    }
    int i = 0;
    while(d[i]!=0){
        i++;
    }
    int leaf = i;
    for (int j=0;j<n; j++) {
        const int v = c[j];
        s->parents[leaf] = v;
        s->children[v][s->nbChildren[v]++] = leaf;
        if(--d[v] == 0 && v<i){
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
    s->children[N-1][s->nbChildren[N-1]++] = leaf;
    free(d);
    free(c);
    return s;
}
*/

// TO TEST
/*
void recursiveObjectiveVector(double * trenches, double * cables, int node, struct solution * s){
    for(int i=0;i<s->nbChildren[node];i++){
        const int child = s->children[node][i];
        *trenches += s->lengths_to_parent[child] = s->problem_instance->distances[child][node];
        *cables += s->paths_length_to_center[child] = s->paths_length_to_center[node] + s->lengths_to_parent[child];
        if(s->nbChildren[child] > 0){
            recursiveObjectiveVector(trenches , cables, child, s);
        }
    }
}
*/

// TO TEST
/*
double * getObjectiveVector(double *objv, struct solution * s){
    const int root = s->problem_instance->n-1;
    double trenches =0.0;
    double cables = 0.0;
    recursiveObjectiveVector(&trenches,&cables,root,s);
    *objv = trenches*s->problem_instance->trench_cost + cables*s->problem_instance->cable_cost;
    return objv;
}
*/


// TO TEST
/*
struct neighborhood{
    int *randomSample;
    int *moves_decomposed; // [from0,to0,from1,to1 ....]
    int *maxSize;
    int *sampleSize;
}



struct move* randomMoveWOR(struct move *v, struct solution *s){
    struct neighborhood * n_view = s->neighborhood;
    if(n_view->sampleSize<1){
        v->node_concerned=-1;
        v->new_parent=-1;
        return v // THIS MOVE MUST NOT BE APPLIED
    }
    const int rand_res = rand()%n_view->sampleSize;
    const int idx = n_view->randomSample[rand_res];
    v->new_parent = n_view->moves_decomposed[idx*2];
    v->node_concerned = n_view->moves_decomposed[idx*2+1];
    n_view->randomSample[idx]=n_view->randomSample[sampleSize--];
    return v;
}
*/

/**
* Allocate the memory for the neighborhood.
*
*/
void allocateNeighborhood(struct solution *s){

}


/**
*
* Free memory allocated for the neighborhood n.
*/
void freeNeighborhood(struct neighborhood * n){

}


/**
*
* Chechs if the solution generated by the move v is
* in the neighborhood of solution s.
*/
bool inNeighbors(struct move *v, struct solution *s){


}

/**
*
* Chechs if the solution generated by the move v is
* in the neighborhood of a solution s.
*/
bool addNeighbor(struct move *v, struct solution *s){


}


/**
* returns a random integer between a and b.
*/
unsigned int randomInt(unsigned int a, unsigned int b){
	return 0;
}


/* randomMove() implements uniform random sampling of the neighbourhood of a
 * given solution, with replacement.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place.
 */
struct move *randomMove(struct move *v, const struct solution *s){
	// edges to delete and to add.
	unsigned int edge_to_delete, edge_to_add = 0;
	
	unsigned solution_length = solution->problem_instance->n;

	edge_to_delete = randomInt(1, solution_length);
	unsigned parent = solution->parents[edge_to_delete];

	x







	return v;
}




  /**************************/
 /* Operations on moves*/
/**************************/

/* copyMove() copies the contents of the second argument to the first
 * argument, which must have been previously allocated with allocMove().
 */
struct move *copyMove(struct move *dest, const struct move *src){
	dest->problem_instance = src->problem_instance;
    dest->node_concerned = src-> node_concerned;
    dest->new_parent = src-> new_parent;
    dest->new_score = src-> new_score;
    
	return dest;
}




/* getObjectiveIncrement() supports single and multiple objective move
 * evaluation with respect to the solution for which it was generated, before
 * it is actually applied to that solution (if it ever is). The result of
 * evaluating a move with respect to a solution other than that for which it
 * was generated (or to a pristine copy of it) is undefined.
 * Once a move is evaluated, results may be cached in the move itself, so
 * that they can be used by applyMove() to update the evaluation state of
 * the solution more efficiently.
 * In addition, results may also be cached in the solution in order
 * to speed up evaluation of future moves. Consequently, neither formal
 * argument is const.
 */
double *getObjectiveIncrement(double *obji, struct move *v, struct solution *s){
	return 0.0;
}



/* randomMoveWOR() implements uniform random sampling of the neighbourhood of
 * a given solution, without replacement.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The function returns this
 * pointer if a new move is generated or NULL if there are no moves left.
 */
struct move *randomMoveWOR(struct move *v, struct solution *s){
	return NULL
}

// JUST FOR DEBUGGING
int main(void) {
    struct problem *p = newProblem("buildings.txt", 1.0, 1.0);
    printProblem(p);
    struct solution *s = allocSolution(p);
    //s = randomSolution(s);
    printSolution(s);    
    struct move* m = allocMove(p);
    printMove(m);
    freeSolution(s);
    freeMove(m);
    freeProblem(p);
    return 0;
}
