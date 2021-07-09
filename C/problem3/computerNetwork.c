#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/s3problem.h"

/**
  * This file implements the functions defined in problem for the 
  * computer network problem.
  */


/**
  * Basic data structures
  */
/*struct coordinate {
    double x, y;
};

struct problem{
    int n;  // number of buildings
    double trench_cost;
    double cable_cost;
    double **distances;  // distance matrix among buildings
    struct coordinate *coordinates;
};

struct move{
    struct problem * problem_instance;
    int source_concerned;
    int target_concerned;
    int new_parent;
    double new_score;
};

struct neighborhood{
    int *randomSample;
    int *moves; // [from0,to0,from1,to1 ....]
    int maxSize;
    int sampleSize;
};

struct solution{
	struct problem *problem_instance;
	struct neighborhood * neighborhood;
	double *paths_length_to_center; // total path to the center from node n
	double *lengths_to_parent;// the distance from each node in the tree to the parent. The parent is the immediate node going up to the root of the tree.
	int * parents; // the array of indexes of the parent for each node.
	double score; // the objective function value for the solution.
    int **children;
    int *nbChildren;
};*/

/**
  * Problem instantiation and inspection
  */
double euclidean_distance(double x1, double x2, double y1, double y2) {
     double x_diff = x1 - x2;
     double y_diff = y1 - y2;
     return sqrt(x_diff * x_diff + y_diff * y_diff);
}


/**
* returns a random integer between a and b.
*/
int randomInt(unsigned int b){
	return rand()%(b+1);
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
    //fclose(fp);

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
    //fclose(new_fp);
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

struct neighborhood *allocNeighborhood(const struct solution* s){
    struct neighborhood * neigh = (struct neighborhood *) malloc(sizeof(struct neighborhood));
    if (neigh == NULL) {
        return neigh;
    }
    const int n = s->problem_instance->n;
    neigh->randomSample = malloc(n*n*sizeof(int));
    neigh->moves = malloc(n*n*sizeof(int));
    return neigh;
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
    sol->neighborhood = allocNeighborhood(sol);
    return sol;
}


void freeNeighborhood(struct neighborhood * n){
    free(n->randomSample);
    free(n->moves);
    free(n);
}

void freeSolution(struct solution *s) {
    return;
    /*free(s->paths_length_to_center);
    free(s->lengths_to_parent);
    free(s->parents);
    for (int i = 0; i < s->problem_instance->n; ++i) {
        free(s->children[i]);
    }
    free(s->children);
    free(s->nbChildren);
    if (s->neighborhood != NULL) {
        freeNeighborhood(s->neighborhood);
    }
    free(s);*/
}

struct move *allocMove(struct problem *p) {
    struct move *m = (struct move *) malloc(sizeof(struct move));
    if (m == NULL) {
        return m;
    }
    m->problem_instance = p;
    m->source_concerned = -1;
    m->target_concerned = -1;
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

void printNeighborhood(struct neighborhood* n) 
{
    if(n == NULL)
        return;
    for (int idx=0; idx<n->sampleSize; idx++) {
        printf("source node: %d, new parent: %d\n", n->moves[idx*2], n->moves[idx*2+1]);
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
    for(i=0; i<n; i++) {
        printf("B%d\t%i\t%f\t%f\n",i,s->parents[i],s->lengths_to_parent[i],s->paths_length_to_center[i]);
    }
    printf("Score: %f\n",s->score);
    printf("Neighborhood :\n");
    printNeighborhood(s->neighborhood);
}


void printMove(struct move* v) 
{
    int n;
    if(v == NULL || v->problem_instance == NULL)
        return;
    n = v->problem_instance->n;
    printf("Move - %d buildings:\n",n);


    printf("source node: %d, target node: %d, new parent: %d, score: %lf\n", v->source_concerned, v->target_concerned, v->new_parent, v->new_score);

    printf("source node: %d, target node: %d, new parent: %d, score: %lf\n", v->source_concerned, v->target_concerned, v->new_parent, v->new_score);


    printf("source node: %d, target node: %d, new parent: %d, score: %lf\n", v->source_concerned, v->target_concerned, v->new_parent, v->new_score);
    printf("source node: %d, target node: %d, new parent: %d, score: %lf\n", v->source_concerned, v->target_concerned, v->new_parent, v->new_score);
}

struct neighborhood *copyNeighborhood(struct neighborhood* dest, const struct neighborhood *src, int n) {
    memcpy(dest->randomSample, src->randomSample, n*n*sizeof(int));
    memcpy(dest->moves, src->moves, n*n*sizeof(int));
    dest->maxSize = src->maxSize;
    dest->sampleSize = src->sampleSize;
    return dest;
}

struct solution *copySolution(struct solution *dest, const struct solution *src) {
    if (dest == NULL || src == NULL) {
        return NULL;
    }
    int n = dest->problem_instance->n;
    dest->problem_instance = src->problem_instance;
    memcpy(dest->paths_length_to_center, src->paths_length_to_center, n * sizeof(double));
    memcpy(dest->lengths_to_parent, src->lengths_to_parent, n * sizeof(double));
    memcpy(dest->parents, src->parents, n * sizeof(int));
    memcpy(dest->children, src->children, n * sizeof(int*));
    for (int i = 0; i < n; ++i) {
        memcpy(&dest->children[i], &src->children[i], n * sizeof(int));
    }
    memcpy(dest->nbChildren, src->nbChildren, n * sizeof(int));
    dest->score = src->score;
    struct neighborhood *ne = allocNeighborhood(src);
    ne = copyNeighborhood(ne, src->neighborhood, n);
    dest->neighborhood = ne;
    return dest;
}

// TO TEST


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

double * getObjectiveVector(double *objv, struct solution * s){
    const int root = s->problem_instance->n-1;
    double trenches =0.0;
    double cables = 0.0;
    recursiveObjectiveVector(&trenches,&cables,root,s);
    *objv = trenches*s->problem_instance->trench_cost + cables*s->problem_instance->cable_cost;
    s->score = *objv;
    return objv;
}

void recursivelyUpdateNode(struct solution *s, int node) {
    s->lengths_to_parent[node] = s->problem_instance->distances[node][s->parents[node]];
    s->paths_length_to_center[node] = s->paths_length_to_center[s->parents[node]] + s->lengths_to_parent[node];    
    for (int i = 0; i < s->nbChildren[node]; ++i) {
        const int child = s->children[node][i];
        if (s->nbChildren[child] > 0){
            recursivelyUpdateNode(s, child);
        }
    }
}

struct solution *applyMove(struct solution *s, const struct move *v) {
    s->parents[v->source_concerned] = v->new_parent;
    s->score = v->new_score;
    int start_copying = 0;
    for (int i = 0; i < s->problem_instance->n; ++i) {
         if (s->children[v->target_concerned][i] == v->source_concerned) {
             start_copying = 1;
         }
         if (start_copying == 1 && i != s->problem_instance->n) {
             s->children[v->target_concerned][i] = s->children[v->target_concerned][i];    
         }
    }
    s->nbChildren[v->target_concerned] -= 1;
    s->nbChildren[v->new_parent] += 1;
    recursivelyUpdateNode(s, v->source_concerned);
    s = resetRandomMoveWOR(s);
    //printSolution(s);
    return s; 
}


// TO TEST 


struct move* randomMoveWOR(struct move *v, struct solution *s){
    struct neighborhood * n_view = s->neighborhood;
    if(n_view->sampleSize<1){
        return NULL;
    }
    const int rand_res = rand()%n_view->sampleSize;
    const int idx = n_view->randomSample[rand_res];
    v->source_concerned = n_view->moves[idx*2];
    v->new_parent = n_view->moves[idx*2+1];
    n_view->randomSample[rand_res] = n_view->randomSample[--n_view->sampleSize];
    return v;
}

struct solution *resetRandomMoveWOR(struct solution *s){
    struct neighborhood * n_view = s->neighborhood;
    int idx =0 ;
    for(int parent=0;parent<s->problem_instance->n-1;parent++){
        for(int from_idx=0;from_idx<s->nbChildren[parent];from_idx++){
            const int from = s->children[parent][from_idx];
            for (int to_idx=0; to_idx<s->nbChildren[parent]; to_idx++) {
                const int to = s->children[parent][to_idx] ;
                if(from!=to){
                    n_view->moves[idx++]=from;
                    n_view->moves[idx++]=to;
                }
            }
            n_view->moves[idx++]=from;
            n_view->moves[idx++]=s->parents[parent];
        }
    }
    const int root = s->problem_instance->n-1;
    for(int from_idx=0;from_idx<s->nbChildren[root];from_idx++){
        const int from = s->children[root][from_idx];
        for (int to_idx=0; to_idx>s->nbChildren[root]; to_idx++) {
            const int to = s->children[root][to_idx] ;
            if(from!=to){
                n_view->moves[idx++]=from;
                n_view->moves[idx++]=to;
            }
        }
    }
    n_view->sampleSize = n_view->maxSize=idx/2;
    for(int i=0;i<n_view->sampleSize;i++){
        n_view->randomSample[i]=i;
    }
    return s;
}

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
    s = resetRandomMoveWOR(s);
    return s;
}

int getNeighbourhoodSize(struct solution *s) {
    int ans = 0;
    for (int node = 0; node < s->problem_instance->n - 1; ++node) {
         ans += s->nbChildren[s->parents[node]] - 1;
         if (s->parents[s->parents[node]] != -1) {
             ans += 1;
         }
    }
    return ans;
}

/**
* Allocate the memory for the neighborhood.
*
*/
void allocateNeighborhood(struct solution *s){

}



/**
*
* Chechs if the solution generated by the move v is
* in the neighborhood of solution s.
*/
int inNeighbors(struct move *v, struct solution *s){


}

/**
*
* Chechs if the solution generated by the move v is
* in the neighborhood of a solution s.
*/
void addNeighbor(struct move *v, struct solution *s){


}





/**
 * randomMove() implements uniform random sampling of the neighbourhood of a
 * given solution, with replacement.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place.
 */
struct move *randomMove(struct move *v, const struct solution *s){

	int solution_length = s->problem_instance->n;
	// edges to delete and to add.
	int edge_to_delete, edge_to_add;
	// the degree of each node in the tree
	int node_degree[solution_length];

	// for each node in the solution computes its degree.
	int sum_degrees = 0;
	int parent;

	for (unsigned int i = 0; i < solution_length - 1 ; ++i){
		parent = s->parents[i];
		node_degree[i] = s->nbChildren[parent];
		sum_degrees += node_degree[i];
	}

	// proportionate selection of the edege to be deleted.
	edge_to_delete = randomInt(sum_degrees-1);
	int node_index = 0;
	int i_sum_degree = node_degree[node_index];

	while(edge_to_delete > i_sum_degree ){
		node_index++;
		i_sum_degree += node_degree[node_index];
	}
	
	edge_to_delete = node_index;
	
	// new parent
	parent = s->parents[edge_to_delete];
	// degree of this parent
	int parent_degree = s->nbChildren[parent];

	int upper_bound = (s->parents[parent] == -1) ? parent_degree : parent_degree + 1;
	int index_edge_to_add = randomInt(upper_bound);
	// 0 is the for representing the index of the 
	// parent. That is the move does not change the 
	// solution.
	while (s->children[parent][index_edge_to_add] == edge_to_delete){
		index_edge_to_add = randomInt(upper_bound);
	}
	
	// new parent is the parent of the parent.
	if (index_edge_to_add == parent_degree + 1 && upper_bound == parent_degree + 1) {
		edge_to_add = s->parents[parent];
	}else{
		edge_to_add = s->children[parent][index_edge_to_add];
	}

	// update the move
	v->source_concerned = edge_to_delete;
	v->target_concerned = parent;
	v->new_parent = edge_to_add;
	// v->new_score = getObjectiveIncrement(obji, v, s);

	return v;
}




  /**************************/
 /* Operations on moves*/
/**************************/

/* copyMove() copies the contents of the second argument to the first
 * argument, which must have been previously allocated with allocMove().
 */
struct move *copyMove(struct move *dest, const struct move *src) {
    dest->problem_instance = src->problem_instance;
    dest->source_concerned = src-> source_concerned;
    dest->target_concerned = src->target_concerned;
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
    struct solution *copy = allocSolution(s->problem_instance);
    copy = copySolution(copy, s);
    copy = applyMove(copy, v);
    double increment = 0.0;
    for (int i = 0; i < copy->problem_instance->n; ++i) {
        increment += copy->problem_instance->trench_cost * copy->lengths_to_parent[i] + copy->problem_instance->cable_cost * copy->paths_length_to_center[i];
    }
    v->new_score = increment;
    *obji = increment;
    return obji;
}

// JUST FOR DEBUGGING
/*int main(void) {
    struct problem *p = newProblem("buildings.txt", 1.0, 1.0);
    printProblem(p);
    struct solution *s = allocSolution(p);
    s = randomSolution(s);
    printSolution(s);
    struct solution *s2 = allocSolution(p);
    s2 = copySolution(s2, s);
//    printf("%d", getNeighbourhoodSize(s));
    double cost;
    double* costPtr = &cost;
    getObjectiveVector(costPtr, s2);
    printSolution(s2);
    printf("score %lf\n", *costPtr);
    struct move* m = allocMove(p);
    //m->source_concerned = 0;
    //m->target_concerned = 2;
    //m->new_parent = 6;
    //obj = *getObjectiveIncrement(&obj, m, s2);
    //printf("%lf\n", obj);
    //s = applyMove(s2, m);
    printSolution(s2);
    //struct move* m2 = allocMove(p);
    //m2 = copyMove(m2, m);
    //printMove(m);
    //printMove(m2);
    //freeSolution(s);
    //freeSolution(s2);
    //freeMove(m);
    //freeMove(m2);
    freeProblem(p);
    return 0;
}*/
