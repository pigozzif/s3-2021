#include "../solver/problem.h"

/**
  * data structures
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
};



/**
  * problem instantiation and inspection
  */

struct problem *newProblem(char* filename, double tc, double cc);
int getNumObjectives(const struct problem *p);



/**
  * memory management
  */

void freeProblem(struct problem *p);
struct solution *allocSolution(struct problem *p);
void freeSolution(struct solution *);
struct move *allocMove(struct problem *p);
void freeMove(struct move *v);


/**
  * reporting
  */

void printProblem(struct problem *p);
void printSolution(struct solution *s);
void printMove(struct move *v);



/**
  * operations on solutions
  */

struct solution *randomSolution(struct solution *s);
struct solution *copySolution(struct solution *dest, const struct solution *src);
double *getObjectiveVector(double *objv, struct solution *s);
struct solution *applyMove(struct solution *s, const struct move *v);
int getNeighbourhoodSize(struct solution *s);
struct solution *resetRandomMoveWOR(struct solution *s);



/**
  * operations on moves
  */

struct move *randomMove(struct move *v, const struct solution *s);
struct move *copyMove(struct move *dest, const struct move *src);
double *getObjectiveIncrement(double *obji, struct move *v, struct solution *s);
struct move *randomMoveWOR(struct move *v, struct solution *s);
