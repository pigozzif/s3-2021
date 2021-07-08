#include "problem.h"

/**
  * data structures
  */

struct problem;
struct solution;
struct move;



/**
  * problem instantiation and inspection
  */

struct problem *newProblem(char* filename);
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
