#include <problem.h>

/**
*This file implements the functions defined in problem for the 
* computer network problem.
*/


struct Move{
    Problem * problem_instance;
    int node_concerned;
    int new_parent;
    double new_score;
};
