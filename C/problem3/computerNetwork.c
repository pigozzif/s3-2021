/**
*This file implements the functions defined in problem for the 
* computer network problem.
*/

struct problem {
    int n;  // number of buildings
    double trench_cost;
    double cable_cost;
    double **distances;  // distance matrix among buildings
};

struct solution{
	problem *problem_instance;
	double *paths_length_to_center; // total path to the center from node n
	double *lengths_to_parent;// the distance from each node in the tree to the parent. The parent is the immediate node going up to the root of the tree.
	int * parents; // the array of indexes of the parent for each node.
	double score; // the objective function value for the solution.
};

struct Move{
    Problem * problem_instance;
    int node_concerned;
    int new_parent;
    double new_score;
};
