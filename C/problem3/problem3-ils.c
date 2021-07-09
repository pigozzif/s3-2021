/* problem3-ils.c
 *
 * (C) 2021 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "s3problem.h"
#include "ils.h"

gsl_rng *rng;    /* The single rng instance used by the whole code */

int main(int argc, char **argv) {
    struct problem *p;
    struct solverState *ss;
    int max_iter, i;
    double mincost;
    double cost;
    double* costPtr = &cost;
    double* mincostPtr = &mincost;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <file name> <max iter>\n", argv[0]);
        return 0;
    }

    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937); 
    gsl_rng_set(rng, time(0));

    /* Input arguments */
    max_iter = atoi(argv[2]);

    /* Problem and solver instantiation */
    p = newProblem(argv[1], 1.0, 1.0);
    if (p != NULL) {
        ss = newSolver(p);

        /* Run */
        mincostPtr = getObjectiveVector(mincostPtr, getSolverSolution(ss));
        printf("iter = 0, obj = %.0f\n", getSolverSolution(ss)->score);
        printSolution(getSolverSolution(ss));
        for (i = 0; i < max_iter; i++) {
            printf("===ITER===\n");
            printf("===BEST SOLUTION SCORE===\n");
            printf("%lf\n", getSolverSolution(ss)->score);
            nextSolverState(ss);
            costPtr = getObjectiveVector(costPtr, getSolverSolution(ss));
            if (*costPtr < *mincostPtr) {
                mincost = cost;
                printf("iter = %d, obj = %.0f\n", i+1, *mincostPtr);
            }
        }

        /* Report result */
        printSolution(getSolverSolution(ss));

        /* Clean up */
        freeSolver(ss);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}

