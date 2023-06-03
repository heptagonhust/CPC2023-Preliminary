#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<sys/time.h>

#include "pcg.h"


#define START_TIME  0
#define STOP_TIME   2.0
#define DELTA_T     0.01
#define MAX_ITERS   300
#define tolerance   1.0e-6
#define normfactor  1.0e-4


int main(int argc, char *argv[]) {
	struct timeval tv;
	double total_time_cost = .0;
	double pcg_time_cost[3] = {0.0};

    INFO("Mesh >> read mesh!\n");
    read_mesh();

    INFO("Step >> Main solve start!\n\n");
	
	for(int mesh = 1; mesh <= 3; mesh++){

		double run_time = START_TIME;
		while(run_time <= STOP_TIME){
			INFO("Step >> mesh : %d, run time: %.2f\n", mesh, run_time);
			// init p equation
			LduMatrix ldu_matrix;
			double *source;
			double *psi;
			
			init_p_equation(ldu_matrix, source, psi, mesh);

			// solve
			gettimeofday(&tv,NULL);
			double pcg_start = (double)(tv.tv_sec)+(double)(tv.tv_usec)*1e-6;
			PCGReturn pcg_return = pcg_solve(ldu_matrix, source, psi, MAX_ITERS, tolerance, normfactor);
			gettimeofday(&tv,NULL);
			double pcg_end = (double)(tv.tv_sec)+(double)(tv.tv_usec)*1e-6;
			pcg_time_cost[mesh - 1] += pcg_end - pcg_start;
			INFO("Step >> run time: %.2f, time cost: %.4lfs\n", run_time, pcg_end - pcg_start);
			
			// check
			if(check_result(pcg_return, run_time, mesh)){
				INFO("Check >> check result correct!\n\n");
			} else {
				INFO("Check >> check result uncorrect!\n\n");
				break;
			}

			// free
			free_p_equation(ldu_matrix, source, psi);

			run_time += DELTA_T;
		}
		total_time_cost += pcg_time_cost[mesh - 1];
		INFO("Time: mesh : %d, pcg solve time: %.4lfs\n\n", mesh, pcg_time_cost[mesh - 1]);
	}

    INFO("Time: total time cost: %.4lfs\n\n", total_time_cost);

    return 0;
}
