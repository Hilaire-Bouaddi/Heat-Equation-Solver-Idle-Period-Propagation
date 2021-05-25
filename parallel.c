#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define INIT_T 100

void showT(double **T, int NL, int NH, int n) {
	for (int i = 0; i < NH; i++) {
		for (int j = 0; j < NL; j++) {
			printf("%.2f\t", T[n][i*NL+j]);
		}
		printf("\n");
	}
}

void writeToFile(char *filename, double **T, int NL, int NH, int N_ITER, double dt) {
	FILE *fp = fopen(filename, "w+");
	fprintf(fp, "%d %d %d %f\n", NL, NH, N_ITER, dt); 
	for (int n = 0; n < N_ITER; n++) {
		for (int i = 0; i < NH; i++) {
			for (int j = 0; j < NL; j++) {
				fprintf(fp, "%f \t", T[n][i*NL+j]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
}

int main(int argc, char **argv) {
	double L = 1.0; // length of the domain (x)
	double H = 1.0; // height of the domain (y)

	double h = 0.01; // distance between 2 points
	int NL = L/h;
	int NH = H/h;
	int N = NL*NH; // number of points 

	double dt = 0.01; // in seconds
	int T_MAX = 5; // time when the simulation ends
	int N_ITER = T_MAX/dt;

	double k = 2E-3; // diffusion coefficient


	// setting up MPI	
	int provided, rank, size;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		printf("---------------\n");
		printf("RUNNING THE SIMULATION WITH:\nL = %fm\nH = %fm\nh = %fm\ndt = %fs\nT_MAX = %ds\nk = %fm^2/s\n", L, H, h, dt, T_MAX, k); 
		printf("---------------\n");
	}

	// let's try to see how we can divide the domain 
	// we will assume that NL == NH == sqrt(N) 

 	// creating the cartesian communicator
	MPI_Comm cart_comm;
	int dims[] = {floor(sqrt(size)), floor(sqrt(size))};
	int periods[] = {0, 0};
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
	

	// showT(T, NL, NH, N_ITER - 1);
	// writeToFile("results_parallel.txt", T, NL, NH, N_ITER, dt);
}


