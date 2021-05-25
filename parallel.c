#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define INIT_T 100

// shows the temperature at the timestep n 
void showT(double **T, int NL, int NH, int n, int *coords) {
	printf("%d %d\n", coords[0], coords[1]);
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

	double h = 0.1; // distance between 2 points
	int NL = L/h;
	int NH = H/h;
	int N = NL*NH; // number of points 

	double dt = 0.01; // in seconds
	int T_MAX = 5; // time when the simulation ends
	int N_ITER = T_MAX/dt;

	double k = 2E-3; // diffusion coefficient


	// setting up MPI	
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	int rank, size;
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
	
	MPI_Comm_rank(cart_comm, &rank);
	int coords[2];
	MPI_Cart_coords(cart_comm, rank, 2, coords);
	
	double **T = malloc(N_ITER * sizeof(double*));
	int mL = floor(NL/floor(sqrt(size))); // the dimension of the tile on the x axis
	int mH = floor(NH/floor(sqrt(size))); // the dimension of the tile on the y axis

	// check if there are leftovers on the right
	if (coords[0] == dims[0] - 1) {
		mL = NL - coords[0] * mL;
	}
	// check same in the bottom
	if (coords[1] == dims[1] - 1) {
		mH = NH - coords[1] * mH;
	}
	// checking that all the tiles are the right size
	//printf("(%d;%d): %d, %d\n", coords[0], coords[1], mL, mH);
	
	// initalizing the borders to INIT_T degrees
	T[0] = malloc(sizeof(double) * mL * mH);
	for (int i = 0; i < mL*mH; i++) 
		T[0][i] = 0;

	if (coords[0] == 0) {
		for (int j = 0; j < mH; j++) {
			T[0][j*mL] = INIT_T;
		}
	}
	if (coords[0] == dims[0]-1) {
		for (int j = 0; j < mH; j++) {
			T[0][(j + 1) * mL - 1] = INIT_T;
		}
	}

	if (coords[1] == 0) {
		for (int i = 0; i < mL; i++) {
			T[0][i] = INIT_T;
		}
	}
	if (coords[1] == dims[1]-1) {
		for (int i = 0; i < mL; i++) {
			T[0][(mH - 1) * mL + i] = INIT_T;
		}
	}
	//showT(T, mL, mH, 0, coords);

	for (int n = 0; n < N_ITER-1; n++) {
		// halo to get the borders from neighboring threads 
		// let's start by shifting horizontally
		int rank_dest;
		int coords_dest[2];
		double *T_from_left = malloc(sizeof(double) * mH);
		MPI_Cart_shift(cart_comm, 0, 1, &rank, &rank_dest);
		MPI_Cart_coords(cart_comm, rank_dest, 2, coords_dest);
		
		   int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
			                          int dest, int sendtag,
				                              void *recvbuf, int recvcount, MPI_Datatype recvtype,
					                           int source, int recvtag,
					                            MPI_Comm comm, MPI_Status *status)			   
		MPI_Sendrecv();
		
		//printf("(%d, %d) wants to send to (%d, %d)\n", coords[0], coords[1], coords_dest[0], coords_dest[1]);
		break;
		// computing	
	}
	// writeToFile("results_parallel.txt", T, NL, NH, N_ITER, dt);
}


