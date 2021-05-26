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
	double L = 5.0; // length of the domain (x)
	double H = 5.0; // height of the domain (y)

	double h = 0.001; // distance between 2 points
	int NL = L/h;
	int NH = H/h;
	int N = NL*NH; // number of points 

	double dt = 0.1; // in seconds
	int T_MAX = 5; // time when the simulation ends
	int N_ITER = T_MAX/dt;

	double k = 2E-3; // diffusion coefficient


	// setting up MPI	
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	int rank_world, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);

	if (rank_world == 0) {
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
	
	int rank;
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
	for (int n = 0; n < N_ITER; n++) 
		T[n] = malloc(sizeof(double) * mL * mH);
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
	// showT(T, mL, mH, 0, coords);

	double *T_from_left = malloc(sizeof(double) * mH);
	double *T_to_left = malloc(sizeof(double) * mH);
	double *T_from_right = malloc(sizeof(double) * mH);
	double *T_to_right = malloc(sizeof(double) * mH);
	double *T_from_top = malloc(sizeof(double) * mL);
	double *T_to_top = malloc(sizeof(double) * mL);
	double *T_from_bottom = malloc(sizeof(double) * mL);
	double *T_to_bottom = malloc(sizeof(double) * mL);

	for (int n = 0; n < N_ITER-1; n++) {
		// halo to get the borders from neighboring threads 
		
		for (int i = 0; i < mL; i++) {
			T_to_top[i] = T[n][i];
			T_to_bottom[i] = T[n][mL * (mH - 1) + i];
		}
		for (int j = 0; j < mH; j++) {
			T_to_left[j] = T[n][j * mL];
			T_to_right[j] = T[n][(j+1) * mL - 1];
		}
		
		
		// let's start by shifting horizontally
		
		MPI_Request *requests = malloc(sizeof(MPI_Request) * 8);
		MPI_Status *statuses = malloc(sizeof(MPI_Status) * 8);

		int request_counter = 0;

		// sending to the right
		
		// if we are not at the very right, we send to the right and receive from the right
		if (coords[0] < dims[0] - 1) {
			// printf("Here are my coords (%d %d) and I want to send to the right to (%d %d)\n", coords[0], coords[1], coords[0]+1, coords[1]);
			int coords_dest[2] = {coords[0] + 1, coords[1]};
			int rank_right;
			MPI_Cart_rank(cart_comm, coords_dest, &rank_right);
			// int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
			// 		                     MPI_Comm comm, MPI_Request *request)
			MPI_Isend(T_to_right, mH, MPI_DOUBLE, rank_right, 0, cart_comm, &requests[request_counter++]); 
			//  int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
			//                       int tag, MPI_Comm comm, MPI_Request *request)
			MPI_Irecv(T_from_right, mH, MPI_DOUBLE, rank_right, 0, cart_comm, &requests[request_counter++]);
		}
		
		if (coords[0] > 0) { // if we are not on the very left
			int coords_dest[2] = {coords[0] - 1, coords[1]};
			int rank_left;
			MPI_Cart_rank(cart_comm, coords_dest, &rank_left);
			MPI_Isend(T_to_left, mH, MPI_DOUBLE, rank_left, 0, cart_comm, &requests[request_counter++]); 
			MPI_Irecv(T_from_left, mH, MPI_DOUBLE, rank_left, 0, cart_comm, &requests[request_counter++]);
		}	

		if (coords[1] > 0) { // if we are not on the very top
			int coords_dest[2] = {coords[0], coords[1] - 1};
			int rank_bottom;
			MPI_Cart_rank(cart_comm, coords_dest, &rank_bottom);
			MPI_Isend(T_to_bottom, mL, MPI_DOUBLE, rank_bottom, 0, cart_comm, &requests[request_counter++]); 
			MPI_Irecv(T_from_bottom, mL, MPI_DOUBLE, rank_bottom, 0, cart_comm, &requests[request_counter++]);
		}

		if (coords[1] < dims[1] - 1) { // if we are not at the very bottom
			int coords_dest[2] = {coords[0], coords[1] + 1};
			int rank_top;
			MPI_Cart_rank(cart_comm, coords_dest, &rank_top);
			MPI_Isend(T_to_top, mL, MPI_DOUBLE, rank_top, 0, cart_comm, &requests[request_counter++]); 
			MPI_Irecv(T_from_top, mL, MPI_DOUBLE, rank_top, 0, cart_comm, &requests[request_counter++]);
		}

		MPI_Waitall(request_counter, requests, statuses);
	
		// COMPUTATION
		double T_cell_left, T_cell_top, T_cell_right, T_cell_bottom;
		for (int j = 1; j < mH - 1; j++) {
			for (int i = 1; i < mL - 1; i++) {		
				T[n+1][i + j*mL] = T[n][i + j*mL] + k * dt * ( \
						 (T[n][i+1 + j*mL] - 2*T[n][i + j*mL] + T[n][i - 1 + j*mL]) / h/h + \
						 (T[n][i + (j+1)*mL] - 2*T[n][i + j*mL] + T[n][i + (j-1)*mL]) / h/h \
						);
			}
		}
		// right border
		if (coords[0] < dims[0] - 1) {	// Here we are trying to update the right border of the tile using T_from_right
			for (int j = 1; j < mH - 1; j++) {
				T[n+1][mL-1 + j*mL] = T[n][mL-1 + j*mL] + k * dt * ( \
					 (T_from_right[j] - 2*T[n][mL-1 + j*mL] + T[n][mL-1 - 1 + j*mL]) / h/h + \
					 (T[n][mL-1 + (j+1)*mL] - 2*T[n][mL-1 + j*mL] + T[n][mL-1 + (j-1)*mL]) / h/h \
					);
			}
		}
		
		if (coords[0] > 0) { // if we are not on the very left
			for (int j = 1; j < mH - 1; j++) {		
				T[n+1][j*mL] = T[n][j*mL] + k * dt * ( \
						 (T[n][1 + j*mL] - 2*T[n][j*mL] + T_from_left[j]) / h/h + \
						 (T[n][(j+1)*mL] - 2*T[n][j*mL] + T[n][(j-1)*mL]) / h/h \
						);
			}
		}	

		if (coords[1] > 0) { // if we are not on the very top: we set j to 0
			for (int i = 1; i < mL - 1; i++) {
				T[n+1][i] = T[n][i] + k * dt * ( \
						(T[n][i+1] - 2*T[n][i] + T[n][i - 1])/h/h + \
						(T[n][i + mL] - 2*T[n][i] + T_from_top[i])/h/h \
					    );
			}
		}

		if (coords[1] < dims[1] - 1) { // if we are not at the very bottom: we set j to mH - 1
			for (int i = 1; i < mL - 1; i++) {
				T[n+1][i + (mH-1)*mL] = T[n][i + (mH-1)* mL] + k * dt * ( \
						(T[n][i+1 + (mH-1)*mL] - 2*T[n][i + (mH-1)*mL] + T[n][i - 1 + (mH-1)*mL])/h/h + \
						(T_from_bottom[i] - 2*T[n][i + (mH-1)*mL] + T[n][i + ((mH-1)-1)*mL])/h/h \
					    );
			}
		}
		//break;
		// computing	
	}
	// writeToFile("results_parallel.txt", T, NL, NH, N_ITER, dt);
}


