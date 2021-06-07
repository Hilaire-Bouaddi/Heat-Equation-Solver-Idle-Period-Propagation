#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define PROCESSOR_CLOCK_SPEED 2.3E9
#define INIT_T 100

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

// Use the preprocessor so we know definitively that these are placed inline
#define RDTSC_START()            \
	__asm__ volatile("CPUID\n\t" \
	                 "RDTSC\n\t" \
       	                 "mov %%edx, %0\n\t" \
      	                 "mov %%eax, %1\n\t" \
    	                 : "=r" (start_hi), "=r" (start_lo) \
      	                 :: "%rax", "%rbx", "%rcx", "%rdx");

#define RDTSC_STOP()              \
       	__asm__ volatile("RDTSCP\n\t" \
	                 "mov %%edx, %0\n\t" \
      	                 "mov %%eax, %1\n\t" \
      	                 "CPUID\n\t" \
      	                 : "=r" (end_hi), "=r" (end_lo) \
       	                 :: "%rax", "%rbx", "%rcx", "%rdx");
// Returns the elapsed time given the high and low bits of the start and stop time.
uint64_t elapsed(uint32_t start_hi, uint32_t start_lo, uint32_t end_hi,   uint32_t end_lo)
{
       	uint64_t start = (((uint64_t)start_hi) << 32) | start_lo;
	uint64_t end   = (((uint64_t)end_hi)   << 32) | end_lo;
	return end-start;
}

// shows the temperature at the timestep n 
void showT(double **T, int NL, int NH, int n, int *coords) {
	printf("----------------------------\n%d %d\n", coords[0], coords[1]);
	for (int i = 0; i < NH; i++) {
		for (int j = 0; j < NL; j++) {
			printf("%.2f\t", T[n][i*NL+j]);
		}
		printf("\n");
	}
}

// this function will only print enough data to generate gifs of 24 fps
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

// this function checks the plausibility of the results: is the temperature increasing and always lower than the thermostat
bool isPlausible(double **T, int N_ITER, int N) {
	for (int n = 0; n < N_ITER; n++) {
		for (int i = 0; i < N; i++) {
			if (T[n][i] > INIT_T) 
				return false;
	
			if (n > 1) {
				if (T[n][i] < T[n-1][i])
					return false;
			}
		}
	}

	return true;
}

int main(int argc, char **argv) {
	if (argc < 7) {
		printf("Usage: parallel <L> <H> <h> <dt> <T_MAX> <k> [save_results] [save_idle_times]\n");
		return 0;
	}
	
	double L = strtod(argv[1], NULL); // length of the domain (x)
	double H = strtod(argv[2], NULL); // height of the domain (y)

	double h = strtod(argv[3], NULL); // distance between 2 points
	int NL = L/h;
	int NH = H/h;
	int N = NL*NH; // number of points 

	double dt = strtod(argv[4], NULL); // in seconds
	int T_MAX = strtol(argv[5], NULL, 10); // time when the simulation ends
	int N_ITER = T_MAX/dt;
	double k = strtod(argv[6], NULL); // diffusion coefficient

	bool save_idle_times = false, save_results = false;

	// decide if we are going to save the results
	for (int i = 7; i < argc; i++) {
		if (strcmp(argv[i], "save_idle_times") == 0) {
			save_idle_times = true;
		}
		if (strcmp(argv[i], "save_results") == 0) {
			save_results = true;
		}
	}
	
	uint64_t *idle_times = malloc(sizeof(uint64_t) * (N_ITER - 1));

	// setting up MPI	
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	
	int rank_world, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);


	if (rank_world == 0) {
		printf("---------------\n");
		printf("RUNNING THE SIMULATION WITH:\nL = %fm\nH = %fm\nh = %fm\ndt = %fs\nT_MAX = %ds\nk = %fm^2/s\n", L, H, h, dt, T_MAX, k); 
		printf("N = %d\n", N);
		printf("final size of matrix T = %.3fMB\n", sizeof(double)*N_ITER*N/1000000.0);
		if (dt < 1/24.0) printf("size of matrix T with dt = 1/24s = %.3fMB\n", sizeof(double)*N_ITER*dt*24*N/1E6);
		printf("---------------\n");
		fflush(stdout);
	}
	
	if (size > NH || size > NL) {
		if (rank_world == 0) printf("grid too small for parallelization\n");
		return 0;
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

	for (int n = 0; n < N_ITER; n++) {
		if (coords[0] == 0) {
			for (int j = 0; j < mH; j++) {
				T[n][j*mL] = INIT_T;
			}
		}
		if (coords[0] == dims[0]-1) {
			for (int j = 0; j < mH; j++) {
				T[n][j*mL + mL - 1] = INIT_T;
			}
		}
	
		if (coords[1] == 0) {
			for (int i = 0; i < mL; i++) {
				T[n][i] = INIT_T;
			}
		}
		if (coords[1] == dims[1]-1) {
			for (int i = 0; i < mL; i++) {
				T[n][(mH - 1) * mL + i] = INIT_T;
			}
		}
	}
	//showT(T, mL, mH, 0, coords);

	double *T_from_left = malloc(sizeof(double) * mH);
	double *T_to_left = malloc(sizeof(double) * mH);
	double *T_from_right = malloc(sizeof(double) * mH);
	double *T_to_right = malloc(sizeof(double) * mH);
	double *T_from_top = malloc(sizeof(double) * mL);
	double *T_to_top = malloc(sizeof(double) * mL);
	double *T_from_bottom = malloc(sizeof(double) * mL);
	double *T_to_bottom = malloc(sizeof(double) * mL);


	uint32_t start_hi=0, start_lo=0;
	uint32_t   end_hi=0,   end_lo=0;

	RDTSC_START();

	uint32_t start_hi_copy=start_hi, start_lo_copy=start_lo;

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
			int rank_top;
			MPI_Cart_rank(cart_comm, coords_dest, &rank_top);
			MPI_Isend(T_to_top, mL, MPI_DOUBLE, rank_top, 0, cart_comm, &requests[request_counter++]); 
			MPI_Irecv(T_from_top, mL, MPI_DOUBLE, rank_top, 0, cart_comm, &requests[request_counter++]);
		}

		if (coords[1] < dims[1] - 1) { // if we are not at the very bottom
			int coords_dest[2] = {coords[0], coords[1] + 1};
			int rank_bottom;
			MPI_Cart_rank(cart_comm, coords_dest, &rank_bottom);
			MPI_Isend(T_to_bottom, mL, MPI_DOUBLE, rank_bottom, 0, cart_comm, &requests[request_counter++]); 
			MPI_Irecv(T_from_bottom, mL, MPI_DOUBLE, rank_bottom, 0, cart_comm, &requests[request_counter++]);
		}
		
		// measure idle time
		
		if (save_idle_times) {
			start_hi=0, start_lo=0;
			end_hi=0,   end_lo=0;
			RDTSC_START();
			MPI_Waitall(request_counter, requests, statuses);
			RDTSC_STOP();
			uint64_t e = elapsed(start_hi, start_lo, end_hi, end_lo);
			idle_times[n] = e;
		}
	
		// COMPUTATION
		// inside of the tile
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
			// top right corner 
			if (coords[1] > 0) {
				T[n+1][mL-1] = T[n][mL-1] + k * dt * ( \
					 (T_from_right[0] - 2*T[n][mL-1] + T[n][mL-1 - 1]) / h/h + \
					 (T[n][mL-1 + 1*mL] - 2*T[n][mL-1] + T_from_top[mL-1]) / h/h \
					);
			}

			// bottom right corner 
			if (coords[1] < dims[1] - 1) {
				T[n+1][mL-1 + (mH-1)*mL] = T[n][mL-1 + (mH-1)*mL] + k * dt * ( \
					 (T_from_right[mH-1] - 2*T[n][mL-1 + (mH-1)*mL] + T[n][mL-1 - 1 + (mH-1)*mL]) / h/h + \
					 (T_from_bottom[(mH-1)] - 2*T[n][mL-1 + (mH-1)*mL] + T[n][mL-1 + ((mH-1)-1)*mL]) / h/h \
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

			// top left corner 
			if (coords[1] > 0) {
				T[n+1][0] = T[n][0] + k * dt * ( \
						 (T[n][1] - 2*T[n][0] + T_from_left[0]) / h/h + \
						 (T[n][mL] - 2*T[n][0] + T_from_top[0]) / h/h \
						);
			}

			// bottom left corner 
			if (coords[1] < dims[1] - 1) {
				T[n+1][(mH-1)*mL] = T[n][(mH-1)*mL] + k * dt * ( \
						 (T[n][1 + (mH-1)*mL] - 2*T[n][(mH-1)*mL] + T_from_left[(mH-1)]) / h/h + \
						 (T_from_bottom[0] - 2*T[n][(mH-1)*mL] + T[n][((mH-1)-1)*mL]) / h/h \
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
	}

	end_hi=0;
     	end_lo=0;
	start_hi = start_hi_copy;
	start_lo = start_lo_copy;

	if (rank == 0) { 
		RDTSC_STOP();
		uint64_t elapsed_time = elapsed(start_hi, start_lo, end_hi, end_lo);
		printf("time elapsed: %lf s\n", elapsed_time / PROCESSOR_CLOCK_SPEED ); 
		fflush(stdout);
	}
	// showT(T, mL, mH, N_ITER-1, coords);	

	// saving idle times
	if (save_idle_times) {
		struct stat st = {0};

		if (stat("/idle_times", &st) == -1) {
	    		mkdir("idle_times", 0700);
		}

		char str_coords[10];
		sprintf(str_coords, "%d,%d", coords[0], coords[1]);
	
		char *filename = malloc(sizeof(char)*30);
		filename = strcpy(filename, "idle_times/coords_");
		strcat(strcat(filename, str_coords), ".txt");
		FILE *fp = fopen(filename, "w+");
		fprintf(fp, "%d %d %d %f\n", NL, NH, N_ITER, dt); 
		for (int n = 0; n < N_ITER - 1; n++) {
			fprintf(fp, "%ld ", idle_times[n]);	
		}
		fclose(fp);
	}

	if (!isPlausible(T, N_ITER, mH * mL)) 
		printf("The simulation looks unstable. Consider lowering timestep\n");

	// merging to be able to write to a file 
	if (!save_results) return 0;	
	if (coords[0] == 0 && coords[1] == 0) {
		double **final_T = malloc(sizeof(double*) * N_ITER);	
		for (int n = 0; n < N_ITER; n++) {
			final_T[n] = malloc(sizeof(double) * N);
			
			for (int i = 0; i < mH; i++) {
				for (int j = 0; j < mL; j++) {
					final_T[n][i * NL + j] = T[n][i * mL + j]; 
				}	
			}

			int m_recv[2];
			MPI_Status status;
			for (int y = 0; y < dims[1]; y++) {
				for (int x = 0; x < dims[0]; x++) {
					if (x != 0 || y!=0) {
						int coords_recv[] = {x, y};
						int rank_recv;
						MPI_Cart_rank(cart_comm, coords_recv, &rank_recv);
					
						MPI_Recv(m_recv, 2, MPI_INT, rank_recv, 0, cart_comm, &status);

						double *T_tile = malloc(sizeof(double) * m_recv[0] * m_recv[1]);
						
						MPI_Recv(T_tile, m_recv[0] * m_recv[1], MPI_DOUBLE, rank_recv, 0, cart_comm, &status);

						int start = coords_recv[0] * floor(NL/floor(sqrt(size))) + coords_recv[1]* floor(NH/floor(sqrt(size))) * NL;
						for (int i = 0; i < m_recv[1]; i++) {
							for (int j = 0; j < m_recv[0]; j++) {
								final_T[n][start + i * NL + j] = T_tile[i * m_recv[0] + j]; 
							}	
						}
					}
				}	
			}
		}

		//showT(final_T, NL, NH, N_ITER-1, coords);	
		if (save_results)
			writeToFile("results_parallel.txt", final_T, NL, NH, N_ITER, dt);
	} 
	else {
		int coords_send[] = {0, 0};
		int rank_send;
		MPI_Cart_rank(cart_comm, coords_send, &rank_send);
		int m[] = {mL, mH};
		for (int n = 0; n < N_ITER; n++) {
			MPI_Send(m, 2, MPI_INT, rank_send, 0, cart_comm);
			MPI_Send(T[n], mL*mH, MPI_DOUBLE, rank_send, 0, cart_comm);
		}
	}

}


