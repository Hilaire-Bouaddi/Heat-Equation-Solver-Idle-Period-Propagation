#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#define INIT_T 100

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
	
void showT(double **T, int NL, int NH, int n) {
	for (int i = 0; i < NH; i++) {
		for (int j = 0; j < NL; j++) {
			printf("%.2f\t", T[n][i*NL+j]);
		}
		printf("\n");
	}
}

// this function will only print enough data to generate gifs of 24 fps
void writeToFile(char *filename, double **T, int NL, int NH, int N_ITER, double dt) {
	bool write = true;
	double count = 0;
	FILE *fp = fopen(filename, "w+");
	fprintf(fp, "%d %d %d %f\n", NL, NH, N_ITER, dt); 
	for (int n = 0; n < N_ITER; n++) {
		if (write) { 
			for (int i = 0; i < NH; i++) {
				for (int j = 0; j < NL; j++) {
					fprintf(fp, "%f \t", T[n][i*NL+j]);
				}
				fprintf(fp, "\n");
			}
			write = false;
			count = 0;
		} 
		count += dt;
		if (count > 1.0/24.0) {
			write = true;
		}
		
	}
	fclose(fp);
}

int main(int argc, char** argv) {
	double L = 2.0; // length of the domain (x)
	double H = 2.0; // height of the domain (y)

	double h = 0.01; // distance between 2 points
	int NL = L/h;
	int NH = H/h;
	int N = NL*NH; // number of points 

	double dt = 0.01; // in seconds
	int T_MAX = 15; // time when the simulation ends
	int N_ITER = T_MAX/dt;

	double k = 2E-3; // diffusion coefficient

	
	printf("---------------\n");
	printf("RUNNING THE SIMULATION WITH:\nL = %fm\nH = %fm\nh = %fm\ndt = %fs\nT_MAX = %ds\nk = %fm^2/s\n", L, H, h, dt, T_MAX, k); 
	printf("---------------\n");
	

	uint32_t start_hi=0, start_lo=0;
	uint32_t   end_hi=0,   end_lo=0;
	RDTSC_START();

	// allocation
	double **T = malloc(N_ITER*sizeof(double)); // theta
	for (int n = 0; n < N_ITER; n++) 
		T[n] = malloc(N * sizeof(double));

	// initializing the borders to INIT_T degrees
	for (int n = 0; n < N_ITER; n++) {
		for (int i=0; i < NL; i++) {
			T[n][i] = INIT_T;
	       		T[n][N - i - 1] = INIT_T;
		}	       
		for (int i=0; i < NH; i++) {
			T[n][i*NL] = INIT_T;
			T[n][(i+1)*NL - 1] = INIT_T;
		}
	}

	// check initialization 
	// showT(T, NL, NH, 0);

	for (int n = 0; n < N_ITER - 1; n++) {
		for (int j = 1; j < NH - 1; j++) {
			for (int i = 1; i < NL - 1; i++) {
				T[n+1][i + j*NL] = T[n][i + j*NL] + k * dt * ( \
				   (T[n][i+1 + j*NL] - 2*T[n][i + j*NL] + T[n][i - 1 + j*NL]) / h/h + \
				   (T[n][i + (j+1)*NL] - 2*T[n][i + j*NL] + T[n][i + (j-1)*NL]) / h/h \
				);
			}		
		}
	}
	

	RDTSC_STOP();
	uint64_t e = elapsed(start_hi, start_lo, end_hi, end_lo);	
	printf("N_ITER * NH * NL = %d\n", N_ITER * NH * NL);  
	printf("Time elapsed %f s\n", e/2.3e9);

	//showT(T, NL, NH, N_ITER - 1);
	if (argc > 1 && strcmp(argv[1], "save") == 0)
		writeToFile("results_serial.txt", T, NL, NH, N_ITER, dt);
}


