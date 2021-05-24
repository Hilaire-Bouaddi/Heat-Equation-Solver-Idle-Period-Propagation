#include <stdio.h>
#include <stdlib.h>

#define INIT_T 100

void showT(double **T, int NL, int NH, int n) {
	for (int i = 0; i < NH; i++) {
		for (int j = 0; j < NL; j++) {
			printf("%.2f\t", T[n][i*NL+j]);
		}
		printf("\n");
	}
}

void writeToFile(char *filename, double **T, int NL, int NH, int N_ITER) {
	FILE *fp = fopen(filename, "w+");
	fprintf(fp, "%d %d %d\n", NL, NH, N_ITER); 
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

int main() {
	double L = 1.0; // length of the domain (x)
	double H = 1.0; // height of the domain (y)

	double h = 0.1; // distance between 2 points
	int NL = L/h;
	int NH = H/h;
	int N = NL*NH; // number of points 

	double dt = 0.1; // in seconds
	int T_MAX = 20; // time when the simulation ends
	int N_ITER = T_MAX/dt;

	double k = 1E-3; // diffusion coefficient

	printf("---------------\n");
	printf("RUNNING THE SIMULATION WITH:\nL = %fm\nH = %fm\nh = %fm\ndt = %fs\nT_MAX = %ds\nk = %fm^2/s\n", L, H, h, dt, T_MAX, k); 
	printf("---------------\n");
	
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
	

	// showT(T, NL, NH, N_ITER - 1);
	writeToFile("results_serial.txt", T, NL, NH, N_ITER);
}


