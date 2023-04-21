#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#pragma warning(disable:4996)

#define N 512
#define M 512

double getGaussRand(double m, double std) { // Marsaglia polar method
	double u, v, s;

	do {
		u = 2 * (rand() / (double)RAND_MAX) - 1; // -1 ~ 1ÀÇ °ª
		v = 2 * (rand() / (double)RAND_MAX) - 1;
		s = u * u + v * v;
	} while (s >= 1 || s == 0);
	s = sqrt((-2 * log(s)) / s);
	return m + std * u * s;
}

int main() {
	FILE* LENA = fopen("lena(512x512).raw", "rb");
	FILE* AGN = fopen("Gaussian_noise(512x512).raw", "wb");
	FILE* SPN = fopen("Salt&pepper_noise(512x512).raw", "wb");;
	srand((int)time(0));

	/* allocate memory */
	unsigned char* buffer = (unsigned char*)malloc(sizeof(unsigned) * N * M * 3);
	unsigned char* saltpepper = (unsigned char*)malloc(sizeof(unsigned) * N * M * 3);

	/* read byte from lena.raw */
	for (int i = 0; i < N * M; i++) {
		fscanf(LENA, "%c", &buffer[i]);
		saltpepper[i] = buffer[i];
	}

	// write additive gaussian noise image
	for (int i = 0; i < N * M; i++) {
		unsigned char gaussnoise = getGaussRand(0, 10);
		fputc(buffer[i] + gaussnoise, AGN);
	}

	// make salt & pepper noise
	for (int i = 0; i < N * M * 0.3; i++) {
		int randidx = rand() * rand() % (N * M); // RAND_MAX is smaller than N*M
		saltpepper[randidx] = (rand() % 2) * 255;
	}

	// write salt & pepper image
	for (int i = 0; i < N * M; i++) {
		fputc(saltpepper[i], SPN);
	}
	
	fclose(LENA);
	fclose(AGN);
	fclose(SPN);
}