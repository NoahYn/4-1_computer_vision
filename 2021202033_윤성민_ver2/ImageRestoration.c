#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>

#pragma warning(disable:4996)

#define W 512 // width of image
#define H 512 // height of image

int cmp(const unsigned char* a, const unsigned char* b) { // compare function for quick sort
	return ((*a > *b) ? 1 : -1); 
}

void getGaussFilter(double(*filter)[5], int size, double std) { // size = window size, std = standard deviation
	double sum = 0;
	double dist = (size - 1) / 2.0; // distance from origin to edge

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double x = i - dist, y = j - dist; 
			filter[i][j] = exp(-(x * x + y * y) / (2 * std * std)); // gaussian filter formular
			sum += filter[i][j];
		}
	}
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) filter[i][j] /= sum;
}

double getGaussRand(double mean, double std) { // Marsaglia polar method
	double u, v, s;

	do {
		u = 2 * (rand() / (double)RAND_MAX) - 1; // value between -1 ~ 1
		v = 2 * (rand() / (double)RAND_MAX) - 1;
		s = u * u + v * v;
	} while (s >= 1 || s == 0);
	s = sqrt((-2 * log(s)) / s);
	return mean + std * u * s;
}

int main() {
	/* open file */
	FILE* LENA = fopen("lena(512x512).raw", "rb");
	FILE* AGN = fopen("Gaussian_noise(512x512).raw", "wb");
	FILE* SPN = fopen("Salt&pepper_noise(512x512).raw", "wb");
	FILE* GG3 = fopen("Gnoise_Gauss33(512x512).raw", "wb"); FILE* GG5 = fopen("Gnoise_Gauss55(512x512).raw", "wb");
	FILE* GM3 = fopen("Gnoise_Median33(512x512).raw", "wb"); FILE* GM5 = fopen("Gnoise_Median55(512x512).raw", "wb");
	FILE* GA23 = fopen("Gnoise_Alpha233(512x512).raw", "wb");	FILE* GA15 = fopen("Gnoise_Alpha155(512x512).raw", "wb");
	FILE* GA43 = fopen("Gnoise_Alpha433(512x512).raw", "wb");	FILE* GA45 = fopen("Gnoise_Alpha455(512x512).raw", "wb");
	FILE* SG3 = fopen("Salt_Gauss33(512x512).raw", "wb"); FILE* SG5 = fopen("Salt_Gauss55(512x512).raw", "wb");
	FILE* SM3 = fopen("Salt_Median33(512x512).raw", "wb"); FILE* SM5 = fopen("Salt_Median55(512x512).raw", "wb");
	FILE* SA23 = fopen("Salt_Alpha233(512x512).raw", "wb"); FILE* SA15 = fopen("Salt_Alpha155(512x512).raw", "wb");
	FILE* SA43 = fopen("Salt_Alpha433(512x512).raw", "wb"); FILE* SA45 = fopen("Salt_Alpha455(512x512).raw", "wb");
	srand((int)time(0));

	double PSNR;
	double SSE;
	double MSE;

	unsigned char origin[H][W];
	unsigned char gaussian[H][W];
	unsigned char saltpepper[H][W];
	double filter[5][5]; // for gaussian filter
	unsigned char filter2[25]; // for median and alpha-trimmed mean filter

	/* read image from lena.raw */
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			fscanf(LENA, "%c", &origin[i][j]);
			saltpepper[i][j] = origin[i][j];
		}
	}
	printf("PSNR TABLE\n");
	printf("\t\t\t\tGaussian noise  |  Salt & pepper noise\n");
	printf("----------------------------------------------------------------------------\n");

	/* step 1 : generating noise */

	/* write additive gaussian noise image and print PSNR */
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			char gaussnoise = getGaussRand(0, 10);
			gaussian[i][j] = origin[i][j] + gaussnoise; // add noise
			fputc(gaussian[i][j], AGN);
			SSE += gaussnoise * gaussnoise;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("Before filtering\t\t%.4f\t\t|\t\t", PSNR);


	// make salt & pepper noise
	for (int i = 0; i < W * H * 0.3; i++) {
		int randi = rand() % H;
		int randj = rand() % W;
		saltpepper[randi][randj] = (rand() % 2) * 255; // 0 or 255
	}

	// write salt & pepper image and print PSNR
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			fputc(saltpepper[i][j], SPN);
			double diff = origin[i][j] - saltpepper[i][j];
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);


	/* step 2. filtering */
	/* generate array : index of filter */
	int dx3[9];
	int dy3[9];
	for (int i = 0; i < 9; i++) {
		dx3[i] = i % 3; // 012 012 012
		dy3[i] = i / 3; // 000 111 222
	}
	int dx5[25];
	int dy5[25];
	for (int i = 0; i < 25; i++) {
		dx5[i] = i % 5; // 01234 01234 ...
		dy5[i] = i / 5; // 00000 11111 ...
	}

	/* 3x3 Gaussian filtering -> Gaussian noise */
	printf("----------------------------------------------------------------------------\n");
	printf("After Gaussian 3x3\t\t");
	SSE = 0;
	getGaussFilter(filter, 3, 0.7);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double conv = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) continue; // zero padding : nothing to add
				conv += gaussian[x][y] * filter[dx3[dir]][dy3[dir]]; // convolution 
			}
			fputc((unsigned char)conv, GG3);
			double diff = origin[i][j] - (unsigned char)conv;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 3x3 Gaussian filtering -> Salt & pepper noise */
	SSE = 0;
	getGaussFilter(filter, 3, 10);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double conv = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) continue; // zero padding : nothing to add
				conv += saltpepper[x][y] * filter[dx3[dir]][dy3[dir]];
			}
			fputc((unsigned char)conv, SG3);
			double diff = origin[i][j] - (unsigned char)conv;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);

	/* 5x5 gaussian filtering -> gaussian noise image */
	printf("----------------------------------------------------------------------------\n");
	printf("After Gaussian 5x5\t\t");
	SSE = 0;
	getGaussFilter(filter, 5, 0.7);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double conv = 0;
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) continue; // zero padding : nothing to add
				conv += gaussian[x][y] * filter[dx5[dir]][dy5[dir]];
			}
			fputc((unsigned char)conv, GG5);
			double diff = origin[i][j] - (unsigned char)conv;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 5x5 Gaussian filtering -> Salt & pepper noise */
	SSE = 0;
	getGaussFilter(filter, 5, 10);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double conv = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) continue; // zero padding : nothing to add
				conv += saltpepper[x][y] * filter[dx5[dir]][dy5[dir]];
			}
			fputc((unsigned char)conv, SG5);
			double diff = origin[i][j] - (unsigned char)conv;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);

	

	/* 3x3 Median filtering -> Gaussian noise */
	printf("----------------------------------------------------------------------------\n");
	printf("After Median 3x3\t\t");
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else filter2[dir] = gaussian[x][y];
			}
			qsort(filter2, 9, sizeof(unsigned char), cmp); // sorting
			fputc(filter2[4], GM3); // median
			double diff = origin[i][j] - filter2[4];
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 3x3 Median filtering -> Salt & pepper noise */
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else filter2[dir] = saltpepper[x][y];
			}
			qsort(filter2, 9, sizeof(unsigned char), cmp);
			fputc(filter2[4], SM3);
			double diff = origin[i][j] - filter2[4];
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);

	/* 5x5 Median filtering -> gaussian noise image */
	printf("----------------------------------------------------------------------------\n");
	printf("After Median 5x5\t\t");
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else filter2[dir] = gaussian[x][y];
			}
			qsort(filter2, 25, sizeof(unsigned char), cmp);
			fputc(filter2[12], GM5);
			double diff = origin[i][j] - filter2[12];
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 5x5 Median filtering -> Salt & pepper noise */
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else filter2[dir] = saltpepper[x][y];
			}
			qsort(filter2, 25, sizeof(unsigned char), cmp);
			fputc(filter2[12], SM5);
			double diff = origin[i][j] - filter2[12];
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);



	/* 3x3 0.2 Alpha-trimmed mean filtering -> Gaussian noise */
	printf("----------------------------------------------------------------------------\n");
	printf("After Alpha 0.2 3x3\t\t");
	SSE = 0;
	// alpha = 1 (floor of 0.2 * 3 * 3) 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else {
					filter2[dir] = gaussian[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 9, sizeof(unsigned char), cmp);
			sum -= (filter2[0] + filter2[8]); // alpha-trim
			sum /= 7; // mean
			fputc((unsigned char)sum, GA23);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 3x3 0.2 Alpha-trimmed mean filtering -> Salt & pepper noise */
	// alpha = 1 (floor of 0.2 * 3 * 3) 
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else {
					filter2[dir] = saltpepper[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 9, sizeof(unsigned char), cmp);
			sum -= (filter2[0] + filter2[8]);
			sum /= 7;
			fputc((unsigned char)sum, SA23);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);


	/* 3x3 0.4 Alpha-trimmed mean filtering -> Gaussian noise */
	printf("----------------------------------------------------------------------------\n");
	printf("After Alpha 0.4 3x3\t\t");
	SSE = 0;
	int alpha = 3; // (floor of 0.4 * 3 * 3) 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else {
					filter2[dir] = gaussian[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 9, sizeof(unsigned char), cmp);
			for (int i = 0; i < alpha; i++) {
				sum -= (filter2[i] + filter2[8 - i]);
			}
			sum /= 3;
			fputc((unsigned char)sum, GA43);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 3x3 0.4 Alpha-trimmed mean filtering -> Salt & pepper noise */
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx3[dir] - 1;
				int y = j + dy3[dir] - 1;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding
				else {
					filter2[dir] = saltpepper[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 9, sizeof(unsigned char), cmp);
			for (int i = 0; i < alpha; i++) {
				sum -= (filter2[i] + filter2[8 - i]);
			}
			sum /= 3;
			fputc((unsigned char)sum, SA43);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);

	/* 5x5 0.1 Alpha-trimmed mean filtering -> gaussian noise image */
	printf("----------------------------------------------------------------------------\n");
	printf("After Alpha 0.1 5x5\t\t");
	SSE = 0;
	alpha = 2; // 5 * 5 * 0.1
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding 
				else {
					filter2[dir] = gaussian[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 25, sizeof(unsigned char), cmp);
			for (int i = 0; i < alpha; i++) {
				sum -= (filter2[i] + filter2[24 - i]);
			}
			sum /= 21;
			fputc((unsigned char)sum, GA15);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 5x5 0.1 Alpha-trimmed mean filtering -> Salt & pepper noise */
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding 
				else {
					filter2[dir] = saltpepper[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 25, sizeof(unsigned char), cmp);
			for (int i = 0; i < alpha; i++) {
				sum -= (filter2[i] + filter2[24 - i]);
			}
			sum /= 21;
			fputc((unsigned char)sum, SA15);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);

	/* 5x5 0.4 Alpha-trimmed mean filtering -> gaussian noise image */
	printf("----------------------------------------------------------------------------\n");
	printf("After Alpha 0.4 5x5\t\t");
	SSE = 0;
	alpha = 10; // 5 * 5 * 0.4
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding 
				else {
					filter2[dir] = gaussian[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 25, sizeof(unsigned char), cmp);
			for (int i = 0; i < alpha; i++) {
				sum -= (filter2[i] + filter2[24 - i]);
			}
			sum /= 5;
			fputc((unsigned char)sum, GA45);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\t\t|\t\t", PSNR);

	/* 5x5 0.4 Alpha-trimmed mean filtering -> Salt & pepper noise */
	SSE = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			float sum = 0;
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx5[dir] - 2;
				int y = j + dy5[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) filter2[dir] = 0; // zero padding 
				else {
					filter2[dir] = saltpepper[x][y];
					sum += filter2[dir];
				}
			}
			qsort(filter2, 25, sizeof(unsigned char), cmp);
			for (int i = 0; i < alpha; i++) {
				sum -= (filter2[i] + filter2[24 - i]);
			}
			sum /= 5;
			fputc((unsigned char)sum, SA45);
			double diff = origin[i][j] - sum;
			SSE += diff * diff;
		}
	}
	MSE = SSE / (W * H);
	PSNR = 10 * log10(255 * 255 / MSE);
	printf("%.4f\n", PSNR);

	fclose(LENA);
	fclose(AGN); fclose(SPN);
	fclose(GG3); fclose(GG5); fclose(GM3); fclose(GM5); fclose(GA23); fclose(GA43); fclose(GA15); fclose(GA45);
	fclose(SG3); fclose(SG5); fclose(SM3); fclose(SM5); fclose(SA23); fclose(SA43); fclose(SA15); fclose(SA45);
}