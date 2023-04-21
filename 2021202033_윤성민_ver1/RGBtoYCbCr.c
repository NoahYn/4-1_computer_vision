#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>

#define N 288
#define M 352

int main() {
	FILE* raw_rgb; // original source
	FILE* raw_ycbcr444, * raw_y444, * raw_cb444, * raw_cr444; // 444 format
	FILE* raw_ycbcr420, * raw_y420, * raw_cb420, * raw_cr420; // 420 format
	unsigned char* r, * g, * b; // buffer to contain
	unsigned char* y, * cb, * cr; // buffer to ycbcr
	//	raw_rgb = fopen("Suzie_CIF_352x288.raw", "rb"); 
	//	raw_ycbcr444 = fopen("Ycbcr444_352x288.raw", "wb"); 
	//	raw_y444 = fopen("Y444_352x288.raw", "wb"); 
	//	raw_cr444 = fopen("cr444_352x288.raw", "wb"); 
	//	raw_cb444 = fopen("cb444_352x288.raw", "wb"); 
	//	raw_ycbcr420 = fopen("Ycbcr420_352x288.raw", "wb"); 
	//	raw_y420 = fopen("Y420_352x288.raw", "wb");
	//	raw_cr420 = fopen("cr420_352x288.raw", "wb"); 
	//	raw_cb420 = fopen("cb420_352x288.raw", "wb"); 
	fopen_s(&raw_rgb, "Suzie_CIF_352x288.raw", "rb"); // open source
	fopen_s(&raw_ycbcr444, "Ycbcr444_352x288.raw", "wb"); // open ycbcr444 
	fopen_s(&raw_y444, "Y444_352x288.raw", "wb"); // open y444 
	fopen_s(&raw_cb444, "cb444_352x288.raw", "wb"); // open cb444 
	fopen_s(&raw_cr444, "cr444_352x288.raw", "wb"); // open cr444 
	fopen_s(&raw_ycbcr420, "Ycbcr420_352x288.raw", "wb"); // open ycbcr420
	fopen_s(&raw_y420, "Y420_352x288.raw", "wb"); // open y444 
	fopen_s(&raw_cb420, "cb420_352x288.raw", "wb"); // open cb444 
	fopen_s(&raw_cr420, "cr420_352x288.raw", "wb"); // open cr444 


	/* allocate memory */
	r = (unsigned char*)malloc(sizeof(unsigned char) * N * M);
	g = (unsigned char*)malloc(sizeof(unsigned char) * N * M);
	b = (unsigned char*)malloc(sizeof(unsigned char) * N * M);
	y = (unsigned char*)malloc(sizeof(unsigned char) * N * M);
	cb = (unsigned char*)malloc(sizeof(unsigned char) * N * M);
	cr = (unsigned char*)malloc(sizeof(unsigned char) * N * M);

	/* read byte from RGB */
	for (int i = 0; i < N * M; i++)
		// fscanf(raw_rgb, "%c", &b[i]);
		fscanf_s(raw_rgb, "%c", &b[i]);
	for (int i = 0; i < N * M; i++)
		// fscanf(raw_rgb, "%c", &g[i]);
		fscanf_s(raw_rgb, "%c", &g[i]);
	for (int i = 0; i < N * M; i++)
		// fscanf(raw_rgb, "%c", &r[i]);
		fscanf_s(raw_rgb, "%c", &r[i]);

	/* get ycbcr from rgb by given formular */
	for (int i = 0; i < N * M; i++) {
		y[i] = (unsigned char)((double)77 / 256 * r[i] + (double)150 / 256 * g[i] + (double)29 / 256 * b[i]);
		cr[i] = (unsigned char)(128 + (double)131 / 256 * r[i] - (double)110 / 256 * g[i] - (double)21 / 256 * b[i]);
		cb[i] = (unsigned char)(128 - (double)44 / 256 * r[i] - (double)87 / 256 * g[i] + (double)131 / 256 * b[i]);
	}

	// ycbcr 444 format
	for (int i = 0; i < N * M; i++)
		fputc(y[i], raw_ycbcr444);
	for (int i = 0; i < N * M; i++)
		fputc(cr[i], raw_ycbcr444);
	for (int i = 0; i < N * M; i++)
		fputc(cb[i], raw_ycbcr444);

	// y 444 
	for (int i = 0; i < N * M; i++)
		fputc(y[i], raw_y444);

	// cr 444
	for (int i = 0; i < N * M; i++)
		fputc(cr[i], raw_cr444);

	// cb 444
	for (int i = 0; i < N * M; i++)
		fputc(cb[i], raw_cb444);

	// ycbcr 420 
	for (int i = 0; i < N * M; i++)
		fputc(y[i], raw_ycbcr420);
	for (int i = 0; i < N; i += 2) {
		for (int j = 0; j < M; j += 2) {
			fputc(cr[i * M + j], raw_ycbcr420);
		}
	}
	for (int i = 0; i < N; i += 2) {
		for (int j = 0; j < M; j += 2) {
			fputc(cb[i * M + j], raw_ycbcr420);
		}
	}

	// y 420
	for (int i = 0; i < N * M; i++)
		fputc(y[i], raw_y420);

	// cb 420
	for (int i = 0; i < N; i += 2) {
		for (int j = 0; j < M; j += 2) {
			fputc(cb[i * M + j], raw_cb420);
		}
	}

	// cr 420
	for (int i = 0; i < N; i += 2) {
		for (int j = 0; j < M; j += 2) {
			fputc(cr[i * M + j], raw_cr420);
		}
	}

	fclose(raw_rgb);
	fclose(raw_ycbcr444);
	fclose(raw_y444);
	fclose(raw_cr444);
	fclose(raw_cb444);
	fclose(raw_ycbcr420);
	fclose(raw_y420);
	fclose(raw_cr420);
	fclose(raw_cb420);
}