#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stack>
#include <queue>
#include <tuple>
#include <vector>
#include <corecrt_math_defines.h>
#include <list>

using namespace std;
#pragma warning(disable:4996)

/*
*/
#define H 340
#define W 340
#define PATH "ball_340x340.raw"
#define HTH 50
#define CLTH 30
#define CHTH 70
#define GAUS 1.2

/*
#define H 192
#define W 256
#define PATH "insulator_256x192.raw"
#define HTH 50
#define CLTH 30
#define CHTH 70
#define GAUS 0.3
*/


/*
#define H 400
#define W 396
#define PATH "Test_img_CV_HW4_396x400.yuv"
#define HTH 60
#define CLTH 30
#define CHTH 70
#define GAUS 1.2
*/



class Canny {
public:
	unsigned int in[H][W]; // input image
	unsigned int out[H][W]; // output image(result)
	double gauss[H][W]; // after gaussian filter
	double gradient[H][W]; // magnitude of gradient 
	double theta[H][W]; // direction of gradient(in angle)
	int lowTH = CLTH, highTH = CHTH; // 2 thresholds
	stack<pair<int, int>> weak; // weak edges
}; 
Canny c;

class Hough {
public:
	unsigned int o[H][W]; // hough space for origin coordinate
	vector<pair<int, int>> ans; // the most accumulated origin position
	int houghTH = HTH; // hough Threshold
};
Hough h;

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


int main() {
	/* open file */
	FILE* fin = fopen(PATH, "rb");
	char outpath[50];
	sprintf(outpath, "%dx%dout.raw", W,H);
	FILE* fout = fopen(outpath, "wb");

	if (fin == NULL || fout == NULL) {
		printf("Failed to open files.\n");
		return 1;
	}

	/* read image */
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			c.in[i][j] = fgetc(fin);
		}
	}


	int dx[25];
	int dy[25];
	for (int i = 0; i < 25; i++) {
		dx[i] = i / 5; // 00000 11111 ...
		dy[i] = i % 5; // 01234 01234 ...
	}
	/* Gaussian filter */
	double GauFlt[5][5];
	getGaussFilter(GauFlt, 5, GAUS);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			for (int dir = 0; dir < 25; dir++) {
				int x = i + dx[dir] - 2;
				int y = j + dy[dir] - 2;
				if (x < 0 || x >= H || y < 0 || y >= W) continue; // zero padding
				c.gauss[i][j] += c.in[x][y] * GauFlt[dx[dir]][dy[dir]];
			}
		}
	}

	/* sobel operation */
	double SobFlt_x[3][3] = {{1, 0, -1},
				     	   {2, 0, -2},
						   {1, 0, -1}};
	double SobFlt_y[3][3] = {{1, 2, 1},
						   {0, 0, 0},
						   {-1, -2, -1}};

	for (int i = 0; i < 9; i++) {
		dx[i] = i / 3; // 000 111 ...
		dy[i] = i % 3; // 012 012 ...
	}

	unsigned int mx = 1;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double gx = 0;
			double gy = 0;
			for (int dir = 0; dir < 9; dir++) {
				int x = i + dx[dir] - 1;
				int y = j + dy[dir] - 1;
				if (x < 0) x = 0; // same_padding
				else if (x >= H) x = H - 1; 
				if (y < 0) y = 0;
				else if (y >= W) y = W - 1;
				gx += (c.gauss[x][y] * SobFlt_x[dx[dir]][dy[dir]]); 
				gy += (c.gauss[x][y] * SobFlt_y[dx[dir]][dy[dir]]);
			}
			c.gradient[i][j] = sqrt(gx * gx + gy * gy); // get gradient magnitude
			mx = max(mx, (unsigned int)c.gradient[i][j]);
			c.theta[i][j] = atan2(gy, gx);
		}
	}
	
/*
	sprintf(outpath, "%dx%dout2.raw", H, W);
	FILE* fout2 = fopen(outpath, "wb");
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			fputc(c.gradient[i][j], fout2);
		}
	}
*/
	/* non maximum supression */
	for (int i = 0; i < H; i++) { 
		for (int j = 0; j < W; j++) {
			int neighbor1, neighbor2;
			int angle = c.theta[i][j] * 180.0 / M_PI;
			if (angle < 0) angle += 180;
			
			if (angle < 22.5 || angle >= 157.5) // compare with neighbor which is in gradient direction  
			{
				neighbor2 = c.gradient[i][j + 1]; 
				neighbor1 = c.gradient[i][j - 1];
			}
			else if (angle >= 22.5 && angle < 67.5) 
			{
				neighbor2 = c.gradient[i - 1][j - 1];
				neighbor1 = c.gradient[i + 1][j + 1];
			}
			else if (angle >= 67.5 && angle < 112.5)
			{
				neighbor2 = c.gradient[i + 1][j];
				neighbor1 = c.gradient[i - 1][j];
			}
			else if (angle >= 112.5 && angle < 157.5)
			{
				neighbor2 = c.gradient[i - 1][j + 1];
				neighbor1 = c.gradient[i + 1][j - 1];
			}

			if (c.gradient[i][j] < max(neighbor1, neighbor2)) { // non maximum : suppression
				c.gradient[i][j] = 0;
			}
			else { // maximum 
				c.out[i][j] = c.gradient[i][j] / mx * 255;
			}
		}	
	}			
	

	/* Two-level Threshold */
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (c.out[i][j] >= c.highTH) { // this coordinate is edge
				c.out[i][j] = 255;
			}
			else if (c.out[i][j] >= c.lowTH) { // weak edge
				c.out[i][j] = 128;
				c.weak.push({ i,j }); // store weak edge to connect 
			}
			else { // not edge
				c.out[i][j] = 0;
				c.gradient[i][j] = 0;
			}
		}
	}


	/* Hysteresis Thresholding */
	pair<int, int> neighbor[8] = { {-1,-1}, {-1, 0}, {-1,1}, {0, -1}, {0, 1}, {1,-1}, {1,0}, {1,1} }; // neighbor direction
	while (!c.weak.empty()) {
		pair<int, int> next = c.weak.top(); c.weak.pop();
		if (c.out[next.first][next.second] == 255) continue;
		for (int i = 0; i < 8; i++) {
			int x = next.first + neighbor[i].first; // get neighbor position
			int y = next.second + neighbor[i].second;
			if (x < 0 || x >= H || y < 0 || y >= W) continue;
			if (c.out[x][y] == 0) continue;
			if (c.out[x][y] == 255) { // connected with strong edge!
				c.out[next.first][next.second] = 255;
			}
			else { // not connected
				c.out[next.first][next.second] = 0;
				c.gradient[next.first][next.second] = 0;
				c.weak.push({ x,y });
			}
		}
	}

	/* vote origin coordinat */
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (c.out[i][j] == 255) {
				double x = i, y	 = j; // set start position
				
				double a = -sin(c.theta[i][j]); // set direction
				double b = -cos(c.theta[i][j]);
				
				int r_range = min(min(x, H - x), min(y, W - y));
				unsigned int reps = r_range / max(abs(a), abs(b));

				for (int k = 0; k < reps; k++) { // draw a line and vote 
					if (x +a*k < 0 || x + a*k >= H || y + b*k < 0 || y + b*k >= W) break;
					h.o[(int)(x + a*k)][(int)(y + b*k)]++;
				}
				a *= -1;
				b *= -1;
				for (int k = 0; k < reps; k++) { // draw a line and vote 
					if (x < 0 || x >= H || y < 0 || y >= W) break;
					h.o[(int)(x + a * k)][(int)(y + b * k)]++;
				}
			}
		}
	}

	
	
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (h.o[i][j] > h.houghTH) {
				h.ans.push_back({i,j}); // candidate for origin coordinate
			}
	//		fputc(c.out[i][j], fout);
		}
	}
	
	mx = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (h.o[i][j])
				fputc(h.o[i][j]*9, fout);
			else
				fputc(0, fout);
			mx = max(mx, h.o[i][j]);
		}
	}
	printf("%d", mx);


	int boundary = 100;//min(W, H) / 2; // boundary for radius space
	vector<vector<int>> R_space(h.ans.size(), vector<int>(boundary, 0)); // hough space for radius
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (c.out[i][j] == 255) { // in edge
				for (int k = 0; k < h.ans.size(); k++) {
					int x, y, r;
					tie(x, y) = h.ans[k];
					int a = i - x;  // get distance from edge to origin
					int b = j - y;
					r = sqrt(a * a + b * b); // get radius
					if (r >= boundary) continue; // radius is too big 
					R_space[k][r]++; // vote in radius space
				}
			}
		}
	}

	/* calculate the answer 
		: if there is overlaped candidate, use average 
		*/
	int n = 0;
	double x_m = 0;
	double y_m = 0;
	double r_m = 0;
	int prev_x = 0, prev_y = 0, prev_r = 0;
	for (int k = 0; k < h.ans.size(); k++) {
		int y, x, r;
		tie(y, x) = h.ans[k];
		int mx = 0;
		for (int i = 0; i < boundary; i++) {
			if (mx < R_space[k][i]) {
				r = i;
				mx = R_space[k][i];
			}
		}
		if (x - prev_x < r && y - prev_y < r) {
			x_m += x;
			y_m += y;
			r_m += r;
			n++;
		}
		else if (x_m > 0) {
			printf("x = %lf, y = %lf, r = %lf\n", x_m / n, y_m / n, r_m / n);
			n = 0;
			x_m = 0;
			y_m = 0;
			r_m = 0;
		}
		prev_x = x;
		prev_y = y;
		prev_r = r;
	}
	if (x_m)
		printf("x = %lf, y = %lf, r = %lf\n", x_m / n, y_m / n, r_m / n);
	
}