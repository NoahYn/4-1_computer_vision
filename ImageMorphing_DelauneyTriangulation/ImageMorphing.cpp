#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <list>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

// Define macros for accessing points of a triangle for source (sTri), destination (dTri), and intermediate (iTri) 
#define X0 sTri[0].x
#define X1 sTri[1].x 
#define X2 sTri[2].x
#define Y0 sTri[0].y
#define Y1 sTri[1].y
#define Y2 sTri[2].y

#define U0 dTri[0].x
#define U1 dTri[1].x
#define U2 dTri[2].x
#define V0 dTri[0].y
#define V1 dTri[1].y
#define V2 dTri[2].y

#define A0 iTri[0].x
#define A1 iTri[1].x
#define A2 iTri[2].x
#define B0 iTri[0].y
#define B1 iTri[1].y
#define B2 iTri[2].y

// mPoint structure is used to store a 2D point along with a match parameter.
struct mPoint : Point2i {
	int match;

	// Different constructors for various initialization scenarios.
	mPoint() : Point2i(), match(-1) {}
	mPoint(int _x, int _y) : Point2i(_x, _y), match(-1) {}
	mPoint(int _x, int _y, int _m) : Point2i(_x, _y), match(_m) {}
	mPoint(Point p) : Point2i(p.x, p.y), match(-1) {}
	mPoint(Point p, int _m) : Point2i(p.x, p.y), match(_m) {}

	// Method to return dot product of the point with itself.
	int dotProd() const {
		return this->ddot(*this);
	}

	// Methods to get the Euclidean distance between the current point and another point or coordinates.
	double getDist(double x, double y) const {
		double xx = this->x - x;
		double yy = this->y - y;

		return sqrt(xx * xx + yy * yy);
	}

	double getDist(const Point& other) const {
		double x = this->x - other.x;
		double y = this->y - other.y;

		return sqrt(x * x + y * y);
	}
};

// mEdge structure to represent an edge of a triangle with mPoints and a isBad flag.
struct mEdge {
	mPoint p1;
	mPoint p2;
	bool isBad;

	mEdge(const mPoint& _p1, mPoint& _p2) {
		this->p1 = _p1;
		this->p2 = _p2;
		this->isBad = false;
	}

	// Overloaded equality operator to allow comparison of mEdge objects.
	bool operator==(const mEdge& other) const {
		return (p1 == other.p1 && p2 == other.p2) || (p1 == other.p2 && p2 == other.p1);
	}
};

// mTriangle structure to represent a triangle composed of mPoints.
struct mTriangle {
	mPoint p1;
	mPoint p2;
	mPoint p3;

	mTriangle() = default;

	mTriangle(const mPoint& v1, const mPoint& v2, const mPoint& v3) {
		this->p1 = v1;
		this->p2 = v2;
		this->p3 = v3;
	}

	// Method checks if the given mPoint lies on the circle circumscribing the mTriangle.
	bool isOnCircle(const mPoint& other) const {
		// dot product each point with itself
		int dot1 = p1.dotProd();
		int dot2 = p2.dotProd();
		int dot3 = p3.dotProd();

		// get coordinate from each point
		double x1 = p1.x;
		double y1 = p1.y;
		double x2 = p2.x;
		double y2 = p2.y;
		double x3 = p3.x;
		double y3 = p3.y;

		// get the origin of circumcircle by formular
		double circum_x = (dot1 * (y3 - y2) + dot2 * (y1 - y3) + dot3 * (y2 - y1)) / (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1));
		double circum_y = (dot1 * (x3 - x2) + dot2 * (x1 - x3) + dot3 * (x2 - x1)) / (y1 * (x3 - x2) + y2 * (x1 - x3) + y3 * (x2 - x1));
		double circum_center_x = circum_x / 2;
		double circum_center_y = circum_y / 2;

		// get the radius : radius is distance between one vertex of triangle and origin
		double circum_radius = p1.getDist(circum_center_x, circum_center_y);

		// get the distance between point and origin of circumcircle
		double dist = other.getDist(circum_center_x, circum_center_y);

		// if distance is in radius, the point is on circle
		return dist <= circum_radius;
	}

	// Method checks if the given mPoint is one of the vertices of the mTriangle.
	bool onPoint(mPoint& other) {
		return p1 == other || p2 == other || p3 == other;
	}
};

inline void matmul3x3(double a[3][3], double b[3][3], double res[3][3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			res[i][j] = 0;
			for (int k = 0; k < 3; k++)
				res[i][j] += a[i][k] * b[k][j];
		}
	}
}

inline void matmul13(double a[3][3], double b[3], double res[3]) {
	for (int i = 0; i < 3; i++) {
		res[i] = 0;
		for (int j = 0; j < 3; j++) {
			res[i] += a[i][j] * b[j];
		}
	}
}

inline void matinv3x3(double a[3][3], double inv[3][3]) {
	double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
		- a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
		+ a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

	if (det == 0) return;

	inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
	inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
	inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
	inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
	inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
	inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
	inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
	inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
	inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
}

// matrix multiply 
inline void matTrs3x3(double a[3][3], double trs[3][3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			trs[i][j] = a[j][i];
		}
	}
}

auto DelaunayTriMesh(vector<mPoint> pts, int row, int col) {
	list<mTriangle> tris;

	// make super tris to start with convex hull(contain all point)
	mPoint superTriA = { -1000, -1000 };
	mPoint superTriB = { 2 * col + 1000, -1000 };
	mPoint superTriC = { -1000, 2 * row + 1000 };
	mTriangle superTri = { superTriA, superTriB, superTriC };

	tris.push_back(superTri);// start with super triangle
	for (auto& point : pts) { 
		vector<mEdge> edges;
		tris.remove_if([&](auto& triangle) {
			auto f = triangle.isOnCircle(point); // remove traingle if point is on circle
			if (f) { // get edge of current triangle if point is on circle
				edges.push_back({ triangle.p1, triangle.p2 }); 
				edges.push_back({ triangle.p2, triangle.p3 });
				edges.push_back({ triangle.p3, triangle.p1 });
			}
			return f;
		});

		for (auto e1 = begin(edges); e1 != end(edges); ++e1) // if two same edges are found, check them as bad 
		{
			for (auto e2 = e1 + 1; e2 != end(edges); ++e2)
			{
				if (*e1 == *e2)
				{
					e1->isBad = true;
					e2->isBad = true;
				}
			}
		}

		edges.erase(std::remove_if(begin(edges), end(edges), [](mEdge& e) { // remove bad edge. 
			return e.isBad;
			}), end(edges));

		for (const auto& e : edges) { // make triangles from remain edges
			tris.push_back(
				mTriangle(e.p1, e.p2, point)
			);
		}
	}

	tris.remove_if([&](auto& triangle) { // remove vertex on super triangle
		return triangle.onPoint(superTriA) ||
			triangle.onPoint(superTriB) ||
			triangle.onPoint(superTriC);
		});
	return tris;
}

void drawLine(cv::Mat img, list<mTriangle>& triangles) {
	for (auto& t : triangles) {
		cv::line(img, t.p1, t.p2, cv::Scalar(255, 255, 0), 1, cv::LINE_AA);
		cv::line(img, t.p2, t.p3, cv::Scalar(255, 255, 0), 1, cv::LINE_AA);
		cv::line(img, t.p3, t.p1, cv::Scalar(255, 255, 0), 1, cv::LINE_AA);
	}
	cv::imshow("", img);
	cv::waitKey(0);
}

inline bool isOnTri(const Point &p, const vector<Point2f> &tri) {
	int s1 = (tri[1].y - tri[0].y) * (p.x - tri[0].x) - (tri[1].x - tri[0].x) * (p.y - tri[0].y);
	int s2 = (tri[2].y - tri[1].y) * (p.x - tri[1].x) - (tri[2].x - tri[1].x) * (p.y - tri[1].y);
	int s3 = (tri[0].y - tri[2].y) * (p.x - tri[2].x) - (tri[0].x - tri[2].x) * (p.y - tri[2].y);
	return (s1 >= 0 && s2 >= 0 && s3 >= 0) || (s1 <= 0 && s2 <= 0 && s3 <= 0);
}

// Opencv version of morphing to test the result
void affineCV(vector<Point2f>& tri1, vector<Point2f>& tri2, Mat& src, Mat& dst)
{
	Mat transMat = getAffineTransform(tri1, tri2);
	warpAffine(src, dst, transMat, dst.size(), INTER_LINEAR, BORDER_REFLECT_101);
}

void triMorph(vector<Point2f>& sTri, vector<Point2f>& dTri, vector<Point2f>& iTri, Mat& imS, Mat& imD, Mat& res, double ratio)
{
	double transMat[2][3][3];
	double coeff[3][3];
	double Mat_inv[3][3];

	double sMat[3][3] = {
		{X0, Y0, 1},
		{X1, Y1, 1},
		{X2, Y2, 1}
	};
	double dMat[3][3] = {
		{U0, V0, 1},
		{U1, V1, 1},
		{U2, V2, 1}
	};
	double iMat[3][3] = {
		{A0, B0, 1},
		{A1, B1, 1},
		{A2, B2, 1}
	};
	// get the affine transMat for s -> i
	matinv3x3(iMat, Mat_inv);
	matmul3x3(Mat_inv, sMat, coeff);
	matTrs3x3(coeff, transMat[0]);

	// transMat for d -> i
	matmul3x3(Mat_inv, dMat, coeff);
	matTrs3x3(coeff, transMat[1]);

	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 3; j++) {
	//		cout << transMat[1][i][j] << ' ';
	//	}cout << '\n';
	//}
	//cout << "\n";

	// get the minimum and maximum of x and y for range
	Point mp(min({ iTri[0].x, iTri[1].x, iTri[2].x }), min({ iTri[0].y, iTri[1].y, iTri[2].y }));
	Point Mp(max({ iTri[0].x, iTri[1].x, iTri[2].x }), max({ iTri[0].y, iTri[1].y, iTri[2].y }));
	if (Mp.x >= res.cols) Mp.x = res.cols - 1; // can not exceed the image size
	if (Mp.y >= res.rows) Mp.y = res.rows - 1;

	for (int y = mp.y; y <= Mp.y; y++) {
		for (int x = mp.x; x <= Mp.x; x++) {
			Point cp(x, y);
			if (isOnTri(cp, iTri)) { // point is on triangle
				double d[3] = { cp.x, cp.y, 1 }; 
				double s[3] = {};
				matmul13(transMat[0], d, s); // use the affine translation matrix
				Point sp(s[0], s[1]);
				matmul13(transMat[1], d, s);
				Point dp(s[0], s[1]);

				Vec3f a = (1 - ratio) * imS.at<Vec3f>(sp.y, sp.x);
				Vec3f b = ratio * imD.at<Vec3f>(dp.y, dp.x);
				res.at<Vec3f>(cp.y, cp.x) = a + b; // get value from two src image 
			}
		}
	}


	//// get rectangles fully containing triangles
	//Rect sRect = boundingRect(sTri);
	//Rect iRect = boundingRect(iTri);
	//Rect dRect = boundingRect(dTri);

	//Mat sRec = img1(sRect);
	//Mat dRec = img2(dRect);
	//Mat iRec = res(iRect);
	//if (!sRec.rows || !sRec.cols || !dRec.rows || !dRec.cols || !iRec.rows || !iRec.cols)
	//	return;

	////initialize mask to use only triangle in rectangle
	//float mx = img1.cols; float my = img1.rows;
	//for (int i = 0; i < sTri.size(); i++) {
	//	mx = min(mx, sTri[i].x);   my = min(my, sTri[i].y);
	//}
	//Point sp(mx, my);

	//mx = img1.cols; my = img1.rows;
	//for (int i = 0; i < dTri.size(); i++) {
	//	mx = min(mx, dTri[i].x);   my = min(my, dTri[i].y);
	//}
	//Point dp(mx, my);

	//mx = img1.cols; my = img1.rows;
	//for (int i = 0; i < iTri.size(); i++) {
	//	mx = min(mx, iTri[i].x);   my = min(my, iTri[i].y);
	//}
	//Point ip(mx, my);

	//vector <Point> maskTri;
	//for (int i = 0; i < 3; i++)	{
	//	iTri[i].x = iTri[i].x - ip.x;      iTri[i].y = iTri[i].y - ip.y;
	//	sTri[i].x = sTri[i].x - sp.x;      sTri[i].y = sTri[i].y - sp.y;
	//	dTri[i].x = dTri[i].x - dp.x;      dTri[i].y = dTri[i].y - dp.y;
	//	maskTri.push_back(Point(iTri[i].x, iTri[i].y));
	//}

	//Mat mask = Mat::zeros(iRec.rows, iRec.cols, CV_32FC3);
	//fillConvexPoly(mask, maskTri, Scalar(1., 1., 1.), 16, 0);

	//// get warped matrices
	//Mat sWarped = Mat::zeros(iRec.rows, iRec.cols, sRec.type());
	//Mat dWarped = Mat::zeros(iRec.rows, iRec.cols, dRec.type());

	//// affine transform 
	//affineTransformTri(sTri, iTri, sRec, sWarped);
	//affineTransformTri(dTri, iTri, dRec, dWarped);
	//// open cv version
	////affineCV(sTri, iTri, sRec, sWarped);
	////affineCV(dTri, iTri, dRec, dWarped);

	//// cross-dissolve two warped triangle with ratio
	//Mat imgRect = ratio * dWarped + (1.0 - ratio) * sWarped;

	//// multiply img and mask
	//multiply(imgRect, mask, imgRect);

	//mask = Scalar(1.0, 1.0, 1.0) - mask;
	//Rect resRec = Rect(ip.x, ip.y, iRec.cols, iRec.rows);

	//// multiply img and mask with offset
	//for (int y = ip.y; y < ip.y + mask.rows; y++) {
	//	for (int x = ip.x; x < ip.x + mask.cols; x++) {
	//		res.at<Vec3f>(y, x)[0] = res.at<Vec3f>(y, x)[0] * mask.at<Vec3f>(y - ip.y, x - ip.x)[0];
	//		res.at<Vec3f>(y, x)[1] = res.at<Vec3f>(y, x)[1] * mask.at<Vec3f>(y - ip.y, x - ip.x)[0];
	//		res.at<Vec3f>(y, x)[2] = res.at<Vec3f>(y, x)[2] * mask.at<Vec3f>(y - ip.y, x - ip.x)[0];
	//	}
	//}

	//// add resRec to img
	//res(resRec) = res(resRec) + imgRect;
}

int main() {
	// file paths
	vector<string> imagefiles = { "1_600x600.png", "2_600x600.png", "3_480x480.png", "4_480x480.png", "5_500x500.png", "6_500x500.png" };
	vector<string> pointfiles = { "1_600x600.txt", "2_600x600.txt", "3_480x480.txt", "4_480x480.txt", "5_500x500.txt", "6_500x500.txt" };
	string imagesPath = "images/";
	string pointsPath = "labels/";
	string outputPath = "output/";

	// vector to contain imgs and points
	vector<Mat> imgs;
	vector<vector<mPoint>> ptss;
	vector<list<mTriangle>> tris;

	// initialize
	for (int i = 0; i < imagefiles.size(); i++) { // input
		int x, y, match;	
		char comma;		
		string temp;
		ifstream ifs(pointsPath + pointfiles[i]); // points file open
		vector<mPoint> pts;
		
		imgs.push_back(imread(imagesPath + imagefiles[i])); // image read

		getline(ifs, temp); // throw away format info 
		while (ifs >> y >> comma >> x >> comma >> match) // read points and matching
			pts.push_back(mPoint(x, y, match-1)); // push to vector
		ptss.push_back(pts);
		
		list<mTriangle> tri = DelaunayTriMesh(pts, imgs[i].rows, imgs[i].cols);
		tris.push_back(tri);
		
		ifs.close();
	}

	// warping
	while(1){
		int choose1, choose2;
		cout << "choose two image by number.(1-6) ex) 1 2 \n";
		cin >> choose1 >> choose2; // get the two of number of image from user
		if (choose1 < 0 || choose2 < 0 || choose1 > 6 || choose2 > 6) {
			cout << "you should choose image in range 1 to 6\n";
			continue;
		}
		choose1--; choose2--; // index range is 0 to 5

		//for (int j = 0; j < points[i].size(); j++) { // draw point on each feature point
		//	Point pt1 = points[i][j], pt2 = points[i + 1][j];
		//	circle(imgs[i], pt1, 3, Scalar(0, 0, 255), -1);
		//	circle(imgs[i + 1], pt2, 3, Scalar(0, 0, 255), -1);
		//}

		vector<mPoint>& sPoints = ptss[choose1]; // source of point
		vector<mPoint>& dPoints = ptss[choose2]; // dest 
		list<mTriangle>& sTris = tris[choose1];

		// get the new points for morphed image based on ratio
		for (double ratio = 0; ratio < 1; ratio += 0.05) {
			Mat &imS = imgs[choose1]; // source of image
			Mat &imD = imgs[choose2]; // destination
			imS.convertTo(imS, CV_32FC3); // convert to float for morphing operations
			imD.convertTo(imD, CV_32FC3);
			int rows = imD.rows * ratio + imS.rows * (1 - ratio);
			int cols = imD.cols * ratio + imS.cols * (1 - ratio);
			Mat morphedImage(rows, cols, CV_32FC3); // allocate space for final image

			//vector<mPoint> iPoints; // middle of two point
			//for (int j = 0; j < sPoints.size(); j++)
			//{
			//	Point pi;
			//	pi.x = ratio * dPoints[j].x + (1 - ratio) * sPoints[j].x;
			//	pi.y = ratio * dPoints[j].y + (1 - ratio) * sPoints[j].y;
			//	iPoints.push_back(mPoint(pi, sPoints[j].match));
			//}

			//draw_line(imgs[choose1], sTris);
			//draw_line(imgs[choose2], dTris);

			//opencv trinagulation 
			/* 
			Subdiv2D sub_s, sub_m, sub_d;
			sub_s.initDelaunay(Rect(0, 0, imS.cols, imS.rows));
			sub_m.initDelaunay(Rect(0, 0, imS.cols, imS.rows));
			sub_d.initDelaunay(Rect(0, 0, imD.cols, imD.rows));
			// add point to subdiv2d
			sub_s.insert(sPoints);
			sub_m.insert(iPoints);
			sub_d.insert(dPoints);
			// triangluation
			vector<Vec6f> tri_s, tri_m, tri_d;
			sub_s.getTriangleList(tri_s);
			sub_m.getTriangleList(tri_m);
			sub_d.getTriangleList(tri_d);
			*/
		
			for (auto& ks : sTris) {
				vector<Point2f> sTri(3), iTri(3), dTri(3);
				sTri = { Point2f(ks.p1), Point2f(ks.p2), Point2f(ks.p3) }; 
				dTri = { dPoints[ks.p1.match], dPoints[ks.p2.match], dPoints[ks.p3.match] };
				iTri = { dTri[0] * ratio + sTri[0] * (1 - ratio),
						dTri[1] * ratio + sTri[1] * (1 - ratio),
						dTri[2] * ratio + sTri[2] * (1 - ratio) };
				triMorph(sTri, dTri, iTri, imS, imD, morphedImage, ratio);
			}
			
			imshow("warped image", morphedImage/255.0);
			waitKey(3);
		}
	}
	return 0;
}