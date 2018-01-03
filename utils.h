#pragma once

#include <string>
#include <vector>

using namespace std;

struct Vertex {
	float x;
	float y;
	float z;
} ;

struct Normal {
	float x;
	float y;
	float z;
} ;

struct Color {
	float r;
	float g;
	float b;
	float a;
} ;

struct Face {
	string materialName;
	int vIndices[3];
	int nIndices[3];
};

struct OBJObject {
	string name;
	vector<Face> faces;
	int numIndices;
};

struct Vector3d {
	float x;
	float y;
	float z;
};

vector<string> split(const string &s, const string &seperator);
int getWordLength(const string &str, int start);
void cross_product(const Vector3d &a, const Vector3d &b, Vector3d &c);
double dot_product(const Vector3d &a, const Vector3d &b);
void normalize(Vector3d &a);
double intersect(double x1, double y1, double x2, double y2, double inter_y);
double cal_plane_z(double a, double b, double c, double d, double x, double y);

void color_add(const Color &a, const Color &b, Color &c);
void color_sub(const Color &a, const Color &b, Color &c);
void color_scale(const Color &a, double scale, Color &c);
void color_assign(const Color &a, Color &c);