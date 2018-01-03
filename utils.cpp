#include "utils.h"
#include <math.h>

using namespace std;

vector<string> split(const string &s, const string &seperator){
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while(i != s.size()){
		//找到字符串中首个不等于分隔符的字母；
		int flag = 0;
		while(i != s.size() && flag == 0){
		flag = 1;
		for(string_size x = 0; x < seperator.size(); ++x)
			if(s[i] == seperator[x]){
				++i;
				flag = 0;
				break;
			}
		}
		
		//找到又一个分隔符，将两个分隔符之间的字符串取出；
		flag = 0;
		string_size j = i;
		while(j != s.size() && flag == 0){
		for(string_size x = 0; x < seperator.size(); ++x)
			if(s[j] == seperator[x]){
				flag = 1;
				break;
			}
		if(flag == 0) 
			++j;
		}
		if(i != j){
			result.push_back(s.substr(i, j-i));
			i = j;
		}
	}
	return result;
}

int getWordLength(const string &str, int start){
	int n = 0;
	int len = str.size();
	int i;
	for(i = start; i < len; i++){
		char c = str.at(i);
		if (c == '\t'|| c == ' ' || c == '(' || c == ')' || c == '"') 
			break;
	}
	return i - start;
}

void cross_product(const Vector3d &a, const Vector3d &b, Vector3d &c) {
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
}
double dot_product(const Vector3d &a, const Vector3d &b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
void normalize(Vector3d &a) {
	double len = sqrt(dot_product(a, a));
	a.x /= len;
	a.y /= len;
	a.z /= len;
}

double intersect(double x1, double y1, double x2, double y2, double inter_y) {
	return x1 + (x2 - x1)*(inter_y - y1) / (y2 - y1);
}

double cal_plane_z(double a, double b, double c, double d, double x, double y) {
	return (-d - a*x - b*y) / c;
}


void color_add(const Color &a, const Color &b, Color &c) {
	c.r = a.r + b.r;
	c.g = a.g + b.g;
	c.b = a.b + b.b;
	c.a = a.a + b.a;
}
void color_sub(const Color &a, const Color &b, Color &c) {
	c.r = a.r - b.r;
	c.g = a.g - b.g;
	c.b = a.b - b.b;
	c.a = a.a - b.a;
}
void color_scale(const Color &a, double scale, Color &c) {
	c.r = a.r * scale;
	c.g = a.g * scale;
	c.b = a.b * scale;
	c.a = a.a * scale;
}
void color_assign(const Color &a, Color &c) {
	c.r = a.r;
	c.g = a.g;
	c.b = a.b;
	c.a = a.a;
}

