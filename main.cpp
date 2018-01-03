#include <string>
#include <iostream>
#include <fstream>
#include "objParser.h"
#include "objViewer.h"

using namespace std;



int main() {
	string objFileName = "bunny.obj";
	float scale = 5;
	string sourceDir = "./resources";

	ifstream fin(sourceDir + "/" + objFileName);
	string fileString((istreambuf_iterator<char>(fin)), istreambuf_iterator<char>());

	// parse obj file
	objParser op = objParser(objFileName);
	op.parse(fileString, scale, 0);

	objViewer ov(1000, 1000, &op);
	
	ov.do_draw();

	cout << "pause" << endl;
}