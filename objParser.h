#pragma once

#include <vector>
#include <string>
#include "utils.h"
#include "stringParser.h"

using namespace std;

class objParser
{
public:
	string fileName;
	vector<OBJObject> objects;
	vector<Vertex> vertices;
	vector<Normal> normals;


	objParser(string fileName);
	objParser();
	~objParser();
	
	void parse(string fileString, float scale, int reverse);
	string parseMtllib(stringParser &sp, string &fileName);
	OBJObject parseObjectName(stringParser &sp);
	Vertex parseVertex(stringParser sp, float scale);
	Normal parseNormal(stringParser sp);
	string parseUsemtl(stringParser sp);
	vector<Face> parseFace(stringParser sp, string materialName, vector<Vertex> &vertices, int reverse);

};