#include "objParser.h"

#include <string>

using namespace std;

objParser::objParser(std::string fileName):fileName(fileName) {

}

objParser::objParser() {
}

objParser::~objParser() {
}

void objParser::parse(string fileString, float scale, int reverse){
	vector<string> lines = split(fileString, "\n");

	OBJObject *currentObject = nullptr;
	string currentMaterialName = "";

	// Parse line by line
	string line;         // A string in the line to be parsed
	stringParser sp;  	 // Create StringParser
	const unsigned long int line_size = lines.size();
	unsigned long int index = 0;
	while (index < line_size) {
		line = lines[index++];
		sp.init(line);                  // init StringParser
		string command = sp.getWord();     // Get command
		if (command.size() == 0)
			continue;  // check null command

		if (command == "#") {
			continue;
		}
		else if (command == "mtllib") {
			string path = this->parseMtllib(sp, this->fileName);
			// var mtl = new MTLDoc();   // Create MTL instance
			// this.mtls.push(mtl);
			// var request = new XMLHttpRequest();
			// request.onreadystatechange = function() {
			// if (request.readyState == 4) {
			//   if (request.status != 404) {
			//     onReadMTLFile(request.responseText, mtl);
			//   }else{
			//     mtl.complete = true;
			//   }
			// }
			// }
			// request.open('GET', path, true);  // Create a request to acquire the file
			// request.send();                   // Send the request
			// continue; // Go to the next line
		}
		else if (command == "o" || command == "g") {
			OBJObject object = this->parseObjectName(sp);
			this->objects.push_back(object);
			currentObject = &this->objects.back();
			continue; // Go to the next line
		}
		else if (command == "v") {
			Vertex vertex = this->parseVertex(sp, scale);
			this->vertices.push_back(vertex);
			continue; // Go to the next line
		}
		else if (command == "vn") {
			Normal normal = this->parseNormal(sp);
			this->normals.push_back(normal);
			continue;
		}
		else if (command == "usemtl") {
			currentMaterialName = this->parseUsemtl(sp);
			continue;
		}
		else if (command == "f") {
			if (currentObject == nullptr) {
				OBJObject object;
				object.name = "default";
				object.numIndices = 0;
				this->objects.push_back(object);
				currentObject = &this->objects.back();
			}

			vector<Face> faces = this->parseFace(sp, currentMaterialName, this->vertices, reverse);
			for (int i = 0; i < faces.size(); ++i) {
				currentObject->faces.push_back(faces[i]);
				currentObject->numIndices += 1;
			}
		}
	}
}

string objParser::parseMtllib(stringParser &sp, string &fileName){
	// Get directory path
	int i = fileName.rfind('/', fileName.size()-1);
	string dirPath = "";
	if(i >= 0)
		dirPath = fileName.substr(0, i+1);

	return dirPath + sp.getWord();   // Get path
}

OBJObject objParser::parseObjectName(stringParser &sp){
	OBJObject obj;
	obj.name = sp.getWord();
	obj.numIndices = 0;
	return obj;
}

Vertex objParser::parseVertex(stringParser sp, float scale) {
	float x = sp.getFloat() * scale;
	float y = sp.getFloat() * scale;
	float z = sp.getFloat() * scale;
	return Vertex{x, y, z};
}

Normal objParser::parseNormal(stringParser sp) {
	float x = sp.getFloat();
	float y = sp.getFloat();
	float z = sp.getFloat();
	return Normal{x, y, z};
}

string objParser::parseUsemtl(stringParser sp) {
	return sp.getWord();
}

vector<Face> objParser::parseFace(stringParser sp, string materialName, vector<Vertex> &vertices, int reverse) {
	float vIndices[4];
	float nIndices[4];

	// ge indices
	unsigned long int index = 0;
	while(true) {
		string word = sp.getWord();
		if (word == "")
			break;
		vector<string> subWords = split(word, "/");
		if (subWords.size() >= 1) {
			int vi = atoi(subWords[0].c_str()) - 1;
			vIndices[index] = vi;
		}
		if (subWords.size() >= 2) {
			int ni = atoi(subWords[1].c_str()) - 1;
			nIndices[index] = ni;
		}
		else {
			nIndices[index] = -1;
		}
		index += 1;
	}

	vector<Face> faces;

	if (index == 3) {
		Face face;
		face.materialName = materialName;

		for (int i = 0; i < 3; ++i) {
			face.nIndices[i] = nIndices[i];
			face.vIndices[i] = vIndices[i];
		}

		faces.push_back(face);
	}
	else if (index == 4) {
		Face face;
		face.materialName = materialName;

		for (int i = 0; i < 3; ++i) {
			face.nIndices[i] = nIndices[i];
			face.vIndices[i] = vIndices[i];
		}

		faces.push_back(face);

		face.nIndices[0] = nIndices[2];
		face.vIndices[0] = vIndices[2];
		face.nIndices[1] = nIndices[3];
		face.vIndices[1] = vIndices[3];
		face.nIndices[2] = nIndices[0];
		face.vIndices[2] = vIndices[0];

		faces.push_back(face);
	}

	return faces;
}








