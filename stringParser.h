#pragma once

#include <string>
using namespace std;

class stringParser
{
private:
	string str;
	int index;
public:
	stringParser(string str);
	stringParser();
	~stringParser();

	void init(string str);

	// Skip delimiters
	void skipDelimiters();

	// Skip to the next word
	void skipToNextWord();

	// Get word
	string getWord();

	// Get integer
	int getInt();

	// Get floating number
	float getFloat();
};
