#include "stringParser.h"
#include "utils.h"

using namespace std;

stringParser::stringParser(string str){
	this->str = str;
	this->index = 0;
}

stringParser::stringParser() {

}
stringParser::~stringParser() {

}

// Initialize StringParser object
void stringParser::init(string str){
	this->str = str;
	this->index = 0;
}

// Skip delimiters
void stringParser::skipDelimiters()  {
	int len = this->str.size();
	int i;
	for(i=this->index; i < len; ++i){
		char c = this->str.at(i);
		if (c == '\t'|| c == ' ' || c == '(' || c == ')' || c == '"') 
			continue;
		break;
	}
	this->index = i;
}

// Skip to the next word
void stringParser::skipToNextWord() {
	this->skipDelimiters();
	int n = getWordLength(this->str, this->index);
	this->index += (n + 1);
}

// Get word
string stringParser::getWord() {
	this->skipDelimiters();
	int n = getWordLength(this->str, this->index);
	if (n == 0) 
		return "";
	string word = this->str.substr(this->index, n);
	this->index += (n + 1);

	return word;
}

// Get integer
int stringParser::getInt() {
	int i;
	i = atoi(this->getWord().c_str());
	return i;
}

// Get floating number
float stringParser::getFloat() {
	float f;
	f = atof(this->getWord().c_str());
	return f;
}