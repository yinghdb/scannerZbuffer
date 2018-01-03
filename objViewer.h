#pragma once
#define GLEW_STATIC

#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "objParser.h"
#include "utils.h"
#include <GL\glew.h>
#include <glfw3.h>


class objViewer {
private:
	objParser *op;
	unsigned long int vertex_num;
	unsigned long int normal_num;
	float *proj_coord;
	float *norm_vec;
	float *vColor;
	unsigned long int face_num;
	Face *face_list;

	int wndWidth;
	int wndHeight;
	double nearZ;
	double farZ;
	
	GLuint bufferObj;
	GLFWwindow *window;

	unsigned char *frame_buffer;

	void gen_face_list();
	void scan_to_image();
	void cal_proj();

public:
	objViewer(int w, int h, objParser *op);
	~objViewer();

	void do_draw();
};