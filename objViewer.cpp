#include "objViewer.h"
#include "scannerZbuffer.h"

#ifndef NDEBUG
#pragma comment(lib, "GLEW_1130.lib")
#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "opengl32.lib")
#else
#pragma comment(lib, "GLEW_1130.lib")
#pragma comment(lib, "glfw3.lib")
#pragma comment(lib, "opengl32.lib")
#endif // NDEBUG

namespace OpenGLUtils {
	bool left_button_pressed = false;
	double last_x, last_y;
	glm::vec3 camera_pos = glm::vec3(0.0, 0.0, 1.0);
	glm::vec3 camera_lookAt = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 obj_pos = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 light_pos = glm::vec3(1.0, 1.0, 1.0);
	float angle_y = 0;
	float angle_x = 0;
	float angle_y_cam = 0;
	float angle_x_cam = 0;

	void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
	void cursor_position_callback(GLFWwindow* window, double x, double y);
}

void OpenGLUtils::mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
	if (action == GLFW_PRESS)
		switch (button)
		{
		case GLFW_MOUSE_BUTTON_LEFT:
			left_button_pressed = true;
			glfwGetCursorPos(window, &last_x, &last_y);
			break;
		case GLFW_MOUSE_BUTTON_MIDDLE:
			break;
		case GLFW_MOUSE_BUTTON_RIGHT:
			break;
		default:
			return;
		}
	else if (action == GLFW_RELEASE)
		switch (button)
		{
		case GLFW_MOUSE_BUTTON_LEFT:
			left_button_pressed = false;
			break;
		case GLFW_MOUSE_BUTTON_MIDDLE:
			break;
		case GLFW_MOUSE_BUTTON_RIGHT:
			break;
		default:
			return;
		}
	return;
}

void OpenGLUtils::cursor_position_callback(GLFWwindow* window, double x, double y) {
	if (left_button_pressed) {
		float ox = x - last_x;
		float oy = y - last_y;

		//angle_y += ox / 100;
		//angle_x += oy / 100;

		angle_y_cam -= ox / 100;
		angle_x_cam += oy / 100;

		camera_pos = glm::vec3(cos(angle_x_cam)*sin(angle_y_cam), sin(angle_x_cam), cos(angle_x_cam)*cos(angle_y_cam));

		last_x = x;
		last_y = y;
	}
}

objViewer::objViewer(int w, int h, objParser *op)
	:wndHeight(h), wndWidth(w), op(op), proj_coord(nullptr), nearZ(0), farZ(2)
{
	gen_face_list();

	vertex_num = op->vertices.size();
	normal_num = op->normals.size();
	proj_coord = new float[vertex_num * 3];
	norm_vec = new float[normal_num * 3];
	vColor = new float[face_num * 3 * 4];
	frame_buffer = new unsigned char[1024 * 1024 * 4];

	// Initialise GLFW
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		getchar();
		return;
	}

	window = glfwCreateWindow(wndWidth, wndHeight, "Scanner Z-buffer", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		getchar();
		glfwTerminate();
		return;
	}
	glfwMakeContextCurrent(window);

	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		return;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	// bind event
	glfwSetCursorPosCallback(window, OpenGLUtils::cursor_position_callback);
	glfwSetMouseButtonCallback(window, OpenGLUtils::mouse_button_callback);

	glGenBuffers(1, &bufferObj);
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, bufferObj);
	glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, wndWidth*wndHeight * 4, frame_buffer, GL_DYNAMIC_DRAW_ARB);
}

objViewer::~objViewer() {
	glDeleteBuffers(1, &bufferObj);
	if (proj_coord != nullptr)
		delete[] proj_coord;
	if (norm_vec != nullptr)
		delete[] norm_vec;
	if (vColor != nullptr)
		delete[] vColor;
}

void objViewer::cal_proj() {
	glm::mat4 R(1.0);
	R = glm::rotate(R, OpenGLUtils::angle_x, glm::vec3(1, 0, 0));
	R = glm::rotate(R, OpenGLUtils::angle_y, glm::vec3(0, 1, 0));

	glm::mat4 M(1.0);
	M = glm::translate(M, OpenGLUtils::obj_pos);
	M = glm::rotate(M, OpenGLUtils::angle_x, glm::vec3(1, 0, 0));
	M = glm::rotate(M, OpenGLUtils::angle_y, glm::vec3(0, 1, 0));

	// calculate vertexes and normals after Model traformation
	for (unsigned long int i = 0; i < normal_num; ++i) {
		Normal norm = op->normals.at(i);
		glm::vec4 n_ori(norm.x, norm.y, norm.z, 1);
		glm::vec4 n_rot = R * n_ori;

		glm::vec3 norm_dir = glm::normalize(glm::vec3(n_rot));

		norm_vec[3 * i + 0] = norm_dir.x;
		norm_vec[3 * i + 1] = norm_dir.y;
		norm_vec[3 * i + 2] = norm_dir.z;
	}

	for (unsigned long int i = 0; i < vertex_num; ++i) {
		Vertex v_ori = op->vertices.at(i);
		glm::vec4 v_glm(v_ori.x, v_ori.y, v_ori.z, 1);
		glm::vec4 v_trans = M * v_glm;

		proj_coord[3 * i + 0] = v_trans.x;
		proj_coord[3 * i + 1] = v_trans.y;
		proj_coord[3 * i + 2] = v_trans.z;
	}

	for (unsigned long int i = 0; i < face_num; ++i) {
		int *vIndices = face_list[i].vIndices;
		int *nIndices = face_list[i].nIndices;

		for (int j = 0; j < 3; ++j) {
			int vid = vIndices[j];
			int nid = nIndices[j];

			glm::vec3 vertex = glm::vec3(proj_coord[3 * vid + 0], proj_coord[3 * vid + 1], proj_coord[3 * vid + 2]);
			glm::vec3 normal = glm::vec3(norm_vec[3 * nid + 0], norm_vec[3 * nid + 1], norm_vec[3 * nid + 2]);

			glm::vec3 light_dir = glm::vec3(0.0, 0.0, 1.0);  // glm::normalize(OpenGLUtils::light_pos - vertex);
			double cos = glm::dot(normal, light_dir);

			if (cos < 0) {
				vColor[12 * i + 4 * j + 0] = 1;
				vColor[12 * i + 4 * j + 1] = 0;
				vColor[12 * i + 4 * j + 2] = 0;
				vColor[12 * i + 4 * j + 3] = 1;
			}
			else {
				vColor[12 * i + 4 * j + 0] = cos;
				vColor[12 * i + 4 * j + 1] = cos;
				vColor[12 * i + 4 * j + 2] = cos;
				vColor[12 * i + 4 * j + 3] = 1;
			}
		}
	}

	glm::mat4 V = glm::lookAt(OpenGLUtils::camera_pos, OpenGLUtils::camera_lookAt, glm::vec3(0, 1, 0));
	glm::mat4 P = glm::ortho(-1.0, 1.0, -1.0, 1.0, nearZ, farZ);
	//glm::mat4 P = glm::frustum(-0.1, 0.1, -0.1, 0.1, 0.1, 10.0);

	glm::mat4 MVP = P * V * M;

	for (unsigned long int i = 0; i < vertex_num; ++i) {
		Vertex v_ori = op->vertices.at(i);
		glm::vec4 v_glm(v_ori.x, v_ori.y, v_ori.z, 1);
		glm::vec4 v_trans = MVP * v_glm;

		proj_coord[3 * i + 0] = v_trans.x / v_trans.w;
		proj_coord[3 * i + 1] = v_trans.y / v_trans.w;
		proj_coord[3 * i + 2] = - v_trans.z / v_trans.w;
	}
}

void objViewer::do_draw() {
	do {
		cal_proj();
		scan_to_image();

		glClearColor(0.0, 0.0, 0.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);

		glDrawPixels(wndWidth, wndHeight, GL_RGBA, GL_UNSIGNED_BYTE, 0);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();
	} while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);
}

void objViewer::gen_face_list() {
	face_num = 0;
	for (auto obj : op->objects) {
		face_num += obj.numIndices;
	}

	face_list = new Face[face_num];
	unsigned long int fid = 0;
	for (auto obj : op->objects) {
		for (auto face : obj.faces) {
			face_list[fid] = face;
			fid += 1;
		}
	}
}

void objViewer::scan_to_image() {
	scannerZbuffer sZ(proj_coord, vColor, vertex_num, face_list, face_num, wndWidth, wndHeight);
	//areaScannerZbuffer sZ(proj_coord, vColor, vertex_num, face_list, face_num, wndWidth, wndHeight);
	sZ.set_cam_z(nearZ, farZ);
	sZ.init();
	sZ.scan_to_image(frame_buffer);

	glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, bufferObj);
	glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, wndWidth*wndHeight * 4, frame_buffer, GL_DYNAMIC_DRAW_ARB);
}

