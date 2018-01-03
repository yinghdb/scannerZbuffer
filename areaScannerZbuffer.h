#pragma once
#include "objParser.h"
#include "utils.h"
#include <memory>

using namespace std;

namespace area_scanner {
	struct classifiedPolygon
	{
		double a, b, c, d;
		unsigned long int id;
		int dy;
		shared_ptr<classifiedPolygon> next;
	};

	struct classifiedEdge
	{
		double x;
		double dx;
		int dy;
		unsigned long int id;

		Color color;
		Color dColor_x;
		Color dColor_y;
		
		shared_ptr<classifiedEdge> next;
	};

	struct activePolygon
	{
		double a, b, c, d;
		unsigned long int id;
		int dy;
		shared_ptr<activePolygon> next;

		activePolygon(classifiedPolygon &cp) {
			this->a = cp.a;
			this->b = cp.b;
			this->c = cp.c;
			this->d = cp.d;
			this->id = cp.id;
			this->dy = cp.dy;
		}
	};

	struct activeEdge
	{
		double x;
		double dx;
		int dy;
		double z;
		double dz_x;
		double dz_y;
		unsigned long int id;

		Color color;
		Color dColor_x;
		Color dColor_y;
		shared_ptr<activeEdge> next;

		activeEdge() {

		}

		activeEdge(classifiedEdge &cE, classifiedPolygon &cP, double y) {
			this->x = cE.x;
			this->dx = cE.dx;
			this->dy = cE.dy;

			this->z = cal_plane_z(cP.a, cP.b, cP.c, cP.d, this->x, y);
			this->dz_x = -cP.a / cP.c;
			this->dz_y = cP.b / cP.c;
			this->id = cP.id;

			this->color = cE.color;
			this->dColor_x = cE.dColor_x;
			this->dColor_y = cE.dColor_y;
		}

		activeEdge(classifiedEdge &cE, activePolygon &aP, double y) {
			this->x = cE.x;
			this->dx = cE.dx;
			this->dy = cE.dy;

			this->z = cal_plane_z(aP.a, aP.b, aP.c, aP.d, this->x, y);
			this->dz_x = -aP.a / aP.c;
			this->dz_y = aP.b / aP.c;
			this->id = aP.id;

			this->color = cE.color;
			this->dColor_x = cE.dColor_x;
			this->dColor_y = cE.dColor_y;
		}
	};

	struct inPolygon
	{
		//double a, b, c, d;
		unsigned long int id;

		double x;
		double z;
		double dz_x;
		Color color;
		Color dColor_x;
		shared_ptr<inPolygon> next;

		inPolygon() {

		}

		inPolygon(activePolygon &aP, activeEdge &aE) {
			//this->a = aP.a;
			//this->b = aP.b;
			//this->c = aP.c;
			//this->d = aP.d;
			this->id = aP.id;

			this->z = aE.z;
			this->x = aE.x;
			this->dz_x = aE.dz_x;
			this->color = aE.color;
			this->dColor_x = aE.dColor_x;
		}
	};

	struct doubleList {
		double val;
		shared_ptr<doubleList> next;

		doubleList(double val) {
			this->val = val;
		}
	};

	struct intList {
		int val;
		shared_ptr<intList> next;

		intList(int val) {
			this->val = val;
		}
	};

	struct colorList {
		Color color;
		shared_ptr<colorList> next;

		colorList(Color color) {
			this->color = color;
		}
	};

	template<class T> void insert_to_list(shared_ptr<T> *list, shared_ptr<T> &p, int position);

	// return the find num, if find pair pos 0 is left, pos 1 is right
	int find_edge_pair(shared_ptr<classifiedEdge> *ptr_cE[2], shared_ptr<classifiedEdge> *cE, int yid, int pid);

	void update_active_edge(activeEdge &aE);

	int find_aP(shared_ptr<activePolygon> *aP_list, shared_ptr<activePolygon> &aP, int pid);

	void update_IPL(shared_ptr<inPolygon> &IPL, shared_ptr<activePolygon> &aP, shared_ptr<activeEdge> &aE);

	void sortActiveEdge(shared_ptr<activeEdge> &head);
}

using namespace area_scanner;

class areaScannerZbuffer
{
public:
	areaScannerZbuffer(float *proj_coord, float *vColor, unsigned long int coord_num, Face *face_list, unsigned long int face_num, int wndWidth, int wndHeight);
	~areaScannerZbuffer();

	void set_cam_z(float nearZ, float farZ);
	void init();
	void scan_to_image(unsigned char *frame_buffer);

private:
	shared_ptr<classifiedPolygon> *cP;
	shared_ptr<classifiedEdge> *cE;
	shared_ptr<activePolygon> *aP;
	shared_ptr<activeEdge> *aE;
	Face *face_list;
	float *proj_coord;
	float *vColor;

	int wndWidth;
	int wndHeight;
	unsigned long int coord_num;
	unsigned long int face_num;

	double cam_nearZ;
	double cam_farZ;

	void gen_classified_data();
	void draw_line_one_color(unsigned char* frame, int x_left, int x_right, Color color);
	void draw_lines(unsigned char* frame, shared_ptr<intList> &inter_x, shared_ptr<colorList> &color, shared_ptr<colorList> &dColor);
	int get_color_dcolor(shared_ptr<inPolygon> &IPL, int x_left, int x_right, shared_ptr<colorList> &color, shared_ptr<colorList> &dColor, shared_ptr<intList> &inter_x);
};


