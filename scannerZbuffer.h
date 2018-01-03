#pragma once
#include "objParser.h"
#include "utils.h"
#include <memory>

using namespace std;

namespace scanner {
	struct classifiedPolygon
	{
		double a, b, c, d;
		unsigned long int id;
		int dy;
		unsigned char color[4];
		shared_ptr<classifiedPolygon> next;
	};

	struct classifiedEdge
	{
		double x;
		double dx;
		int dy;
		unsigned long int id;

		float fColor[4];
		float dfColor[4];
		
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
		double x_l;
		double dx_l;
		double dy_l;
		double x_r;
		double dx_r;
		double dy_r;
		double z_l;
		double dz_x;
		double dz_y;
		unsigned long int id;

		float fColor_l[4];
		float dfColor_l[4];
		float fColor_r[4];
		float dfColor_r[4];
		shared_ptr<activeEdge> next;

		activeEdge(classifiedEdge &cE_left, classifiedEdge &cE_right, classifiedPolygon &cP, double y) {
			this->x_l = cE_left.x;
			this->dx_l = cE_left.dx;
			this->dy_l = cE_left.dy;
			this->x_r = cE_right.x;
			this->dx_r = cE_right.dx;
			this->dy_r = cE_right.dy;

			this->z_l = cal_plane_z(cP.a, cP.b, cP.c, cP.d, this->x_l, y);
			this->dz_x = -cP.a / cP.c;
			this->dz_y = cP.b / cP.c;
			this->id = cP.id;

			this->fColor_l[0] = cE_left.fColor[0];
			this->fColor_l[1] = cE_left.fColor[1];
			this->fColor_l[2] = cE_left.fColor[2];
			this->fColor_l[3] = cE_left.fColor[3];
			this->dfColor_l[0] = cE_left.dfColor[0];
			this->dfColor_l[1] = cE_left.dfColor[1];
			this->dfColor_l[2] = cE_left.dfColor[2];
			this->dfColor_l[3] = cE_left.dfColor[3];

			this->fColor_r[0] = cE_right.fColor[0];
			this->fColor_r[1] = cE_right.fColor[1];
			this->fColor_r[2] = cE_right.fColor[2];
			this->fColor_r[3] = cE_right.fColor[3];
			this->dfColor_r[0] = cE_right.dfColor[0];
			this->dfColor_r[1] = cE_right.dfColor[1];
			this->dfColor_r[2] = cE_right.dfColor[2];
			this->dfColor_r[3] = cE_right.dfColor[3];
		}

		activeEdge(classifiedEdge &cE_left, classifiedEdge &cE_right, activePolygon &aP, double y) {
			this->x_l = cE_left.x;
			this->dx_l = cE_left.dx;
			this->dy_l = cE_left.dy;
			this->x_r = cE_right.x;
			this->dx_r = cE_right.dx;
			this->dy_r = cE_right.dy;

			this->z_l = cal_plane_z(aP.a, aP.b, aP.c, aP.d, this->x_l, y);
			this->dz_x = -aP.a / aP.c;
			this->dz_y = aP.b / aP.c;
			this->id = aP.id;

			this->fColor_l[0] = cE_left.fColor[0];
			this->fColor_l[1] = cE_left.fColor[1];
			this->fColor_l[2] = cE_left.fColor[2];
			this->fColor_l[3] = cE_left.fColor[3];
			this->dfColor_l[0] = cE_left.dfColor[0];
			this->dfColor_l[1] = cE_left.dfColor[1];
			this->dfColor_l[2] = cE_left.dfColor[2];
			this->dfColor_l[3] = cE_left.dfColor[3];

			this->fColor_r[0] = cE_right.fColor[0];
			this->fColor_r[1] = cE_right.fColor[1];
			this->fColor_r[2] = cE_right.fColor[2];
			this->fColor_r[3] = cE_right.fColor[3];
			this->dfColor_r[0] = cE_right.dfColor[0];
			this->dfColor_r[1] = cE_right.dfColor[1];
			this->dfColor_r[2] = cE_right.dfColor[2];
			this->dfColor_r[3] = cE_right.dfColor[3];
		}
	};

	template<class T> void insert_to_list(shared_ptr<T> *list, shared_ptr<T> &p, int position);
	// return the find num, if find pair pos 0 is left, pos 1 is right
	int find_edge_pair(shared_ptr<classifiedEdge> *ptr_cE[2], shared_ptr<classifiedEdge> *cE, int yid, int pid);

	void active_edge_replace_left(classifiedEdge &cE, activeEdge &aE, activePolygon &aP, int yid);

	void active_edge_replace_right(classifiedEdge &cE, activeEdge &aE);

	void update_active_edge(activeEdge &aE);

	int find_cP(shared_ptr<classifiedPolygon> &cP_list, shared_ptr<classifiedPolygon> &cP, int pid);

	int find_aP(shared_ptr<activePolygon> *aP_list, shared_ptr<activePolygon> &aP, int pid);
}

using namespace scanner;

class scannerZbuffer
{
public:
	scannerZbuffer(float *proj_coord, float *vColor, unsigned long int coord_num, Face *face_list, unsigned long int face_num, int wndWidth, int wndHeight);
	~scannerZbuffer();

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
};


