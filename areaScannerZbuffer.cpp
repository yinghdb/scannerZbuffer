#include "areaScannerZbuffer.h"

areaScannerZbuffer::areaScannerZbuffer(float *proj_coord, float *vColor, unsigned long int coord_num,
	Face *face_list, unsigned long int face_num, int wndWidth, int wndHeight)
	:face_list(face_list), proj_coord(proj_coord), wndWidth(wndWidth), wndHeight(wndHeight), coord_num(coord_num), face_num(face_num), vColor(vColor) {

	cP = new shared_ptr<classifiedPolygon>[wndHeight];
	cE = new shared_ptr<classifiedEdge>[wndHeight];
	aP = new shared_ptr<activePolygon>();
	aE = new shared_ptr<activeEdge>();
}

areaScannerZbuffer::~areaScannerZbuffer()
{
	delete[] cP;
	delete[] cE;
	delete aP;
	delete aE;
}

void areaScannerZbuffer::set_cam_z(float nearZ, float farZ) {
	this->cam_nearZ = nearZ;
	this->cam_farZ = farZ;
}

void areaScannerZbuffer::init() {
	gen_classified_data();
}

void areaScannerZbuffer::gen_classified_data() {
	unsigned long int fid;
	unsigned long int pid = 0;

	for (fid = 0; fid < face_num; ++fid) {
		Face *face = &face_list[fid];
		float vertex[3][3];
		Color color[3];
		float minmax[3][2];

		for (int i = 0; i < 3; ++i) {
			if (i == 0) {
				for (int j = 0; j < 3; j++) {
					float val = proj_coord[3 * face->vIndices[i] + j];
					vertex[i][j] = val;
					minmax[j][0] = val;
					minmax[j][1] = val;
				}
			}
			else {
				for (int j = 0; j < 3; j++) {
					float val = proj_coord[3 * face->vIndices[i] + j];
					vertex[i][j] = val;
					if (val < minmax[j][0]) {
						minmax[j][0] = val;
					}
					else if (val > minmax[j][1]) {
						minmax[j][1] = val;
					}
				}
			}
		}

		for (int i = 0; i < 3; ++i) {
			color[i].r = vColor[12 * fid + 4 * i + 0];
			color[i].g = vColor[12 * fid + 4 * i + 1];
			color[i].b = vColor[12 * fid + 4 * i + 2];
			color[i].a = vColor[12 * fid + 4 * i + 3];
		}

		//判断是否在界内
		bool notInBound = false;
		for (int i = 0; i < 3; ++i) {
			if (minmax[i][0] > 1 || minmax[i][1] < -1) {
				notInBound = true;
				break;
			}
		}
		if (notInBound)
			continue;

		// 将坐标变换到整数
		int vertex_int[3][3];
		for (int i = 0; i < 3; ++i) {
			vertex_int[i][0] = round((vertex[i][0] + 1) * (wndWidth - 1) / 2);
			vertex_int[i][1] = round((vertex[i][1] + 1) * (wndHeight - 1) / 2);
			vertex_int[i][2] = round((vertex[i][2] + 1) * (wndWidth - 1) / 2);
		}

		//生成分类多边形表
		Vector3d AB, AC, norm_ABC, A;
		//AB.x = (vertex[1][0] - vertex[0][0]) * (wndWidth - 1) / 2;
		//AB.y = (vertex[1][1] - vertex[0][1]) * (wndHeight - 1) / 2;
		//AB.z = (vertex[1][2] - vertex[0][2]) * (wndWidth - 1) / 2;
		//AC.x = (vertex[2][0] - vertex[0][0]) * (wndWidth - 1) / 2;
		//AC.y = (vertex[2][1] - vertex[0][1]) * (wndHeight - 1) / 2;
		//AC.z = (vertex[2][2] - vertex[0][2]) * (wndWidth - 1) / 2;
		//A.x = (vertex[0][0] + 1) * (wndWidth - 1) / 2;
		//A.y = (vertex[0][1] + 1) * (wndHeight - 1) / 2;
		//A.z = (vertex[0][2] + 1) * (wndWidth - 1) / 2;
		AB.x = (vertex_int[1][0] - vertex_int[0][0]);
		AB.y = (vertex_int[1][1] - vertex_int[0][1]);
		AB.z = (vertex_int[1][2] - vertex_int[0][2]);
		AC.x = (vertex_int[2][0] - vertex_int[0][0]);
		AC.y = (vertex_int[2][1] - vertex_int[0][1]);
		AC.z = (vertex_int[2][2] - vertex_int[0][2]);

		A.x = vertex_int[0][0];
		A.y = vertex_int[0][1];
		A.z = vertex_int[0][2];
		cross_product(AB, AC, norm_ABC);
		normalize(norm_ABC);

		if (norm_ABC.z == 0 || isnan(norm_ABC.z))	// 平面垂直于视平面
			continue;
		if (AB.x*AC.y == AB.y*AC.x)
			continue;

		int minY = round((minmax[1][0] + 1) * (wndHeight - 1) / 2);
		int maxY = round((minmax[1][1] + 1) * (wndHeight - 1) / 2);
		if (minY < 0)
			minY = 0;
		if (maxY > wndHeight - 1)
			maxY = wndHeight - 1;
		if (maxY <= minY)		// 平行于x轴
			continue;
		shared_ptr<classifiedPolygon> p_cP = make_shared<classifiedPolygon>();
		p_cP->id = pid;							//id
		p_cP->dy = maxY - minY;					//dy
		p_cP->a = norm_ABC.x;
		p_cP->b = norm_ABC.y;
		p_cP->c = norm_ABC.z;
		p_cP->d = -dot_product(norm_ABC, A);	//a,b,c,d

		insert_to_list<classifiedPolygon>(cP, p_cP, maxY);

		Color dColor_x, dColor_y;
		Color dColor_AB, dColor_AC;
		color_sub(color[1], color[0], dColor_AB);
		color_sub(color[2], color[0], dColor_AC);
		double xy_sub = 1 / (AB.x*AC.y - AB.y*AC.x);

		Color x1c2, x2c1;
		Color y2c1, y1c2;

		color_scale(dColor_AB, AC.y, y2c1);
		color_scale(dColor_AC, AB.y, y1c2);
		color_scale(dColor_AC, AB.x, x1c2);
		color_scale(dColor_AB, AC.x, x2c1);

		color_sub(y2c1, y1c2, dColor_x);
		color_scale(dColor_x, xy_sub, dColor_x);
		color_sub(x1c2, x2c1, dColor_y);
		color_scale(dColor_y, xy_sub, dColor_y);

		//生成分类边表
		for (int eid = 0; eid < 3; ++eid) {
			Vector3d a, b;
			Color color_a, color_b;
			a.x = vertex_int[eid][0];
			a.y = vertex_int[eid][1];
			a.z = vertex_int[eid][2];
			b.x = vertex_int[(eid + 1) % 3][0];
			b.y = vertex_int[(eid + 1) % 3][1];
			b.z = vertex_int[(eid + 1) % 3][2];
			color_a = color[eid];
			color_b = color[(eid + 1) % 3];

			//if ((a.x < -1 && b.x < -1) || (a.x > 1 && b.x > 1)
			//	|| (a.y < -1 && b.y < -1) || (a.y > 1 && b.y > 1)
			//	|| (a.z < nearZ && b.z < nearZ) || (a.z > farZ && b.z > farZ))
			//	continue;

			int topY, botY;
			Color topColor, botColor, dColor;
			if (a.y > b.y) {
				topY = a.y;
				botY = b.y;
				topColor = color_a;
				botColor = color_b;
			}
			else {
				topY = b.y;
				botY = a.y;
				topColor = color_b;
				botColor = color_a;
			}
			int minY_edge = botY;
			int maxY_edge = topY;
			if (maxY_edge - minY_edge <= 0)
				continue;
			double dy = maxY_edge - minY_edge;
			color_sub(botColor, topColor, dColor);
			color_scale(dColor, 1 / dy, dColor);
			if (maxY_edge >= wndHeight - 1) {
				int cy = maxY_edge - wndHeight + 1;
				Color sfColor;
				color_scale(dColor, cy, sfColor);
				color_add(topColor, sfColor, topColor);
				maxY_edge = wndHeight - 1;
			}
			if (maxY_edge != maxY) {
				color_add(topColor, dColor, topColor);
				maxY_edge -= 1;
			}
			if (minY_edge < 0) {
				minY_edge = 0;
			}
			shared_ptr<classifiedEdge> p_cE = make_shared<classifiedEdge>();
			p_cE->id = pid;										//id
			p_cE->x = intersect(a.x, a.y, b.x, b.y, maxY_edge);	//x
			p_cE->dx = -(a.x - b.x) / (a.y - b.y);				//dx 
			p_cE->dy = maxY_edge - minY_edge;					//dy
			p_cE->color = topColor;
			p_cE->dColor_x = dColor_x;
			p_cE->dColor_y = dColor_y;

			insert_to_list<classifiedEdge>(cE, p_cE, maxY_edge);
		}

		++pid;
	}
}



void areaScannerZbuffer::scan_to_image(unsigned char *frame_buffer) {
	int yid;
	Color background_color;
	background_color.r = 0.0;
	background_color.g = 0.0;
	background_color.b = 0.0;
	background_color.a = 1.0;

	for (yid = wndHeight - 1; yid >= 0; --yid) {
		// delete active polygon
		shared_ptr<activePolygon> *ptr_aP;
		ptr_aP = aP;
		while (*ptr_aP) {
			if ((*ptr_aP)->dy < 0) {
				// delete polygon
				(*ptr_aP) = (*ptr_aP)->next;
			}
			else {
				ptr_aP = &(*ptr_aP)->next;
			}
		}

		// delete active edge
		shared_ptr<activeEdge> *ptr_aE;
		ptr_aE = aE;
		while (*ptr_aE) {
			if ((*ptr_aE)->dy < 0) {
				// delete
				(*ptr_aE) = (*ptr_aE)->next;
			}
			else {
				ptr_aE = &(*ptr_aE)->next;
			}
		}

		// check classified polygon table to add
		while (cP[yid]) {
			// add to active polygon table
			shared_ptr<activePolygon> p_aP = make_shared<activePolygon>(*cP[yid]);
			insert_to_list<activePolygon>(aP, p_aP, 0);

			cP[yid] = cP[yid]->next;
		}

		// check classified edge table to add
		while (cE[yid]) {
			shared_ptr<activePolygon> p_aP;
			int res = area_scanner::find_aP(aP, p_aP, cE[yid]->id);
			shared_ptr<activeEdge> p_aE = make_shared<activeEdge>(*cE[yid], *p_aP, yid);
			area_scanner::insert_to_list<activeEdge>(aE, p_aE, 0);

			cE[yid] = cE[yid]->next;
		}

		// sort activeEdge table
		area_scanner::sortActiveEdge(*aE);

		// update frame buffer
		if (!(*aE)) {
			// draw background color
			draw_line_one_color(frame_buffer + yid * wndWidth * 4, 0, wndWidth, background_color);
		}
		else {
			shared_ptr<colorList> color;
			shared_ptr<colorList> dColor;
			shared_ptr<intList> inter_x;
			shared_ptr<activeEdge> p_aE = *aE;
			shared_ptr<activeEdge> p_aEnext = p_aE->next;
			shared_ptr<inPolygon> IPL;
			draw_line_one_color(frame_buffer + yid * wndWidth * 4, 0, p_aE->x, background_color);
			while (p_aEnext) {
				int pid = p_aE->id;
				shared_ptr<activePolygon> p_aP;
				find_aP(aP, p_aP, pid);
				area_scanner::update_IPL(IPL, p_aP, p_aE);

				if (p_aE->x == p_aEnext->x) {
					p_aE = p_aE->next;
					p_aEnext = p_aEnext->next;
					continue;
				}
				int res = get_color_dcolor(IPL, p_aE->x, p_aEnext->x, color, dColor, inter_x);
				if (res != -1)
					draw_lines(frame_buffer + yid * wndWidth * 4, inter_x, color, dColor);
				else
					draw_line_one_color(frame_buffer + yid * wndWidth * 4, p_aE->x, p_aEnext->x, background_color);

				p_aE = p_aE->next;
				p_aEnext = p_aEnext->next;
			}
			draw_line_one_color(frame_buffer + yid * wndWidth * 4, p_aE->x, wndWidth, background_color);
			if ((IPL)->next) {
				printf("error IPL end, yid: %d\n", yid);
			}
		}

		// update
		ptr_aP = aP;
		while (*ptr_aP) {
			(*ptr_aP)->dy -= 1;
			ptr_aP = &(*ptr_aP)->next;
		}
		ptr_aE = aE;
		while (*ptr_aE) {
			update_active_edge(**ptr_aE);
			ptr_aE = &(*ptr_aE)->next;
		}

	}
}

void areaScannerZbuffer::draw_line_one_color(unsigned char* frame, int x_left, int x_right, Color color) {
	for (int x = x_left; x < x_right; ++x) {
		if (color.r > 1)
			frame[4 * x + 0] = 255;
		else if(color.r < 0)
			frame[4 * x + 0] = 0;
		else
			frame[4 * x + 0] = color.r * 255;

		if (color.g > 1)
			frame[4 * x + 1] = 255;
		else if (color.g < 0)
			frame[4 * x + 1] = 0;
		else
			frame[4 * x + 1] = color.g * 255;

		if (color.b > 1)
			frame[4 * x + 2] = 255;
		else if (color.b < 0)
			frame[4 * x + 2] = 0;
		else
			frame[4 * x + 2] = color.b * 255;

		if (color.a > 1)
			frame[4 * x + 3] = 255;
		else if (color.a < 0)
			frame[4 * x + 3] = 0;
		else
			frame[4 * x + 3] = color.a * 255;
	}
}

void areaScannerZbuffer::draw_lines(unsigned char* frame, shared_ptr<intList> &inter_x, shared_ptr<colorList> &color, shared_ptr<colorList> &dColor) {
	shared_ptr<intList> p_inter_x = inter_x;
	shared_ptr<colorList> p_color = color;
	shared_ptr<colorList> p_dColor = dColor;
	while (p_inter_x->next) {
		Color m_color = p_color->color;
		Color d_color = p_dColor->color;
		int x_left = p_inter_x->val;
		int x_right = p_inter_x->next->val;

		for (int x = x_left; x < x_right; ++x) {
			if (m_color.r > 1)
				frame[4 * x + 0] = 255;
			else if (m_color.r < 0)
				frame[4 * x + 0] = 0;
			else
				frame[4 * x + 0] = m_color.r * 255;

			if (m_color.g > 1)
				frame[4 * x + 1] = 255;
			else if (m_color.g < 0)
				frame[4 * x + 1] = 0;
			else
				frame[4 * x + 1] = m_color.g * 255;

			if (m_color.b > 1)
				frame[4 * x + 2] = 255;
			else if (m_color.b < 0)
				frame[4 * x + 2] = 0;
			else
				frame[4 * x + 2] = m_color.b * 255;

			if (m_color.a > 1)
				frame[4 * x + 3] = 255;
			else if (m_color.a < 0)
				frame[4 * x + 3] = 0;
			else
				frame[4 * x + 3] = m_color.a * 255;

			color_add(m_color, d_color, m_color);
		}

		p_inter_x = p_inter_x->next;
		p_color = p_color->next;
		p_dColor = p_dColor->next;
	}
}

int areaScannerZbuffer::get_color_dcolor(shared_ptr<inPolygon> &IPL, int x_left, int x_right, shared_ptr<colorList> &color, shared_ptr<colorList> &dColor, shared_ptr<intList> &inter_x) {
	if (!IPL) {
		return -1;
	}
	
	shared_ptr<doubleList> inter_z;
	shared_ptr<inPolygon> pIPL = IPL;


	shared_ptr<intList> xl = make_shared<intList>(x_left);
	inter_x = xl;
	shared_ptr<intList> xr = make_shared<intList>(x_right);
	inter_x->next = xr;

	double z_l, z_r;
	z_l = IPL->z + IPL->dz_x*(x_left - IPL->x);
	z_r = IPL->z + IPL->dz_x*(x_right - IPL->x);
	shared_ptr<doubleList> zl = make_shared<doubleList>(z_l);
	inter_z = zl;
	shared_ptr<doubleList> zr = make_shared<doubleList>(z_r);
	inter_z->next = zr;

	Color color_l;
	color_scale(IPL->dColor_x, x_left - IPL->x, color_l);
	color_add(IPL->color, color_l, color_l);
	color = make_shared<colorList>(color_l);

	dColor = make_shared<colorList>(IPL->dColor_x);

	while (pIPL) {
		shared_ptr<doubleList> inter_z_now;
		shared_ptr<doubleList> inter_z_now_end;
		shared_ptr<intList> p_inter_x = inter_x;
		while (p_inter_x) {
			int x = p_inter_x->val;
			double z = pIPL->z + pIPL->dz_x*(x - pIPL->x);
			shared_ptr<doubleList> zl = make_shared<doubleList>(z);
			if (inter_z_now) {
				inter_z_now_end->next = zl;
				inter_z_now_end = zl;
			}
			else {
				inter_z_now = zl;
				inter_z_now_end = zl;
			}

			p_inter_x = p_inter_x->next;
		}

		shared_ptr<doubleList> p_inter_z = inter_z;
		shared_ptr<colorList> p_color = color;
		shared_ptr<colorList> p_dColor = dColor;
		p_inter_x = inter_x;
		while (p_inter_z->next && inter_z_now->next && p_inter_x->next) {
			double dz_l = inter_z_now->val - p_inter_z->val;
			double dz_r = inter_z_now->next->val - p_inter_z->next->val;

			//if (dz_l > 0 && dz_r > 0) { // edge now is more near
			//	p_inter_z->val = inter_z_now->val;
			//	p_inter_z->next->val = inter_z_now->next->val;

			//	color_scale(pIPL->dColor_x, p_inter_x->val - pIPL->x, color_l);
			//	color_add(pIPL->color, color_l, color_l);
			//	p_color->color = color_l;
			//	p_dColor->color = pIPL->dColor_x;
			//}
			//else if (dz_l <= 0 && dz_r <= 0) { // edge now is more far

			//}

			if (dz_l + dz_r > 0) { // edge now is more near, 不考虑穿透
				p_inter_z->val = inter_z_now->val;
				p_inter_z->next->val = inter_z_now->next->val;

				color_scale(pIPL->dColor_x, p_inter_x->val - pIPL->x, color_l);
				color_add(pIPL->color, color_l, color_l);
				p_color->color = color_l;
				p_dColor->color = pIPL->dColor_x;
			}

			p_inter_z = p_inter_z->next;
			inter_z_now = inter_z_now->next;
			p_inter_x = p_inter_x->next;
			p_color = p_color->next;
			p_dColor = p_dColor->next;
		}

		pIPL = pIPL->next;
	}

	return 0;
}


template<class T> void area_scanner::insert_to_list(shared_ptr<T> *list, shared_ptr<T> &p, int position) {
	p->next = list[position];
	list[position] = p;
}
// return the find num,
// when find 2 edge, return positive if the ptr_cE order is right, else return negetive
int area_scanner::find_edge_pair(shared_ptr<classifiedEdge> *ptr_cE[2], shared_ptr<classifiedEdge> *cE, int yid, int pid) {
	shared_ptr<classifiedEdge> *ptr = &(cE[yid]);
	int p_num = 0;
	while (*ptr) {
		if ((*ptr)->id == pid) {
			ptr_cE[p_num] = ptr;
			++p_num;
		}
		if (p_num == 2) {
			break;
		}
		ptr = &(*ptr)->next;
	}
	if (p_num != 2) {
		return p_num;
	}
	if (abs((*ptr_cE[0])->x - (*ptr_cE[1])->x) < 0.01) {
		if ((*ptr_cE[0])->dx < (*ptr_cE[1])->dx) {
			return 2;
		}
		else if ((*ptr_cE[0])->dx >(*ptr_cE[1])->dx) {
			return -2;
		}
		else {
			return 0;
		}
	}
	else if ((*ptr_cE[0])->x < (*ptr_cE[1])->x) {
		return 2;
	}
	else {
		return -2;
	}
}


void area_scanner::update_active_edge(activeEdge &aE) {
	aE.x += aE.dx;
	aE.z += aE.dz_x*aE.dx + aE.dz_y;

	aE.dy -= 1;

	Color dColor_x;
	color_scale(aE.dColor_x, aE.dx, dColor_x);
	color_add(aE.color, dColor_x, aE.color);
	color_sub(aE.color, aE.dColor_y, aE.color);
}

int area_scanner::find_aP(shared_ptr<activePolygon> *aP_list, shared_ptr<activePolygon> &aP, int pid) {
	aP = *aP_list;
	while (aP) {
		if ((aP)->id == pid)
			return 1;
		aP = (aP)->next;
	}
	return -1;
}

void area_scanner::sortActiveEdge(shared_ptr<activeEdge> &head) {
	if (!(head) || !(head)->next)
		return;

	shared_ptr<activeEdge> p = head->next;
	shared_ptr<activeEdge> pstart = make_shared<activeEdge>();
	shared_ptr<activeEdge> pend = head;
	pstart->next = head; //为了操作方便，添加一个头结点
	while (p)
	{
		shared_ptr<activeEdge> tmp = pstart->next;
		shared_ptr<activeEdge> pre = pstart;
		while (tmp != p && p->x >= tmp->x) //找到插入位置
		{
			tmp = tmp->next;
			pre = pre->next;
		}
		if (tmp == p)
			pend = p;
		else
		{
			pend->next = p->next;
			p->next = tmp;
			pre->next = p;
		}
		p = pend->next;
	}
	head = pstart->next;
}

void area_scanner::update_IPL(shared_ptr<inPolygon> &IPL, shared_ptr<activePolygon> &aP, shared_ptr<activeEdge> &aE) {
	shared_ptr<inPolygon> p_start = make_shared<inPolygon>();
	p_start->next = IPL;
	shared_ptr<inPolygon> p_pre = p_start;

	shared_ptr<inPolygon> p = IPL;
	while (p) {
		if (p->id == aP->id) {
			p_pre->next = p->next;
			IPL = p_start->next;

			return;
		}
		p_pre = p_pre->next;
		p = p->next;
	}

	shared_ptr<inPolygon> pIPL = make_shared<inPolygon>(*aP, *aE);
	pIPL->next = IPL;
	IPL = pIPL;
}

