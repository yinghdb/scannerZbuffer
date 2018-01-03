#include "scannerZbuffer.h"

scannerZbuffer::scannerZbuffer(float *proj_coord, float *vColor, unsigned long int coord_num, Face *face_list, unsigned long int face_num, int wndWidth, int wndHeight)
	:face_list(face_list), proj_coord(proj_coord), wndWidth(wndWidth), wndHeight(wndHeight), coord_num(coord_num), face_num(face_num), vColor(vColor) {

	cP = new shared_ptr<classifiedPolygon>[wndHeight];
	cE = new shared_ptr<classifiedEdge>[wndHeight];
	aP = new shared_ptr<activePolygon>();
	aE = new shared_ptr<activeEdge>();
}

scannerZbuffer::~scannerZbuffer()
{
	delete[] cP;
	delete[] cE;
	delete aP;
	delete aE;
}

void scannerZbuffer::set_cam_z(float nearZ, float farZ) {
	this->cam_nearZ = nearZ;
	this->cam_farZ = farZ;
}

void scannerZbuffer::init() {
	gen_classified_data();
}

void scannerZbuffer::gen_classified_data() {
	unsigned long int fid;
	unsigned long int pid = 0;

	for (fid = 0; fid < face_num; ++fid) {
		Face *face = &face_list[fid];
		float vertex[3][3];
		float color[3][4];
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
			for (int j = 0; j < 4; ++j) {
				color[i][j] = vColor[12 * fid + 4 * i + j];
			}
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
		//if (AB.x*AC.y == AB.y*AC.x)
		//	continue;

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
			color_a.r = color[eid][0];
			color_a.g = color[eid][1];
			color_a.b = color[eid][2];
			color_a.a = color[eid][3];
			color_b.r = color[(eid + 1) % 3][0];
			color_b.g = color[(eid + 1) % 3][1];
			color_b.b = color[(eid + 1) % 3][2];
			color_b.a = color[(eid + 1) % 3][3];

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
			int dy = maxY_edge - minY_edge;
			dColor.r = (botColor.r - topColor.r) / dy;
			dColor.g = (botColor.g - topColor.g) / dy;
			dColor.b = (botColor.b - topColor.b) / dy;
			dColor.a = (botColor.a - topColor.a) / dy;
			if (maxY_edge >= wndHeight - 1) {
				int cy = maxY_edge - wndHeight + 1;
				topColor.r = topColor.r + dColor.r * cy;
				topColor.g = topColor.g + dColor.g * cy;
				topColor.b = topColor.b + dColor.b * cy;
				topColor.a = topColor.a + dColor.a * cy;
				maxY_edge = wndHeight - 1;
			}
			if (maxY_edge != maxY) {
				topColor.r = topColor.r + dColor.r;
				topColor.g = topColor.g + dColor.g;
				topColor.b = topColor.b + dColor.b;
				topColor.a = topColor.a + dColor.a;
				maxY_edge -= 1;
			}
			if (minY_edge < 0) {
				minY_edge = 0;
			}
			shared_ptr<classifiedEdge> p_cE = make_shared<classifiedEdge>();
			p_cE->id = pid;										//id
			p_cE->x = intersect(a.x, a.y, b.x, b.y, maxY_edge);	//x
			p_cE->dx = - (a.x - b.x) / (a.y - b.y);				//dx 
			p_cE->dy = maxY_edge - minY_edge;					//dy
			p_cE->fColor[0] = topColor.r;
			p_cE->fColor[1] = topColor.g;
			p_cE->fColor[2] = topColor.b;
			p_cE->fColor[3] = topColor.a;
			p_cE->dfColor[0] = dColor.r;
			p_cE->dfColor[1] = dColor.g;
			p_cE->dfColor[2] = dColor.b;
			p_cE->dfColor[3] = dColor.a;

			insert_to_list<classifiedEdge>(cE, p_cE, maxY_edge);
		}

		++pid;
	}
}



void scannerZbuffer::scan_to_image(unsigned char *frame_buffer) {
	int yid;
	double *zbuffer = new double[wndWidth];

	//// debug: 判断同一行的同一个多边形的边是否为两个
	//for (yid = wndHeight - 1; yid >= 0; --yid) {
	//	while (cP[yid]) {
	//		// add to active polygon table
	//		shared_ptr<activePolygon> p_aP = make_shared<activePolygon>(*cP[yid]);

	//		// check and add edge table
	//		shared_ptr<classifiedEdge> *ptr_cE[2];
	//		int find_num = find_edge_pair(ptr_cE, cE, yid, p_aP->id);
	//		shared_ptr<activeEdge> p_aE;
	//		if (find_num == 2) {
	//			//p_aE = make_shared<activeEdge>(*(*ptr_cE[0]), *(*ptr_cE[1]), *cP[yid], yid);
	//			//insert_to_list<activeEdge>(aE, p_aE, 0);
	//			//(*ptr_cE[1]) = (*ptr_cE[1])->next;
	//			//(*ptr_cE[0]) = (*ptr_cE[0])->next;
	//		}
	//		else if (find_num == -2) {
	//			//p_aE = make_shared<activeEdge>(*(*ptr_cE[1]), *(*ptr_cE[0]), *cP[yid], yid);
	//			//insert_to_list<activeEdge>(aE, p_aE, 0);
	//			//(*ptr_cE[1]) = (*ptr_cE[1])->next;
	//			//(*ptr_cE[0]) = (*ptr_cE[0])->next;
	//		}
	//		else if (find_num == 0){
	//			printf("error new polygon edge pair, yid: %d, pid: %d\n", yid, p_aP->id);
	//		}

	//		cP[yid] = cP[yid]->next;
	//	}
	//}
	for (yid = wndHeight - 1; yid >= 0; --yid) {
		// set frame buffer as background color, default black
		for (int i = 0; i < wndWidth; ++i) {
			for (int j = 0; j < 4; ++j) {
				frame_buffer[yid * wndWidth * 4 + 4 * i + j] = 0;
			}
		}
		// set z buffer as scale_farZ
		for (int i = 0; i < wndWidth; ++i) {
			zbuffer[i] = cam_nearZ;
		}

		// replace active polygon
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

		// replace active edge
		shared_ptr<activeEdge> *ptr_aE;
		ptr_aE = aE;
		while (*ptr_aE) {
			bool if_delete = false;
			// check if end
			if ((*ptr_aE)->dy_l < 0 && (*ptr_aE)->dy_r < 0) {
				shared_ptr<classifiedEdge> *ptr_cE[2];
				int find_num = find_edge_pair(ptr_cE, cE, yid, (*ptr_aE)->id);
				if (find_num == 2) {
					shared_ptr<activePolygon> p_aP;
					int res = scanner::find_aP(aP, p_aP, (*ptr_aE)->id);
					shared_ptr<activeEdge> p_aE = make_shared<activeEdge>(*(*ptr_cE[0]), *(*ptr_cE[1]), *p_aP, yid);
					insert_to_list<activeEdge>(aE, p_aE, 0);

					(*ptr_cE[0]) = (*ptr_cE[0])->next;
					(*ptr_cE[1]) = (*ptr_cE[1])->next;
				}
				else if (find_num == -2) {
					shared_ptr<activePolygon> p_aP;
					int res = scanner::find_aP(aP, p_aP, (*ptr_aE)->id);
					shared_ptr<activeEdge> p_aE = make_shared<activeEdge>(*(*ptr_cE[1]), *(*ptr_cE[0]), *p_aP, yid);
					insert_to_list<activeEdge>(aE, p_aE, 0);

					(*ptr_cE[0]) = (*ptr_cE[0])->next;
					(*ptr_cE[1]) = (*ptr_cE[1])->next;
				}
				else if (find_num == 1) {
					printf("error old polygon new edge pair, yid: %d, pid: %d\n", yid, (*ptr_aE)->id);
				}

				//delete
				if_delete = true;
			}
			else if ((*ptr_aE)->dy_l < 0) {
				shared_ptr<classifiedEdge> *ptr_cE[2];
				int find_num = find_edge_pair(ptr_cE, cE, yid, (*ptr_aE)->id);
				if (find_num == 1) {
					shared_ptr<activePolygon> p_aP;
					int res = scanner::find_aP(aP, p_aP, (*ptr_aE)->id);
					if (res == 1) {
						active_edge_replace_left(**ptr_cE[0], **ptr_aE, *p_aP, yid);
						(*ptr_cE[0]) = (*ptr_cE[0])->next;
					}
					else {
						printf("error old polygon new edge left, yid: %d, pid: %d\n", yid, (*ptr_aE)->id);
					}
				}
				else {
					printf("error old polygon new edge left, yid: %d, pid: %d\n", yid, (*ptr_aE)->id);
				}

			}
			else if ((*ptr_aE)->dy_r < 0) {
				shared_ptr<classifiedEdge> *ptr_cE[2];
				int find_num = find_edge_pair(ptr_cE, cE, yid, (*ptr_aE)->id);
				if (find_num == 1) {
					active_edge_replace_right(**ptr_cE[0], **ptr_aE);
					(*ptr_cE[0]) = (*ptr_cE[0])->next;
				}
				else {
					printf("error old polygon new edge right, yid: %d, pid: %d\n", yid, (*ptr_aE)->id);
				}

			}

			if (if_delete) {
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

			// check and add edge table
			shared_ptr<classifiedEdge> *ptr_cE[2];
			int find_num = find_edge_pair(ptr_cE, cE, yid, p_aP->id);
			shared_ptr<activeEdge> p_aE;
			if (find_num == 2) {
				p_aE = make_shared<activeEdge>(*(*ptr_cE[0]), *(*ptr_cE[1]), *cP[yid], yid);
				insert_to_list<activeEdge>(aE, p_aE, 0);
				(*ptr_cE[1]) = (*ptr_cE[1])->next;
				(*ptr_cE[0]) = (*ptr_cE[0])->next;
			}
			else if (find_num == -2) {
				p_aE = make_shared<activeEdge>(*(*ptr_cE[1]), *(*ptr_cE[0]), *cP[yid], yid);
				insert_to_list<activeEdge>(aE, p_aE, 0);
				(*ptr_cE[1]) = (*ptr_cE[1])->next;
				(*ptr_cE[0]) = (*ptr_cE[0])->next;
			}
			else if (find_num != 0) {
				printf("error new polygon edge pair, yid: %d, pid: %d\n", yid, p_aP->id);
			}
			else {
				printf("same edge pair, yid: %d, pid: %d\n", yid, p_aP->id);
			}

			cP[yid] = cP[yid]->next;
		}

		// update zbuffer
		auto p_aE = *aE;
		while (p_aE) {
			int x_left = ceil(p_aE->x_l);
			int x_right = floor(p_aE->x_r);
			double zx = p_aE->z_l + p_aE->dz_x * (x_left - p_aE->x_l);
			float *color_left = p_aE->fColor_l;
			float *color_right = p_aE->fColor_r;
			float dColor[4];
			float color[4];
			int dx = x_right - x_left;
			dColor[0] = (color_right[0] - color_left[0]) / dx;
			dColor[1] = (color_right[1] - color_left[1]) / dx;
			dColor[2] = (color_right[2] - color_left[2]) / dx;
			dColor[3] = (color_right[3] - color_left[3]) / dx;
			color[0] = color_left[0];
			color[1] = color_left[1];
			color[2] = color_left[2];
			color[3] = color_left[3];
			for (int x = x_left; x <= x_right; ++x) {
				// compare with the zbuffer
				double z_buf = zbuffer[x];
				if (zx > z_buf) {
					zbuffer[x] = zx;
					// temp white color
					for (int j = 0; j < 4; ++j) {
						if (color[j] < 0)
							frame_buffer[yid * wndWidth * 4 + 4 * x + j] = 0;
						else if (color[j] > 1)
							frame_buffer[yid * wndWidth * 4 + 4 * x + j] = 255;
						else
							frame_buffer[yid * wndWidth * 4 + 4 * x + j] = (int)(color[j] * 255);
					}
				}
				zx += p_aE->dz_x;
				color[0] += dColor[0];
				color[1] += dColor[1];
				color[2] += dColor[2];
				color[3] += dColor[3];
			}
			p_aE = p_aE->next;
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

	delete[] zbuffer;
}

template<class T> void scanner::insert_to_list(shared_ptr<T> *list, shared_ptr<T> &p, int position) {
	p->next = list[position];
	list[position] = p;
}
// return the find num,
// when find 2 edge, return positive if the ptr_cE order is right, else return negetive
int scanner::find_edge_pair(shared_ptr<classifiedEdge> *ptr_cE[2], shared_ptr<classifiedEdge> *cE, int yid, int pid) {
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
		else if ((*ptr_cE[0])->dx > (*ptr_cE[1])->dx) {
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

void scanner::active_edge_replace_left(classifiedEdge &cE, activeEdge &aE, activePolygon &aP, int yid) {
	aE.x_l = cE.x;
	aE.dx_l = cE.dx;
	aE.dy_l = cE.dy;
	aE.z_l = cal_plane_z(aP.a, aP.b, aP.c, aP.d, aE.x_l, yid);

	aE.fColor_l[0] = cE.fColor[0];
	aE.fColor_l[1] = cE.fColor[1];
	aE.fColor_l[2] = cE.fColor[2];
	aE.fColor_l[3] = cE.fColor[3];

	aE.dfColor_l[0] = cE.dfColor[0];
	aE.dfColor_l[1] = cE.dfColor[1];
	aE.dfColor_l[2] = cE.dfColor[2];
	aE.dfColor_l[3] = cE.dfColor[3];
}

void scanner::active_edge_replace_right(classifiedEdge &cE, activeEdge &aE) {
	aE.x_r = cE.x;
	aE.dx_r = cE.dx;
	aE.dy_r = cE.dy;

	aE.fColor_r[0] = cE.fColor[0];
	aE.fColor_r[1] = cE.fColor[1];
	aE.fColor_r[2] = cE.fColor[2];
	aE.fColor_r[3] = cE.fColor[3];

	aE.dfColor_r[0] = cE.dfColor[0];
	aE.dfColor_r[1] = cE.dfColor[1];
	aE.dfColor_r[2] = cE.dfColor[2];
	aE.dfColor_r[3] = cE.dfColor[3];
}

void scanner::update_active_edge(activeEdge &aE) {
	aE.x_l += aE.dx_l;
	aE.x_r += aE.dx_r;
	aE.z_l += aE.dz_x*aE.dx_l + aE.dz_y;

	aE.dy_l -= 1;
	aE.dy_r -= 1;

	aE.fColor_l[0] += aE.dfColor_l[0];
	aE.fColor_l[1] += aE.dfColor_l[1];
	aE.fColor_l[2] += aE.dfColor_l[2];
	aE.fColor_l[3] += aE.dfColor_l[3];

	aE.fColor_r[0] += aE.dfColor_r[0];
	aE.fColor_r[1] += aE.dfColor_r[1];
	aE.fColor_r[2] += aE.dfColor_r[2];
	aE.fColor_r[3] += aE.dfColor_r[3];
}

int scanner::find_cP(shared_ptr<classifiedPolygon> &cP_list, shared_ptr<classifiedPolygon> &cP, int pid) {
	cP = make_shared<classifiedPolygon>(*cP_list);
	while (cP) {
		if ((cP)->id == pid)
			return 1;
		cP = (cP)->next;
	}
	return -1;
}

int scanner::find_aP(shared_ptr<activePolygon> *aP_list, shared_ptr<activePolygon> &aP, int pid) {
	aP = *aP_list;
	while (aP) {
		if ((aP)->id == pid)
			return 1;
		aP = (aP)->next;
	}
	return -1;
}