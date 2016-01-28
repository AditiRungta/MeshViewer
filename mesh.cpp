#include <GL/glui.h>
#include <string.h>
#include <iostream>
#include<string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <math.h>

using namespace std;

const int MAX_CHARS_PER_LINE = 1024;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = " ";

struct Vertex
{
	GLfloat x, y, z;
	GLfloat vnormal_vec[3];
};
int v_count = 0;
Vertex *m_vertices;
struct TriangleFace
{
	GLint v_indices[3];
	GLfloat fnormal_vec[3];
};
int f_count = 0;
TriangleFace *m_faces;

//Bounding Box
Vertex center;
float radius, maxwidth;
float xmax, xmin, ymax, ymin, zmax, zmin;

struct HE_vert;
struct HE_face;

struct HE_edge {
	HE_vert* vert; // vertex at the end of the half-edge
	HE_edge* pair; // oppositely oriented half-edge
	HE_face* face; // the incident face
	HE_edge* prev; // previous half-edge around the face
	HE_edge* next; // next half-edge around the face
};
HE_edge *m_HE_edge;
int he_e_count;

struct HE_vert {
	float x, y, z; // the vertex coordinates
	HE_edge* edge; // one of the half-edges emanating from the vertex
	float vnx, vny, vnz; // vertex normal
};
HE_vert *m_HE_vert;
int he_v_count;

struct HE_face {
	HE_edge* edge; // one of the half-edges bordering the face
	float fnx, fny, fnz; // surface normal
	float area; // surface area
};
HE_face *m_HE_face;
int he_f_count;

int *edgepair;

float xy_aspect;
int   old_x, old_y;
float rotationX = 0.0, rotationY = 0.0, rotationZ = 0.0, rotationXZ=0.0;
float moveX = 0.0, moveY = 0.0, moveZ = 0.0;
float rotobjX = 0.0, rotobjY = 0.0, rotobjZ = 0.0;
float moveobjX = 0.0, moveobjY = 0.0, moveobjZ = 0.0;

// These are the live variables passed into GLUI
int bunny = 0, cap = 1, gargoyle = 0, knot = 0, eight = 0;
int   wireframe = 0;
int	  pointcloud = 0;
int   flat = 0;
int   smooth = 1;
int   obj_type = 0;
int   rmode_type = 3;
int   proj_type = 0;
int   segments = 8;
int   segments2 = 8;
int   light0_enabled = 1;
int   light1_enabled = 1;
float light0_intensity = 1.0;
float light1_intensity = 1.0;
int   main_window;
float scale = 1.0;
float zoom = 1.0;
int leftbutton = 0;
int rightbutton = 0;
int middlebutton = 0;
float fov = 135.0;
float Oleft = -2.5, Oright = 2.5, Obottom = -2.5, Otop = 2.5;
float scaleXYZ = 1.0;
int showaxes = 1, showbbox = 1, showplane = 1;

// Pointers to the windows and some of the controls in Glui
GLUI *cmd_line_glui = 0, *glui;
GLUI_Checkbox    *checkbox1, *checkbox2, *checkbox3;
GLUI_Spinner     *light0_spinner, *light1_spinner, *scale_spinner;
GLUI_Spinner     *moveX_spinner, *moveY_spinner, *moveZ_spinner, *rotX_spinner, *rotY_spinner, *rotZ_spinner;
GLUI_RadioGroup  *radio, *radio1, *radio2;
GLUI_EditText    *edittext;
GLUI_CommandLine *cmd_line;
GLUI_Panel       *obj_panel;
GLUI_Button      *open_console_btn;

// Glui User IDs for callbacks
#define OPEN_CONSOLE_ID      100
#define CMD_HIST_RESET_ID    101
#define CMD_CLOSE_ID         102
#define LIGHT0_ENABLED_ID    200
#define LIGHT1_ENABLED_ID    201
#define LIGHT0_INTENSITY_ID  250
#define LIGHT1_INTENSITY_ID  251
#define OBJECT_TYPE_ID       4
#define RENDERING_MODE_ID    5
#define PROJECTION_TYPE_ID   6
#define	RESET_CAMERA_ID		 99

GLfloat light0_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
GLfloat light0_diffuse[] = { .6f, .6f, 1.0f, 1.0f };
GLfloat light0_position[] = { .5f, .5f, 1.0f, 0.0f };
GLfloat light1_ambient[] = { 0.1f, 0.1f, 0.3f, 1.0f };
GLfloat light1_diffuse[] = { .9f, .6f, 0.0f, 1.0f };
GLfloat light1_position[] = { -1.0f, -1.0f, 1.0f, 0.0f };

void swapEdgeList(int a, int b) // Used during sorting the edge pair list before identifying edge pairs
{
	int temp[3];
	temp[0] = edgepair[3 * a + 0];
	temp[1] = edgepair[3 * a + 1];
	temp[2] = edgepair[3 * a + 2];
	edgepair[3 * a + 0] = edgepair[3 * b + 0];
	edgepair[3 * a + 1] = edgepair[3 * b + 1];
	edgepair[3 * a + 2] = edgepair[3 * b + 2];
	edgepair[3 * b + 0] = temp[0];
	edgepair[3 * b + 1] = temp[1];
	edgepair[3 * b + 2] = temp[2];
}

int partitionEdgeList(int low, int high)// Used during quick sort of the edge pair list before identifying edge pairs
{
	int pivot = edgepair[3 * low + 1];
	int i;
	int leftwall = low;
	for (i = low + 1; i <= high; i++)
	{
		if (edgepair[3 * i + 1] < pivot)
		{
			leftwall = leftwall + 1;
			swapEdgeList(i, leftwall);
		}
	}
	swapEdgeList(low, leftwall);
	return leftwall;
}

void sortEdgeList(int low, int high) // Quick sort for the edge pair list before identifying edge pairs
{
	int pivot_location;
	if (low < high)
	{
		pivot_location = partitionEdgeList(low, high);
		sortEdgeList(low, pivot_location);
		sortEdgeList(pivot_location + 1, high);
	}
}

int searchEdgeList(int start, int end) // Binary search of edge pair list while identifying edge pairs
{
	int low = 0, high = he_e_count - 1, mid, temp;
	//cout << "Searching for (" << start << "," << end << ")" << endl;
	while (low <= high)
	{
		mid = low + (high - low) / 2;
		//cout << "Mid=" << mid << endl;
		
		if (start == edgepair[3 * mid + 1])
		{
			temp = mid;
			do
			{
				//cout << "Temp=" << temp << endl;
				
				if (end == edgepair[3 * temp + 2])
					return edgepair[3 * temp + 0];
				else
					temp++;
			} while (start == edgepair[3 * temp + 1]);
			temp = mid;
			do
			{
				//cout << "Temp=" << temp << endl;
				
				if (end == edgepair[3 * temp + 2])
					return edgepair[3 * temp + 0];
				else
					temp--;
			} while (start == edgepair[3 * temp + 1]);

			return -1;
		}
		else if (edgepair[3 * mid + 1] < start)
			low = mid + 1;
		else
			high = mid - 1;
	}
	return -1;
}

void halfEdgeInit() // Initializing half-edge data structure with 3D model
{
	int i, j, e;
	he_v_count = v_count;
	he_f_count = f_count;
	he_e_count = f_count * 3;
	m_HE_face = new HE_face[he_f_count];
	m_HE_vert = new HE_vert[he_v_count];
	m_HE_edge = new HE_edge[he_e_count];

	edgepair = new int[he_e_count * 3];

	for (i = 0; i < he_v_count; i++)
	{
		m_HE_vert[i].x = m_vertices[i].x;
		m_HE_vert[i].y = m_vertices[i].y;
		m_HE_vert[i].z = m_vertices[i].z;
		m_HE_vert[i].edge = NULL;
	}

	for (i = 0; i < he_f_count; i++)
	{
		m_HE_face[i].edge = NULL;
	}

	for (i = 0; i < he_e_count * 3; i++)
	{
		edgepair[i] = -1;
	}
	for (i = 0; i < f_count; i++)
	{
		m_HE_face[i].edge = &m_HE_edge[i * 3];

		for (j = 0; j <= 2; j++)
		{
			int vstart, vend;
			vstart = m_faces[i].v_indices[j] - 1;
			vend = m_faces[i].v_indices[(j + 1) % 3] - 1;
			if (!m_HE_vert[vstart].edge)
			{
				m_HE_vert[vstart].edge = &m_HE_edge[3 * i + j];
			}

			m_HE_edge[3 * i + j].vert = &m_HE_vert[vstart];
			m_HE_edge[3 * i + j].face = &m_HE_face[i];
			m_HE_edge[3 * i + j].prev = &m_HE_edge[(3 * i) + (j == 0 ? 2 : (j - 1))];
			m_HE_edge[3 * i + j].next = &m_HE_edge[(3 * i) + (j == 2 ? 0 : (j + 1))];
			m_HE_edge[3 * i + j].pair = NULL;
			e = 3 * i + j;
			edgepair[3 * e + 0] = 3 * i + j;
			edgepair[3 * e + 1] = vstart;
			edgepair[3 * e + 2] = vend;

		}
	}

	sortEdgeList(0, he_e_count - 1);
	
	int start, end;
	for (i = 0; i < he_e_count; i++)
	{
		e = edgepair[3 * i + 0];
		if (e >= 0)
		{
			start = edgepair[3 * i + 1];
			end = edgepair[3 * i + 2];

			if (!m_HE_edge[e].pair)
			{
				j = searchEdgeList(end, start);
				if (j >= 0)
				{
					m_HE_edge[e].pair = &m_HE_edge[j];
					m_HE_edge[j].pair = &m_HE_edge[e];
				}
			}
		}
	}
	delete edgepair;
}

void faceNormal() // Computing face normals every face in the 3D model
{
	int i;
	HE_vert a, b, c;

	float vect_ab[3], vect_ac[3];
	HE_edge *curr;
	float fnorm;
	for (i = 0; i < he_f_count; i++)
	{
		curr = m_HE_face[i].edge;
		a.x = curr->vert->x;
		a.y = curr->vert->y;
		a.z = curr->vert->z;

		b.x = curr->next->vert->x;
		b.y = curr->next->vert->y;
		b.z = curr->next->vert->z;

		c.x = curr->prev->vert->x;
		c.y = curr->prev->vert->y;
		c.z = curr->prev->vert->z;

		vect_ab[0] = b.x - a.x;
		vect_ab[1] = b.y - a.y;
		vect_ab[2] = b.z - a.z;

		vect_ac[0] = c.x - a.x;
		vect_ac[1] = c.y - a.y;
		vect_ac[2] = c.z - a.z;

		curr->face->fnx = (vect_ab[1] * vect_ac[2]) - (vect_ac[1] * vect_ab[2]);
		curr->face->fny = (vect_ab[2] * vect_ac[0]) - (vect_ac[2] * vect_ab[0]);
		curr->face->fnz = (vect_ab[0] * vect_ac[1]) - (vect_ac[0] * vect_ab[1]);

		fnorm = sqrt(curr->face->fnx*curr->face->fnx + curr->face->fny*curr->face->fny + curr->face->fnz*curr->face->fnz);

		curr->face->fnx /= fnorm;
		curr->face->fny /= fnorm;
		curr->face->fnz /= fnorm;
		curr->face->area = 0.5*fnorm;

	}
}

void vertexNormal() // Computing vertex normal for every vertex of the 3D model
{
	int i;
	HE_edge *curr;
	HE_vert *vert;

	float vnorm;
	for (i = 0; i < he_v_count; i++)
	{
		vert = &m_HE_vert[i];
		curr = vert->edge;
		vert->vnx = (curr->face->fnx*curr->face->area);
		vert->vny = (curr->face->fny*curr->face->area);
		vert->vnz = (curr->face->fnz*curr->face->area);
	
		while (1)
		{
			if (curr->pair)
				if (curr->pair->next == vert->edge)
					break;
				else
					curr = curr->pair->next;
			else
			{
				if (!curr->prev->pair)
					break;

				while (curr->prev->pair) { curr = curr->prev->pair; }
				if ((curr->prev == vert->edge) || (curr == vert->edge))
					break;
			}


			vert->vnx += (curr->face->fnx*curr->face->area);
			vert->vny += (curr->face->fny*curr->face->area);
			vert->vnz += (curr->face->fnz*curr->face->area);
		
		}

		vnorm = sqrt(vert->vnx*vert->vnx + vert->vny*vert->vny + vert->vnz*vert->vnz);

		vert->vnx /= vnorm;
		vert->vny /= vnorm;
		vert->vnz /= vnorm;
	}
}

void boundingBox() // Computing Axis-Aligned Bounding Box for 3D model
{
	int i;
	
	xmax = m_vertices[0].x;
	xmin = m_vertices[0].x;
	for (i = 0; i < v_count; i++)
	{
		if (m_vertices[i].x>xmax)
			xmax = m_vertices[i].x;
		if (m_vertices[i].x < xmin)
			xmin = m_vertices[i].x;
	}
	ymax = m_vertices[0].y;
	ymin = m_vertices[0].y;
	for (i = 0; i < v_count; i++)
	{
		if (m_vertices[i].y>ymax)
			ymax = m_vertices[i].y;
		if (m_vertices[i].y < ymin)
			ymin = m_vertices[i].y;
	}
	zmax = m_vertices[0].z;
	zmin = m_vertices[0].z;
	for (i = 0; i < v_count; i++)
	{
		if (m_vertices[i].z>zmax)
			zmax = m_vertices[i].z;
		if (m_vertices[i].z < zmin)
			zmin = m_vertices[i].z;
	}

	center.x = xmin + (xmax - xmin) / 2;
	center.y = ymin + (ymax - ymin) / 2;
	center.z = zmin + (zmax - zmin) / 2;

	if (abs(xmax - xmin)>abs(ymax - ymin))
	{
		if (abs(xmax - xmin) >abs(zmax - zmin))
			maxwidth = abs(xmax - xmin);
		else
			maxwidth = abs(zmax - zmin);

	}
	else
	{
		if (abs(ymax - ymin) > abs(zmax - zmin))
			maxwidth = abs(ymax - ymin);
		else
			maxwidth = abs(zmax - zmin);
	}
}

void moveToOrigin() // Moving centre of 3D model to the origin
{
	int i;
	for (i = 0; i < v_count; i++)
	{
		m_vertices[i].x = m_vertices[i].x - center.x;
		m_vertices[i].y = m_vertices[i].y - center.y;
		m_vertices[i].z = m_vertices[i].z - center.z;
	}

	xmax -= center.x;
	xmin -= center.x;
	ymax -= center.y;
	ymin -= center.y;
	zmax -= center.z;
	zmin -= center.z;
}

void scaleToUnitCube() // Scaling the 3D object to fit inside a unit cube
{
	int i;
	float scale;
	scale = 1 / maxwidth;
	for (i = 0; i < v_count; i++)
	{
		m_vertices[i].x = m_vertices[i].x *scale;
		m_vertices[i].y = m_vertices[i].y *scale;
		m_vertices[i].z = m_vertices[i].z *scale;
	}

	xmax *= scale;
	xmin *= scale;
	ymax *= scale;
	ymin *= scale;
	zmax *= scale;
	zmin *= scale;
}

int readFile(char* filename) //Parsing the M-file
{
	int line = 0;
	f_count = 0;
	v_count = 0;

	ifstream fin;
	fin.open(filename); // open a file
	if (!fin.good())
	{
		cout << "File not found!" << endl;
		return 0;
	}

	while (!fin.eof())
	{
		char buf[MAX_CHARS_PER_LINE];
		fin.getline(buf, MAX_CHARS_PER_LINE);
		line++;

	}
	fin.close();

	if (m_vertices) delete m_vertices;
	if (m_faces) delete m_faces;

	m_vertices = new Vertex[line];
	m_faces = new TriangleFace[line];
	fin.open(filename);

	while (!fin.eof())
	{
		char buf[MAX_CHARS_PER_LINE];
		fin.getline(buf, MAX_CHARS_PER_LINE);

		const char* token[MAX_TOKENS_PER_LINE] = {};

		int n = 0, i;
		token[0] = strtok(buf, DELIMITER);
		if (token[0])
		{
			for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
			{
				token[n] = strtok(0, DELIMITER);
				if (!token[n]) break;
			}
			if (strcmp(token[0], "#") != 0)
			{
				if (strcmp(token[0], "Vertex") == 0)
				{
					i = atoi(token[1]) - 1;
					m_vertices[i].x = stof(token[2]);
					m_vertices[i].y = stof(token[3]);
					m_vertices[i].z = stof(token[4]);
					v_count++;
				}
				if (strcmp(token[0], "Face") == 0)
				{
					i = atoi(token[1]) - 1;
					m_faces[i].v_indices[0] = atoi(token[2]);
					m_faces[i].v_indices[1] = atoi(token[3]);
					m_faces[i].v_indices[2] = atoi(token[4]);
					f_count++;
				}
			}
		}
	}

	boundingBox();
	moveToOrigin();
	scaleToUnitCube();

	halfEdgeInit();
	faceNormal();
	vertexNormal();
	return 1;
}

void control_cb(int control) // Callback function for Glui controls
{
	if (control == RENDERING_MODE_ID)
	{
		switch (rmode_type)
		{
		case 0:
			wireframe = 1;
			pointcloud = 0;
			flat = 0;
			smooth = 0;
			break;
		case 1:
			wireframe = 0;
			pointcloud = 1;
			flat = 0;
			smooth = 0;
			break;
		case 2:
			wireframe = 0;
			pointcloud = 0;
			flat = 1;
			smooth = 0;
			break;
		case 3:
			wireframe = 0;
			pointcloud = 0;
			flat = 0;
			smooth = 1;
			break;
		}


	}
	else if (control == LIGHT0_ENABLED_ID) {
		if (light0_enabled) {
			glEnable(GL_LIGHT0);
			light0_spinner->enable();
		}
		else {
			glDisable(GL_LIGHT0);
			light0_spinner->disable();
		}
	}
	else if (control == LIGHT1_ENABLED_ID) {
		if (light1_enabled) {
			glEnable(GL_LIGHT1);
			light1_spinner->enable();
		}
		else {
			glDisable(GL_LIGHT1);
			light1_spinner->disable();
		}
	}
	else if (control == LIGHT0_INTENSITY_ID) {
		float v[] = { light0_diffuse[0],  light0_diffuse[1],
			light0_diffuse[2],  light0_diffuse[3] };

		v[0] *= light0_intensity;
		v[1] *= light0_intensity;
		v[2] *= light0_intensity;

		glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
	}
	else if (control == LIGHT1_INTENSITY_ID) {
		float v[] = { light1_diffuse[0],  light1_diffuse[1],
			light1_diffuse[2],  light1_diffuse[3] };

		v[0] *= light1_intensity;
		v[1] *= light1_intensity;
		v[2] *= light1_intensity;

		glLightfv(GL_LIGHT1, GL_DIFFUSE, v);

	}
	else if (control == RESET_CAMERA_ID)
	{
		rotationX = 0.0; rotationY = 0.0; rotationZ = 0.0; rotationXZ = 0.0;
		moveX = 0.0; moveY = 0.0; moveZ = 0.0;
		zoom = 1.0;
	}
	if (rotobjX >= 360)
		rotobjX -= 360;
	if (rotobjY >= 360)
		rotobjY -= 360;
	if (rotobjZ >= 360)
		rotobjZ -= 360;
	if (rotobjX <= -360)
		rotobjX += 360;
	if (rotobjY <= -360)
		rotobjY += 360;
	if (rotobjZ <= -360)
		rotobjZ += 360;

	if (moveobjX >= 2)
		moveobjX = 2;
	if (moveobjY >= 2)
		moveobjY = 2;
	if (moveobjZ >= 2)
		moveobjZ = 2;
	if (moveobjX <= -2)
		moveobjX = -2;
	if (moveobjY <= -2)
		moveobjY = -2;
	if (moveobjZ <= -2)
		moveobjZ = -2;

}

void myGlutIdle(void)
{
	if (glutGetWindow() != main_window)
		glutSetWindow(main_window);

	glutPostRedisplay();

	glui->sync_live();

}

void myGlutMouse(int button, int button_state, int x, int y)
{
	old_x = x;
	old_y = y;
	if (button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN)
	{
		leftbutton = 1;
		rightbutton = 0;
		middlebutton = 0;
	}
	else if (button == GLUT_RIGHT_BUTTON && button_state == GLUT_DOWN)
	{
		leftbutton = 0;
		rightbutton = 1;
		middlebutton = 0;
	}
	else if (button == GLUT_MIDDLE_BUTTON && button_state == GLUT_DOWN)
	{
		leftbutton = 0;
		rightbutton = 0;
		middlebutton = 1;
	}
}

void myGlutMotion(int x, int y)
{
	int dx = old_x - x;
	int dy = old_y - y;

	//Camera Control
	if (rightbutton == 1) //Zoom in - Zoom out
	{		
			if (dx*dy > 0)
				zoom += (0.01*(dx + dy));
			else
				zoom += (0.01*((-1 * dx) + dy));
	
			if (proj_type == 0)
			{
				if (zoom > 2.5)
					zoom = 2.5;
				if (zoom < 0.2)
					zoom = 0.2;
			}
			else
			{
				if (zoom > 3.75)
					zoom = 3.75;
				if (zoom < 0.35)
					zoom = 0.35;
			}
	}
	else if (leftbutton == 1) // Rotate
	{
		if(abs(dy)>=10*abs(dx))
			rotationXZ -= 0.1*dy;
		else if (abs(dx)>10*abs(dy))
			rotationY += 0.1*dx;
		else 
		{
			if(dx*dy>0)
				rotationZ += 0.05*(dx + dy);
			else
				rotationX -= 0.05*((-1*dx) + dy);
		}

		if (rotationX >= 360)
			rotationX -= 360;
		if (rotationY >= 360)
			rotationY -= 360;
		if (rotationZ >= 360)
			rotationZ -= 360;
		if (rotationX <= -360)
			rotationX += 360;
		if (rotationY <= -360)
			rotationY += 360;
		if (rotationZ <= -360)
			rotationZ += 360;

	}
	else if (middlebutton == 1) // Translate
	{
		
		if (abs(dx) >= abs(2 * dy))
		{
			moveX -= 0.1*dx;
			moveZ += 0.1*dx;

		}
		else if (abs(dy) >= abs(2 * dx))
		{
			moveY += 0.1*dy;
		}
		else if (abs(abs(dy) - abs(dx)) <= 5)
		{
			if (dx*dy > 0)
				moveX -= 0.1*(dx + dy);
			else
				moveZ -= 0.1*((-1 * dx) + dy);
		}


		if (moveX >= 5)
			moveX = 5;
		if (moveY >= 5)
			moveY = 5;
		if (moveZ >= 5)
			moveZ = 5;
		if (moveX <= -5)
			moveX = -5;
		if (moveY <= -5)
			moveY = -5;
		if (moveZ <= -5)
			moveZ = -5;			
		
		
	}
	old_x = x;
	old_y = y;

	glutPostRedisplay();
}

void myGlutReshape(int x, int y)
{
	xy_aspect = (float)x / (float)y;
	GLUI_Master.auto_set_viewport();

	glutPostRedisplay();
}

void drawMesh() // Drawing the 3D model
{
	int i, j;
	glColor3d(1., 0., 0.);

	HE_edge *curr;

	glBegin(GL_TRIANGLES);
	for (int i = 0; i < he_f_count; i++)
	{
		curr = m_HE_face[i].edge;
		glColor3f(curr->prev->vert->x, curr->prev->vert->y, curr->prev->vert->z);
		glNormal3f(curr->prev->vert->vnx, curr->prev->vert->vny, curr->prev->vert->vnz);
		glVertex3f(curr->prev->vert->x, curr->prev->vert->y, curr->prev->vert->z);
		
		glColor3f(curr->vert->x, curr->vert->y, curr->vert->z);
		glNormal3f(curr->vert->vnx, curr->vert->vny, curr->vert->vnz);
		glVertex3f(curr->vert->x, curr->vert->y, curr->vert->z);

		glColor3f(curr->next->vert->x, curr->next->vert->y, curr->next->vert->z);
		glNormal3f(curr->next->vert->vnx, curr->next->vert->vny, curr->next->vert->vnz);
		glVertex3f(curr->next->vert->x, curr->next->vert->y, curr->next->vert->z);
	}
	glEnd();
}

void drawPlane() // Drawing XZ plane
{
	int i;
	float x, offset;
	x = -2;
	offset = 0.2;
	
	while (x < 2.2)
	{
		glColor3f(1.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex3f(x, 0.0, -2.0);
		glVertex3f(x, 0.0, 2.0);
		glVertex3f(-2.0, 0.0, x);
		glVertex3f(2.0, 0.0, x);
		x += offset;
		glEnd();
	}

}

void drawBoundingBox() // Drawing Axis-Aligned Bounding Box
{
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glColor3d(1., 0., 1.);

	glBegin(GL_LINES);

	glVertex3f(xmin, ymin, zmin);
	glVertex3f(xmax, ymin, zmin);
	glVertex3f(xmin, ymin, zmax);
	glVertex3f(xmax, ymin, zmax);
	glVertex3f(xmin, ymin, zmin);
	glVertex3f(xmin, ymin, zmax);
	glVertex3f(xmax, ymin, zmin);
	glVertex3f(xmax, ymin, zmax);

	glVertex3f(xmin, ymax, zmin);
	glVertex3f(xmax, ymax, zmin);
	glVertex3f(xmin, ymax, zmax);
	glVertex3f(xmax, ymax, zmax);
	glVertex3f(xmin, ymax, zmin);
	glVertex3f(xmin, ymax, zmax);
	glVertex3f(xmax, ymax, zmin);
	glVertex3f(xmax, ymax, zmax);

	glVertex3f(xmin, ymin, zmin);
	glVertex3f(xmin, ymax, zmin);
	glVertex3f(xmax, ymin, zmin);
	glVertex3f(xmax, ymax, zmin);
	glVertex3f(xmin, ymin, zmax);
	glVertex3f(xmin, ymax, zmax);
	glVertex3f(xmax, ymin, zmax);
	glVertex3f(xmax, ymax, zmax);

	glEnd();

}

void drawAxes() // Drawing coordinate axes
{
	GLUquadricObj *quadratic;
	quadratic = gluNewQuadric();
	
	//z-axis
	glPushMatrix();
	glColor3f(0, 0, 1);
	gluCylinder(quadratic, 0.01f, 0.01f, 1.0f, 32, 32);
	glTranslatef(0, 0, 1);
	glutSolidCone(0.02f, 0.06f, 8, 8);
	glPopMatrix();

	//x-axis
	glPushMatrix();
	glColor3f(1, 0, 0);
	glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
	gluCylinder(quadratic, 0.01f, 0.01f, 1.0f, 32, 32);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(1, 0, 0);
	glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
	glutSolidCone(0.02f, 0.06f, 8, 8);
	glPopMatrix();

	//y-axis
	glPushMatrix();
	glColor3f(0, 1, 0);
	glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
	gluCylinder(quadratic, 0.01f, 0.01f, 1.0f, 32, 32);
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, 1, 0);
	glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
	glutSolidCone(0.02f, 0.06f, 8, 8);
	glPopMatrix();

	glFlush();

}

void myGlutDisplay(void)
{
	glClearColor(.9f, .9f, .9f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//Determining projection type
	if (proj_type == 0)
	{
		gluPerspective(fov, 1, 0.1, 200);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(0.3f,1.2f, 0.5f, 0, 0, 0, 0, 1, 0);
	}
	else
	{
		glOrtho(Oleft, Oright, Obottom, Otop, 0.1, 50);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(10, 7, 10, 0, 0, 0, 0, 1, 0);
	}
	
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	//Camera Control

	glScalef(zoom, zoom, zoom);
	glTranslatef(moveX, moveY, moveZ);
	glRotatef(rotationXZ, 1, 0, -1);
	glRotatef(rotationX, 1, 0, 0);
	glRotatef(rotationY, 0, 1, 0);
	glRotatef(rotationZ, 0, 0, 1);

	if (showplane)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		drawPlane();
	}
	if (showaxes)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glShadeModel(GL_SMOOTH);
		drawAxes();
	}

	//Object transformations
	glScalef(scaleXYZ, scaleXYZ, scaleXYZ);
	glTranslatef(moveobjX, moveobjY, moveobjZ);
	glRotatef(rotobjX, 1, 0, 0);
	glRotatef(rotobjY, 0, 1, 0);
	glRotatef(rotobjZ, 0, 0, 1);

	//Determining and loading 3D object to be displayed
	if (obj_type == 0)
	{
		if (!bunny)
		{
			if (readFile("bunny.m") == 0) return;
			bunny = 1; cap = 0; gargoyle = 0; knot = 0; eight = 0;
		}
	}
	else if (obj_type == 1)
	{
		if (!cap)
		{
			if (readFile("cap.m") == 0) return;
			bunny = 0; cap = 1; gargoyle = 0; knot = 0; eight = 0;
		}
	}
	else if (obj_type == 2)
	{
		if (!eight)
		{
			if (readFile("eight.m") == 0) return;
			bunny = 0; cap = 0; gargoyle = 0; knot = 0; eight = 1;
		}
	}
	else if (obj_type == 3)
	{
		if (!gargoyle)
		{
			if (readFile("gargoyle.m") == 0) return;
			bunny = 0; cap = 0; gargoyle = 1; knot = 0; eight = 0;
		}
	}
	else if (obj_type == 4)
	{
		if (!knot)
		{
			if (readFile("knot.m") == 0) return;
			bunny = 0; cap = 0; gargoyle = 0; knot = 1; eight = 0;
		}
	}

	// Determining display mode
	if (wireframe)
	{
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else if (pointcloud)
		glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
	else if (flat)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glShadeModel(GL_FLAT);
	}
	else
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glShadeModel(GL_SMOOTH);
	}

	drawMesh();
	
	if (showbbox)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		drawBoundingBox();
	}

	glutSwapBuffers();
}

int main(int argc, char* argv[])
{

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(10, 10);
	glutInitWindowSize(1000, 750);

	main_window = glutCreateWindow("3D Mesh Viewer - Rungta Aditi - G1502070B");
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutReshape);
	glutMotionFunc(myGlutMotion);
	glutMouseFunc(myGlutMouse);

	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

	
	glEnable(GL_DEPTH_TEST);
	
	printf("GLUI version: %3.2f\n", GLUI_Master.get_version());

	glui = GLUI_Master.create_glui_subwindow(main_window, GLUI_SUBWINDOW_LEFT);

	new GLUI_StaticText(glui, "Toolbar");

	obj_panel = new GLUI_Panel(glui, "Object");

	radio = new GLUI_RadioGroup(obj_panel, &obj_type, OBJECT_TYPE_ID, control_cb);
	new GLUI_RadioButton(radio, "Bunny");
	new GLUI_RadioButton(radio, "Cap");
	new GLUI_RadioButton(radio, "Eight");
	new GLUI_RadioButton(radio, "Gargoyle");
	new GLUI_RadioButton(radio, "Knot");
	
	GLUI_Panel *rmode_panel = new GLUI_Panel(glui, "Rendering Mode");
	radio1 = new GLUI_RadioGroup(rmode_panel, &rmode_type, RENDERING_MODE_ID, control_cb);
	new GLUI_RadioButton(radio1, "Wireframe");
	new GLUI_RadioButton(radio1, "Point Cloud");
	new GLUI_RadioButton(radio1, "Flat Shading");
	new GLUI_RadioButton(radio1, "Smooth Shading");

	GLUI_Panel *proj_panel = new GLUI_Panel(glui, "Projection");
	radio2 = new GLUI_RadioGroup(proj_panel, &proj_type, PROJECTION_TYPE_ID, control_cb);
	new GLUI_RadioButton(radio2, "Perspective");
	new GLUI_RadioButton(radio2, "Orthogonal");

	GLUI_Panel *trans_panel = new GLUI_Panel(glui, "Object Transformations");
	scale_spinner =
		new GLUI_Spinner(trans_panel, "Scale    ", &scaleXYZ);
	scale_spinner->set_float_limits(.2f, 4.0);
	scale_spinner->set_alignment(GLUI_ALIGN_LEFT);

	GLUI_Panel *move_panel = new GLUI_Panel(trans_panel, "Translate");
	moveX_spinner =
		new GLUI_Spinner(move_panel, "X", &moveobjX, 10, control_cb);
	moveX_spinner->set_float_limits(-10.0, 10.0);
	moveX_spinner->set_alignment(GLUI_ALIGN_LEFT);
	moveY_spinner =
		new GLUI_Spinner(move_panel, "Y", &moveobjY, 11, control_cb);
	moveY_spinner->set_float_limits(-10.0, 10.0);
	moveY_spinner->set_alignment(GLUI_ALIGN_LEFT);
	moveZ_spinner =
		new GLUI_Spinner(move_panel, "Z", &moveobjZ, 12, control_cb);
	moveZ_spinner->set_float_limits(-10.0, 10.0);
	moveZ_spinner->set_alignment(GLUI_ALIGN_LEFT);

	GLUI_Panel *rot_panel = new GLUI_Panel(trans_panel, "Rotate");
	rotX_spinner =
		new GLUI_Spinner(rot_panel, "X", &rotobjX, 20, control_cb);
	rotX_spinner->set_alignment(GLUI_ALIGN_LEFT);
	rotY_spinner =
		new GLUI_Spinner(rot_panel, "Y", &rotobjY, 21, control_cb);
	rotY_spinner->set_alignment(GLUI_ALIGN_LEFT);
	rotZ_spinner =
		new GLUI_Spinner(rot_panel, "Z", &rotobjZ, 22, control_cb);
	rotZ_spinner->set_alignment(GLUI_ALIGN_LEFT);

	checkbox1 =
		new GLUI_Checkbox(glui, "Show Coordinate Axes", &showaxes, 1, control_cb);
	checkbox2 =
		new GLUI_Checkbox(glui, "Show XZ-Plane", &showplane, 2, control_cb);
	checkbox3 =
		new GLUI_Checkbox(glui, "Show Bounding Box", &showbbox, 3, control_cb);

	GLUI_Panel *light0 = new GLUI_Panel(glui, "Light 1");
	GLUI_Panel *light1 = new GLUI_Panel(glui, "Light 2");

	new GLUI_Checkbox(light0, "Enabled", &light0_enabled,
		LIGHT0_ENABLED_ID, control_cb);
	light0_spinner =
		new GLUI_Spinner(light0, "Intensity:",
			&light0_intensity, LIGHT0_INTENSITY_ID,
			control_cb);
	light0_spinner->set_float_limits(0.0, 1.0);

	new GLUI_Checkbox(light1, "Enabled", &light1_enabled,
		LIGHT1_ENABLED_ID, control_cb);
	light1_spinner =
		new GLUI_Spinner(light1, "Intensity:",
			&light1_intensity, LIGHT1_INTENSITY_ID,
			control_cb);
	light1_spinner->set_float_limits(0.0, 1.0);
	
	new GLUI_Button(glui, "Reset Camera", RESET_CAMERA_ID, control_cb);
	new GLUI_Button(glui, "Quit", 0, (GLUI_Update_CB)exit);

	glui->set_main_gfx_window(main_window);

	GLUI_Master.set_glutIdleFunc(myGlutIdle);
	
	glutMainLoop();

	return EXIT_SUCCESS;
}
