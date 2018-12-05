#include <Windows.h>
#include <math.h> 
#include <iostream>
#include "resource.h"
#include <fstream>
#include <stack>
#include <vector>
#include <cstring>
using namespace std;


struct Shape 
{
	char name[20];
	int xs;
	int ys;
	int xe;
	int ye;
	int x1;
	int y1;
	int x2;
	int y2;
	int prev;
};
struct Line
{
	int xs, xe, ys, ye;
};
vector<Line> lines2;
int width = 1250;
int height = 650;
int shapeSize = 0;
Shape shapes[1000];

// Hermite Curves 
struct Vector2 { 
	double x, y;  
	Vector2(double a = 0, double b = 0)
	{
		x = a;
		y = b; 
	} 
}; 
class Vector4 
{
	double v[4];
public:  
	Vector4(double a = 0, double b = 0, double c = 0, double d = 0) 
	{
		v[0] = a;
		v[1] = b;
		v[2] = c;
		v[3] = d; 
	}
	Vector4(double a[])
	{ 
		memcpy(v, a, 4 * sizeof(double));
	}  
	double& operator[](int i) 
	{ 
		return v[i];
	}
}; 
class Matrix4 
{ 
	Vector4 M[4]; 
public:  
	Matrix4(double A[]) 
	{ 
		memcpy(M, A, 16 * sizeof(double)); 
	}  
	Vector4& operator[](int i)
	{ 
		return M[i]; 
	} 
};
Vector4 operator*(Matrix4 M, Vector4& b) // right multiplication of M by b 
{  
	Vector4 res;  
	for (int i = 0; i < 4; i++) 
	{
		for (int j = 0; j < 4; j++)
		{
			res[i] += M[i][j] * b[j];
		}
	}
	return res; 
} 
double DotProduct(Vector4& a, Vector4& b) //multiplying a raw vector by a column vector 
{  
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]; 
} 
Vector4 GetHermiteCoeff(double x0, double s0, double x1, double s1) 
{ 
	static double H[16] = { 2,1,-2,1,-3,-2,3,-1,0,1,0,0,1,0,0,0 }; 
	static Matrix4 basis(H); 
	Vector4 v(x0, s0, x1, s1); 
	return basis * v;
}
void DrawHermiteCurve(HDC hdc, Vector2& P0, Vector2& T0, Vector2& P1, Vector2& T1, int numpoints) 
{ 
	Vector4 xcoeff = GetHermiteCoeff(P0.x, T0.x, P1.x, T1.x); 
	Vector4 ycoeff = GetHermiteCoeff(P0.y, T0.y, P1.y, T1.y); 
	if (numpoints<2)return;  double dt = 1.0 / (numpoints - 1); 
	for (double t = 0; t <= 1; t += dt)
	{ 
		Vector4 vt;  
		vt[3] = 1;   
		for (int i = 2; i >= 0; i--)
		{
			vt[i] = vt[i + 1] * t;
		} 
		int x = round(DotProduct(xcoeff, vt));  
		int y = round(DotProduct(ycoeff, vt));   
		if (t == 0)
		{
			MoveToEx(hdc, x, y, NULL);
		}
		else
		{
			LineTo(hdc, x, y);
		}
	} 
}
// Bezier Curve
void DrawBezierCurve(HDC hdc, Vector2& P0, Vector2& P1, Vector2& P2, Vector2& P3, int numpoints) 
{
	Vector2 T0(3 * (P1.x - P0.x), 3 * (P1.y - P0.y));  
	Vector2 T1(3 * (P3.x - P2.x), 3 * (P3.y - P2.y));  
	DrawHermiteCurve(hdc, P0, T0, P3, T1, numpoints); 
}
// Polygon Clipping
struct Vertex1
{
	double x, y;
	Vertex1(int x1 = 0, int y1 = 0)
	{
		x = x1;
		y = y1;
	}
};
vector<Line> lines;
vector<Vertex1> points;
typedef vector<Vertex1> VertexList;
typedef bool(*IsInFunc)(Vertex1& v, int edge);
typedef Vertex1(*IntersectFunc)(Vertex1& v1, Vertex1& v2, int edge);
VertexList ClipWithEdge(VertexList p, int edge, IsInFunc In, IntersectFunc Intersect) 
{ 
	VertexList OutList;   
	Vertex1 v1 = p[p.size() - 1];   
	bool v1_in = In(v1, edge);   
	for (int i = 0; i<(int)p.size(); i++) 
	{ 
		Vertex1 v2 = p[i];    
		bool v2_in = In(v2, edge);    
		if (!v1_in && v2_in) 
		{ 
			OutList.push_back(Intersect(v1, v2, edge));      
			OutList.push_back(v2); 
		}
		else if (v1_in && v2_in) 
		{
			OutList.push_back(v2);
		}  
		else if (v1_in) 
		{
			OutList.push_back(Intersect(v1, v2, edge));
		} 
		v1 = v2;   
		v1_in = v2_in;
	}
	return OutList; 
}
bool InLeft(Vertex1& v, int edge) 
{ 
	return v.x >= edge; 
} 
bool InRight(Vertex1& v, int edge)
{ 
	return v.x <= edge; 
}
bool InTop(Vertex1& v, int edge) 
{ 
	return v.y >= edge; 
} 
bool InBottom(Vertex1& v, int edge) 
{ 
	return v.y <= edge; 
}
Vertex1 VIntersect(Vertex1& v1, Vertex1& v2, int xedge) 
{ 
	Vertex1 res;  
	res.x = xedge;  
	res.y = v1.y + (xedge - v1.x)*(v2.y - v1.y) / (v2.x - v1.x);  
	return res; 
}
Vertex1 HIntersect(Vertex1& v1, Vertex1& v2, int yedge) 
{ 
	Vertex1 res; 
	res.y = yedge; 
	res.x = v1.x + (yedge - v1.y)*(v2.x - v1.x) / (v2.y - v1.y);  
	return res; 
}
void PolygonClip(HDC hdc, vector<Vertex1> &p, int n, int xleft, int ytop, int xright, int ybottom)
{ 
	VertexList vlist;  
	for (int i = 0; i < n; i++) 
	{
		vlist.push_back(Vertex1(p[i].x, p[i].y));
	}
	vlist = ClipWithEdge(vlist, xleft, InLeft, VIntersect); 
	vlist = ClipWithEdge(vlist, ytop, InTop, HIntersect); 
	vlist = ClipWithEdge(vlist, xright, InRight, VIntersect); 
	vlist = ClipWithEdge(vlist, ybottom, InBottom, HIntersect); 
	Vertex1 v1 = vlist[vlist.size() - 1];  
	for (int i = 0; i<(int)vlist.size(); i++) 
	{ 
		Vertex1 v2 = vlist[i];  
		MoveToEx(hdc, round(v1.x), round(v1.y), NULL);  
		LineTo(hdc, round(v2.x), round(v2.y));
		Line line;
		line.xs = round(v1.x);
		line.ys = round(v1.y);
		line.xe = round(v2.x);
		line.ye = round(v2.y);
		lines2.push_back(line);
		v1 = v2; 
	} 
}
void swap(int &x, int &y)
{
	int temp = x;
	x = y;
	y = temp;
}
void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
	SetPixel(hdc, xc + a, yc + b, color);
	SetPixel(hdc, xc - a, yc + b, color);
	SetPixel(hdc, xc - a, yc - b, color);
	SetPixel(hdc, xc + a, yc - b, color);
	SetPixel(hdc, xc + b, yc + a, color);
	SetPixel(hdc, xc - b, yc + a, color);
	SetPixel(hdc, xc - b, yc - a, color);
	SetPixel(hdc, xc + b, yc - a, color);
}
void CircleDirect(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = 0, y = R;
	int R2 = R * R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x<y)
	{
		x++;
		y = round(sqrt((double)(R2 - x * x)));
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}
void CirclePolar(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = R, y = 0;
	double theta = 0, dtheta = 1.0 / R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x>y)
	{
		theta += dtheta;
		x = round(R*cos(theta));
		y = round(R*sin(theta));
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}
void CircleIterativePolar(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	double x = R, y = 0;
	double dtheta = 1.0 / R;
	double cdtheta = cos(dtheta), sdtheta = sin(dtheta);
	Draw8Points(hdc, xc, yc, R, 0, color);
	while (x>y)
	{
		double x1 = x * cdtheta - y * sdtheta;
		y = x * sdtheta + y * cdtheta;
		x = x1;
		Draw8Points(hdc, xc, yc, round(x), round(y), color);
	}
}
void CircleBresenham(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = 0, y = R;
	int d = 1 - R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x<y)
	{
		if (d < 0)
		{
			d += 2 * x + 3;
		}
		else
		{
			d += 2 * (x - y) + 5;
			y--;
		}
		x++;
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}
void CircleFasterBresenham(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = 0, y = R;
	int d = 1 - R;
	int c1 = 3, c2 = 5 - 2 * R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x<y)
	{
		if (d<0)
		{
			d += c1;
			c2 += 2;
		}
		else
		{
			d += c2;
			c2 += 4;
			y--;
		}
		c1 += 2;
		x++;
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}
void MidPoint(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{
	if (xs>xe)
	{
		swap(xs, xe);
		swap(ys, ye);
	}
	/*if (ys>ye)
	{
	swap(xs, xe);
	swap(ys, ye);
	}*/

	int dx = xe - xs;
	int dy = ye - ys;
	int x = xs;
	int y = ys;
	int d, ch1, ch2;
	int increment;
	if (abs(dx) >= abs(dy))
	{
		if (dy < 0) {
			increment = -1;
			dy = -dy;
		}
		else {
			increment = 1;
		}
		d = dx - 2 * dy;
		ch1 = 2 * (dx - dy);
		ch2 = -2 * dy;
		SetPixel(hdc, x, y, color);
		while (x < xe)
		{
			if (d <= 0)
			{
				d += ch1;
				y += increment;
			}
			else
			{
				d += ch2;
			}
			x++;
			SetPixel(hdc, x, y, color);
		}
	}
	else
	{
		if (dx < 0) {
			increment = -1;
			dx = -dx;
		}
		else {
			increment = 1;
		}
		d = 2 * dx - dy;
		ch1 = 2 * (dx - dy);
		ch2 = 2 * dx;
		SetPixel(hdc, x, y, color);
		while (y < ye)
		{
			if (d >= 0)
			{
				d += ch1;
				x += increment;
			}
			else
			{
				d += ch2;
			}
			y++;
			SetPixel(hdc, x, y, color);
		}
	}
}
void SimpleDDA(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{
	int dx = xe - xs;
	int dy = ye - ys;
	SetPixel(hdc, xs, ys, color);
	if (abs(dx) >= abs(dy))
	{
		int x = xs, xinc = dx>0 ? 1 : -1;
		double y = ys, yinc = (double)dy / dx * xinc;
		while (x != xe)
		{
			x += xinc;
			y += yinc;
			SetPixel(hdc, x, round(y), color);
		}
	}
	else
	{
		int y = ys, yinc = dy>0 ? 1 : -1;
		double x = xs, xinc = (double)dx / dy * yinc;
		while (y != ye)
		{
			x += xinc;
			y += yinc;
			SetPixel(hdc, round(x), y, color);
		}
	}
}
void DrawLine(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color)
{
	int dx = xe - xs;
	int dy = ye - ys;
	if (abs(dy) <= abs(dx))
	{
		double slope = (double)dy / dx;
		if (xs>xe)
		{
			swap(xs, xe);
			swap(ys, ye);
		}
		for (int x = xs; x <= xe; x++)
		{
			int y = round(ys + (x - xs)*slope);
			SetPixel(hdc, x, y, color);
		}
	}
	else
	{
		double islope = (double)dx / dy;
		if (ys>ye)
		{
			swap(xs, xe);
			swap(ys, ye);
		}
		for (int y = ys; y <= ye; y++)
		{
			int x = round(xs + (y - ys)*islope);
			SetPixel(hdc, x, y, color);
		}

	}
}
int CalculateDistance(int xs, int ys, int xe, int ye)
{
	return sqrt(pow(abs(xe - xs), 2) + pow(abs(ye - ys), 2));
}
// Flood Filling
struct Vertex
{
	int x, y;
	Vertex(int X, int Y)
	{
		x = X;
		y = Y;
	}
};
void NRFloodFill(HDC hdc, int x, int y, COLORREF Cb, COLORREF Cf)
{
	stack<Vertex> S;
	S.push(Vertex(x, y));
	while (!S.empty()) {
		Vertex v = S.top();
		S.pop();
		COLORREF c = GetPixel(hdc, v.x, v.y);
		if (c == Cb || c == Cf)
		{
			continue;
		}
		SetPixel(hdc, v.x, v.y, Cf);
		S.push(Vertex(v.x + 1, v.y));
		S.push(Vertex(v.x - 1, v.y));
		S.push(Vertex(v.x, v.y + 1));
		S.push(Vertex(v.x, v.y - 1));
	}
}
// Line Clipping
union OutCode {
	unsigned All : 4;
	struct
	{
		unsigned left : 1, top : 1, right : 1, bottom : 1;
	};
};
OutCode GetOutCode(double x, double y, int xleft, int ytop, int xright, int ybottom)
{
	OutCode out;
	out.All = 0;
	if (x < xleft)
	{
		out.left = 1;
	}
	else if (x > xright)
	{
		out.right = 1;
	}
	if (y < ytop)
	{
		out.top = 1;
	}
	else if (y > ybottom)
	{
		out.bottom = 1;
	}
	return out;
}
void VIntersect(double xs, double ys, double xe, double ye, int x, double *xi, double *yi)
{
	*xi = x;
	*yi = ys + (x - xs)*(ye - ys) / (xe - xs);
}
void HIntersect(double xs, double ys, double xe, double ye, int y, double *xi, double *yi)
{
	*yi = y;
	*xi = xs + (y - ys)*(xe - xs) / (ye - ys);
}
void CohenSuth(HDC hdc, int xs, int ys, int xe, int ye, int xleft, int ytop, int xright, int ybottom)
{
	double x1 = xs, y1 = ys, x2 = xe, y2 = ye;
	OutCode out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
	OutCode out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
	while ((out1.All || out2.All) && !(out1.All & out2.All))
	{
		double xi, yi;
		if (out1.All)
		{
			if (out1.left)
			{
				VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			}
			else if (out1.top)
			{
				HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			}
			else if (out1.right)
			{
				VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			}
			else
			{
				HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			}
			x1 = xi;
			y1 = yi;
			SetPixel(hdc, x1, y1, RGB(255, 255, 255));
			out1 = GetOutCode(x1, y1, xleft, ytop, xright, ybottom);
		}
		else
		{
			if (out2.left)
			{
				VIntersect(x1, y1, x2, y2, xleft, &xi, &yi);
			}
			else if (out2.top)
			{
				HIntersect(x1, y1, x2, y2, ytop, &xi, &yi);
			}
			else if (out2.right)
			{
				VIntersect(x1, y1, x2, y2, xright, &xi, &yi);
			}
			else
			{
				HIntersect(x1, y1, x2, y2, ybottom, &xi, &yi);
			}
			x2 = xi;
			y2 = yi;
			SetPixel(hdc, x2, y2, RGB(255, 255, 255));
			out2 = GetOutCode(x2, y2, xleft, ytop, xright, ybottom);
		}
	}
	if (!out1.All && !out2.All)
	{
		MoveToEx(hdc, round(x1), round(y1), NULL);
		LineTo(hdc, round(x2), round(y2));
		Line line;
		line.xs = round(x1);
		line.ys = round(y1);
		line.xe = round(x2);
		line.ye = round(y2);
		lines2.push_back(line);
	}
}
void SaveShape(string name, int xs, int ys, int xe, int ye, int x1, int y1, int x2, int y2, int  prev)
{
	Shape shape;
	strcpy_s(shape.name, name.c_str());
	shape.xs = xs;
	shape.ys = ys;
	shape.xe = xe;
	shape.ye = ye;
	shape.x1 = x1;
	shape.y1 = y1;
	shape.x2 = x2;
	shape.y2 = y2;
	shape.prev = prev;
	if (strcmp(shape.name, "LineClipping") == 0 || strcmp(shape.name, "PolygonClipping") == 0)
	{
		//shapeSize -= prev;
		shapes[shapeSize-prev] = shape;
		//shapeSize++;
		
	}
	else
	{
		shapes[shapeSize++] = shape;
	}
}

void LoadShape(Shape shape, HDC hdc)
{
	if (strcmp(shape.name, "DirectLine") == 0)
	{
		DrawLine(hdc, shape.xs, shape.ys, shape.xe, shape.ye, RGB(0, 0, 0));
	}
	else if (strcmp(shape.name, "DDALine") == 0)
	{
		SimpleDDA(hdc, shape.xs, shape.ys, shape.xe, shape.ye, RGB(0, 0, 0));
	}
	else if (strcmp(shape.name, "MidPointLine") == 0)
	{
		MidPoint(hdc, shape.xs, shape.ys, shape.xe, shape.ye, RGB(0, 0, 0));
	}
	else if (strcmp(shape.name, "QuadraticCircle") == 0)
	{
		CircleDirect(hdc, shape.xs, shape.ys, CalculateDistance(shape.xs, shape.ys, shape.xe, shape.ye), RGB(0, 0, 0));
	}
	else if (strcmp(shape.name, "PolarCircle") == 0)
	{
		CircleIterativePolar(hdc, shape.xs, shape.ys, CalculateDistance(shape.xs, shape.ys, shape.xe, shape.ye), RGB(0, 0, 0));
	}
	else if (strcmp(shape.name, "MidPointCircle") == 0)
	{
		CircleFasterBresenham(hdc, shape.xs, shape.ys, CalculateDistance(shape.xs, shape.ys, shape.xe, shape.ye), RGB(0, 0, 0));
	}
	else if (strcmp(shape.name, "FlatFill") == 0)
	{
		NRFloodFill(hdc, shape.xs, shape.ys, RGB(0, 0, 0), RGB(0, 0, 0));
	}
	else if (strcmp(shape.name, "LineClipping") == 0)
	{
		//CohenSuth(hdc, shape.xs, shape.ys, shape.xe, shape.ye, shape.x1, shape.y2, shape.x2, shape.y1);
		for (int i = 0; i < lines2.size(); i++)
		{
			Line line = lines2.at(i);
			MoveToEx(hdc, line.xs, line.ys, NULL);
			LineTo(hdc, line.xe, line.ye);
		}
	}
	else if (strcmp(shape.name, "PolygonClipping") == 0)
	{
		//PolygonClip(hdc, points, points.size(), shape.xs, shape.ye, shape.xe, shape.ys);
		for (int i = 0; i < lines2.size(); i++)
		{
			Line line = lines2.at(i);
			MoveToEx(hdc, line.xs, line.ys, NULL);
			LineTo(hdc, line.xe, line.ye);
		}
	}
	else if (strcmp(shape.name, "BezeirCurve") == 0)
	{
		Vector2 v1;
		v1.x = shape.xs;
		v1.y = shape.ys;
		Vector2 v2;
		v2.x = shape.x1;
		v2.y = shape.y1;
		Vector2 v3;
		v3.x = shape.x2;
		v3.y = shape.y2;
		Vector2 v4;
		v2.x = shape.xe;
		v2.y = shape.ye;
		DrawBezierCurve(hdc, v1, v4, v2, v3, 100);
	}
	else if (strcmp(shape.name, "HermitCurve") == 0)
	{
		Vector2 v1;
		v1.x = shape.xs;
		v1.y = shape.ys;
		Vector2 v2;
		v2.x = shape.x1;
		v2.y = shape.y1;
		Vector2 v3;
		v3.x = shape.x2;
		v3.y = shape.y2;
		Vector2 v4;
		v2.x = shape.xe;
		v2.y = shape.ye;
		DrawHermiteCurve(hdc, v1, v2, v3, v4, 100);
	}
}
void SaveFile()
{
	fstream file;
	file.open("draw.dr", ios::out | ios::binary | ios::trunc);
	int size = lines2.size();
	file.write((char *)& size, sizeof(int));
	for (int i = 0; i < size; i++) {
		file.write((char *)&lines2.at(i), sizeof(Line));
	}
	file.write((char *)&shapeSize, sizeof(int));
	for (int i = 0; i < shapeSize; i++) {
		file.write((char *)&shapes[i], sizeof(Shape));
	}
	file.close();
}
void LoadFile(HDC hdc)
{
	fstream file;
	file.open("draw.dr", ios::in | ios::binary);
	int size;
	file.read((char *)&size, sizeof(int));
	for (int i = 0; i < size; i++) {
		Line line;
		file.read((char *)&line, sizeof(Line));
		lines2.push_back(line);
	}
	file.read((char *)&shapeSize, sizeof(int));
	for (int i = 0; i < shapeSize; i++) {
		Shape shape;
		file.read((char *)&shape, sizeof(Shape));
		shapes[i] = shape;
		LoadShape(shape, hdc);
	}
	file.close();
}
LRESULT WINAPI WndProc(HWND hWnd/*msg reciever*/, UINT mcode, WPARAM wp/*the key i pressed for example with the mouse click*/, LPARAM lp/* position of click for exaple(x,y) cooerdinate*/)
{

	PAINTSTRUCT ps;
	HDC hdc;
	static int wmId, wmEvent;
	static int xs = 0, ys = 0, xe = 0, ye = 0, x1 = 0, y1 = 0, x2 = 0, y2 = 0, x3 = 0, y3 = 0;
	static bool saved = true;
	static bool load = false;
	static bool clip = false;
	static bool flag = false;
	static bool flag2 = false;
	static bool draw = false;
	static bool DirectLine = false;
	static bool DDALine = false;
	static bool MidPointLine = false;
	static bool QuadraticCircle = false;
	static bool PolarCircle = false;
	static bool MidPointCircle = false;
	static bool BezeirCurve = false;
	static bool HermitCurve = false;
	static bool FlatFill = false;
	static bool LineClipping = false;
	static bool PolygonClipping = false;
	switch (mcode)
	{
	case WM_COMMAND:
		wmId = LOWORD(wp);
		wmEvent = HIWORD(wp);
		switch (wmId)
		{
		case ID_FILE_SAVE:
			SaveFile();
			saved = true;
			MessageBox(hWnd, L"File Saved Completed", L"Save", MB_OK);
			break;
		case ID_FILE_LOAD:
			load = true;
			RedrawWindow(hWnd, NULL, NULL, RDW_INVALIDATE);
			break;
		case ID_LINE_DIRECT:
			DirectLine = true;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_LINE_DDA:
			DirectLine = false;
			DDALine = true;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_LINE_MIDPOINT:
			DirectLine = false;
			DDALine = false;
			MidPointLine = true;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_CIRCLE_QUADRATIC:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = true;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_CIRCLE_POLAR:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = true;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_CIRCLE_MIDPOINT:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = true;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_CURVES_BEZEIR:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = true;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_CURVES_HERMIT:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = true;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_FILLING_FLATFILLING:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = true;
			LineClipping = false;
			PolygonClipping = false;
			break;
		case ID_CLIPPING_LINE:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = true;
			PolygonClipping = false;
			RedrawWindow(hWnd, NULL, NULL, RDW_INVALIDATE);
			break;
		case ID_CLIPPING_POLYGON:
			DirectLine = false;
			DDALine = false;
			MidPointLine = false;
			QuadraticCircle = false;
			PolarCircle = false;
			MidPointCircle = false;
			BezeirCurve = false;
			HermitCurve = false;
			FlatFill = false;
			LineClipping = false;
			PolygonClipping = true;
			break;
		case ID_FILE_EXIT:
			if (!saved) {
				if (MessageBox(hWnd, L"Do You Want To Save Changes", L"Save", MB_YESNO) == IDYES) {
					SaveFile();
				}
			}
			DestroyWindow(hWnd);
			break;
		default:
			return DefWindowProc(hWnd, mcode, wp, lp);
		}
		break;

	case WM_PAINT:
		hdc = BeginPaint(hWnd, &ps);
		if (load)
		{
			LoadFile(hdc);
			load = false;
		}
		if (FlatFill) {
			NRFloodFill(hdc, xs, ys, RGB(0, 0, 0), RGB(0, 0, 0));
			SaveShape("FlatFill", xs, ys, 0, 0, x1, y1, x2, y2, 0);
			saved = false;
		}
		if (draw) {
			if (DirectLine) 
			{
				Line line;
				line.xs = xs;
				line.ys = ys;
				line.xe = xe;
				line.ye = ye;
				Vertex1 p;
				Vertex1 p2;
				p.x = xs;
				p.y = ys;
				p2.x = xe;
				p2.y = ye;
				lines.push_back(line);
				points.push_back(p);
				points.push_back(p2);
				DrawLine(hdc, xs, ys, xe, ye, RGB(0, 0, 0));
				SaveShape("DirectLine", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
			else if (DDALine) 
			{
				Line line;
				line.xs = xs;
				line.ys = ys;
				line.xe = xe;
				line.ye = ye;
				Vertex1 p;
				Vertex1 p2;
				p.x = xs;
				p.y = ys;
				p2.x = xe;
				p2.y = ye;
				lines.push_back(line);
				points.push_back(p);
				points.push_back(p2);
				SimpleDDA(hdc, xs, ys, xe, ye, RGB(0, 0, 0));
				SaveShape("DDALine", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
			else if (MidPointLine)
			{
				Line line;
				line.xs = xs;
				line.ys = ys;
				line.xe = xe;
				line.ye = ye;
				Vertex1 p;
				Vertex1 p2;
				p.x = xs;
				p.y = ys;
				p2.x = xe;
				p2.y = ye;
				lines.push_back(line);
				points.push_back(p);
				points.push_back(p2);
				MidPoint(hdc, xs, ys, xe, ye, RGB(0, 0, 0));
				SaveShape("MidPointLine", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
			else if (QuadraticCircle) 
			{
				CircleDirect(hdc, xs, ys, CalculateDistance(xs, ys, xe, ye), RGB(0, 0, 0));
				SaveShape("QuadraticCircle", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
			else if (PolarCircle)
			{
				CircleIterativePolar(hdc, xs, ys, CalculateDistance(xs, ys, xe, ye), RGB(0, 0, 0));
				SaveShape("PolarCircle", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
			else if (MidPointCircle) 
			{
				CircleFasterBresenham(hdc, xs, ys, CalculateDistance(xs, ys, xe, ye), RGB(0, 0, 0));
				SaveShape("MidPointCircle", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
			else if (LineClipping) 
			{
				CohenSuth(hdc, lines.at(lines.size()-1).xs, lines.at(lines.size() - 1).ys, lines.at(lines.size() - 1).xe, lines.at(lines.size() - 1).ye, xs, ye, xe, ys);
				SaveShape("LineClipping", lines.at(0).xs, lines.at(0).ys, lines.at(0).xe, lines.at(0).ye, xs, ys, xe, ye, 1);
				saved = false;
			}
			else if (PolygonClipping) 
			{
				PolygonClip(hdc, points, points.size(), xs, ye, xe, ys);
				SaveShape("PolygonClipping", xs, ys, xe, ye, x1, y1, x2, y2, lines.size());
				saved = false;
			}
			else if (HermitCurve) 
			{
				Vector2 v1;
				v1.x = xs;
				v1.y = ys;
				Vector2 v2;
				v2.x = x1;
				v2.y = y1;
				Vector2 v3;
				v3.x = x2;
				v3.y = y2;
				Vector2 v4;
				v2.x = xe;
				v2.y = ye;
				DrawHermiteCurve(hdc, v1, v2, v3, v4, 100);
				SaveShape("HermitCurve", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
			else if (BezeirCurve) 
			{
				Vector2 v1;
				v1.x = xs;
				v1.y = ys;
				Vector2 v2;
				v2.x = x1;
				v2.y = y1;
				Vector2 v3;
				v3.x = x2;
				v3.y = y2;
				Vector2 v4;
				v2.x = xe;
				v2.y = ye;
				DrawBezierCurve(hdc, v1, v4, v2, v3, 100);
				SaveShape("BezeirCurve", xs, ys, xe, ye, x1, y1, x2, y2, 0);
				saved = false;
			}
		}
		EndPaint(hWnd, &ps);
		break;

	case WM_LBUTTONDOWN:
		if (HermitCurve || BezeirCurve)
		{
			if (!flag && !flag2)
			{
				xs = LOWORD(lp);
				ys = HIWORD(lp);
				flag = true;
			}
			else if (flag && !flag2)
			{
				xe = LOWORD(lp);
				ye = HIWORD(lp);
				flag2 = true;
			}
			else if (flag && flag2)
			{
				x1 = LOWORD(lp);
				y1 = HIWORD(lp);
				flag = false;
			}
			else if (!flag && flag2)
			{
				x2 = LOWORD(lp);
				y2 = HIWORD(lp);
				flag2 = false;
				draw = true;
				RedrawWindow(hWnd, NULL, NULL, RDW_INVALIDATE);
			}
			
		}
		else 
		{
			if (!flag)
			{
				xs = LOWORD(lp);		//x-coordinate LOWWORD function to get x
				ys = HIWORD(lp);		//y-coordinate HIWORD function to get y
				if (FlatFill)
				{
					flag = false;
					RedrawWindow(hWnd, NULL, NULL, RDW_INVALIDATE);
				}
				else
				{
					flag = true;
				}
			}
			else
			{
				xe = LOWORD(lp);		//x-coordinate LOWWORD function to get x
				ye = HIWORD(lp);		//y-coordinate HIWORD function to get y
				draw = true;
				flag = false;
				RedrawWindow(hWnd, NULL, NULL, RDW_INVALIDATE);
			}
		}
		break;

	case WM_CLOSE:
		if (!saved) {
			if (MessageBox(hWnd, L"Do You Want To Save Changes", L"Save", MB_YESNO) == IDYES) {
				SaveFile();
			}
		}
		DestroyWindow(hWnd);
		break;

	case WM_DESTROY:
		PostQuitMessage(0);
		break;

	default:
		return DefWindowProc(hWnd, mcode, wp, lp);
	}
	return 0;
}
int APIENTRY WinMain(HINSTANCE hi, HINSTANCE, LPSTR cmd, int n)
{
	/*
	two types of string exist
	-ascii string (1 byte) "Graphics"
	-unicode (2 byte) L"Graphics"
	*/
	//first step define window class with its properties
	WNDCLASS wc;
	wc.hInstance = hi;
	wc.cbClsExtra = wc.cbWndExtra = 0;
	wc.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.lpszClassName = L"MyClass";		//window name
	wc.lpszMenuName = MAKEINTRESOURCE(IDR_MENU1);
	wc.hIcon = LoadIcon(GetModuleHandle(NULL), MAKEINTRESOURCE(IDI_ICON1));
	wc.style = CS_HREDRAW | CS_VREDRAW;
	//wc.lpfnWndProc = DefWindowProc;  //long pointer to function(lpfn) handles received messages
	wc.lpfnWndProc = WndProc;
	//second step regist the class
	RegisterClass(&wc);		//regist the class in the operating system


							//third step create window,show,update

							//HWND hwnd = CreateWindow(L"Edit", L"Hello Windows", WS_OVERLAPPEDWINDOW, 0, 0, 800, 600, NULL, NULL, hi, NULL);	//TextFieldWindow
							//HWND hwnd = CreateWindow(L"MyClass", L"Hello Windows", WS_OVERLAPPEDWINDOW, 0, 0, 800, 600, NULL, NULL, hi, NULL);
	HWND hwnd = CreateWindow(L"MyClass", L"Drawing", WS_OVERLAPPEDWINDOW, 50, 50, width, height, NULL, NULL, hi, NULL);
	ShowWindow(hwnd, n);
	UpdateWindow(hwnd);
	MSG msg;
	while (GetMessage(&msg, NULL, 0, 0) > 0)
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);		//call lpfnWndProc inside it, distribute messages to the one need it
	}
	return 0;
}