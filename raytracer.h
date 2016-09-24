#include <iostream>
#include <vector>

using namespace std;

class Camera;
class Position;
class Intensity;
class Material;
class Vertex;
class Triangle;
class Vertices;
class Sphere;
class PointLight;
class Plane;


class Material
{
	public:
		int mid;
		Intensity *ambient;
		Intensity *diffuse;
		Intensity *specular;
		float specExp;
		Intensity *reflectance;
};

class Vertex
{
	public:
		int vid;
		Position *coordinates;
};

class Vertices
{
	public:
		int a;
		int b;
		int c;
};

class Triangle
{
	public:
		int tid;
		Vertices *vertice;
		int mid;
		Position *normalN;
		Position *aVert;
		Position *bVert;
		Position *cVert;
		
};

class Sphere
{
	public:
		int sid;
		int vid;
		int mid;
		float radious;
		Position *center;
};

class PointLight
{
	public:
		int plid;
		Position *pos;
		Intensity *intens;
};


class Position
{
	public:
		float x;
		float y;
		float z;
		friend Position operator* (const Position& v1, const Position& v2);
		friend Position operator- (const Position& v1, const Position& v2);
		void dot(Position* p1, Position* p2);
};


class Intensity
{
	public:
		float r;
		float g;
		float b;
};

class Plane
{
	public:
		float left;
		float right;
		float bottom;
		float top;
};

class Camera
{
	public:
		int cid;
		Position *position;
		Position *gaze;
		Position *up;
		Plane *plane;
		float distance;
		int horRes;
		int verRes;
		string outputName;

};
