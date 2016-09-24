#include <iostream>
#include "parseScene.h"
#include "writePPM.h"
#include "raytracer.h"
#include <vector>
#include <math.h>
#include <limits>
#include <map>


Intensity rayColor(Intensity color, const Position &o ,const Position &d,const Position &v) ;

vector<Material *> vectorMaterial;
map<int,Vertex * > mapVertex;
vector<Triangle *> vectorTriangle;
vector<Sphere *> vectorSphere;
vector<PointLight *> vectorPointLight;
vector<Camera *> vectorCamera;

Intensity *ambientLight = new Intensity;
Intensity *backGroundColor = new Intensity;


float tmin;
float* data;
int rayReflectionCount;
int countRay = rayReflectionCount;


float dot(Position p1, Position p2)
{
	return (p1.x * p2.x) + (p1.y * p2.y) + (p1.z * p2.z);
}

bool sgnTest(float num1, float num2, float num3)
{
		if (num1 > 0 && num2 > 0 && num3 > 0)
		{
			return true;
		}
		else if(num1 < 0 && num2 < 0 && num3 < 0)
		{
			return true;
		}
		else 
		{
			return false;
		}
}

Intensity compoMult(const Intensity &p1,const Intensity &p2)
{
	Intensity result;
	result.r = p1.r + p2.r;
	result.g = p1.g + p2.g;
	result.b = p1.b + p2.b;
	return result;
}

Position addition(const Position &p1,const Position &p2)
{
	Position result;
	result.x = p1.x + p2.x;
	result.y = p1.y + p2.y;
	result.z = p1.z + p2.z;
	return result;
}

Intensity additionInt(const Intensity &p1,const Intensity &p2)
{
	Intensity result;
	result.r = p1.r + p2.r;
	result.g = p1.g + p2.g;
	result.b = p1.b + p2.b;
	return result;
}

Position makeNormalized(const Position &p1)
{
	Position Np1;
	float total = sqrt(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);
	Np1.x = p1.x / total;
	Np1.y = p1.y / total;
	Np1.z = p1.z / total;
	return Np1;
}

float maxN(float num1, float num2)
{
	if (num1 > num2)
	{
		return num1;
	}
	else
	{
		return num2;
	}
}

Intensity shading(Intensity color ,const Position &hitPoint,const Position &normal,const Material &material, const Position &camera)
{
	color.r = ambientLight -> r * material.ambient -> r;
	color.g = ambientLight -> g * material.ambient -> g;
	color.b = ambientLight -> b * material.ambient -> b;
	Position o;
	Position v;
	
	v = camera - hitPoint;
	v = makeNormalized(v);
		
	for(int i = 0 ; i < vectorPointLight.size() ; i++)
	{
		Position lh;
		lh.x = vectorPointLight[i]->pos->x - hitPoint.x;
		lh.y = vectorPointLight[i]->pos->y - hitPoint.y;
		lh.z = vectorPointLight[i]->pos->z - hitPoint.z;
		
		o.x = hitPoint.x + lh.x * 0.001;
		o.y = hitPoint.y + lh.y * 0.001;
		o.z = hitPoint.z + lh.z * 0.001;
		
		int flagSh = 0;
		for(int j = 0 ; j < vectorSphere.size() ; j++)
		{	
			Position c;
			
			c.x = vectorSphere[j] -> center -> x;
			c.y = vectorSphere[j] -> center -> y;
			c.z = vectorSphere[j] -> center -> z;
			
			float A = lh.x * lh.x + lh.y * lh.y + lh.z * lh.z;

			Position OminusC;
			OminusC = (o) - c;
					
			float B = 2 * ((lh.x) *  (OminusC.x) + (lh.y) *  (OminusC.y) + (lh.z) *  (OminusC.z));
			float C = (OminusC.x * OminusC.x + OminusC.y * OminusC.y + OminusC.z * OminusC.z) - (vectorSphere[j] -> radious) * (vectorSphere[j] -> radious);
					
			float delta = (B * B - (4 * A * C));
			//cout << delta << endl;
			if (delta >= 0)
			{
				// calculate t
				float tSphere1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
				float tSphere2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
						
				if(tSphere1 < 1 && tSphere1 > 0)
				{
					flagSh = 1;
					break;
							
				}
				if (tSphere2 < 1 && tSphere2 > 0)
				{
					flagSh = 1;
					break;
				}
			}	
		}
		if (flagSh == 1)
		{
			continue;
		}
		Position aVertice;
		Position bVertice;
		Position cVertice;
		for(int l = 0 ; l < vectorTriangle.size() ; l++)
		{
			aVertice.x = vectorTriangle[l]->aVert->x;
			aVertice.y = vectorTriangle[l]->aVert->y;
			aVertice.z = vectorTriangle[l]->aVert->z;
			
			bVertice.x = vectorTriangle[l]->bVert->x;
			bVertice.y = vectorTriangle[l]->bVert->y;
			bVertice.z = vectorTriangle[l]->bVert->z;
			
			cVertice.x = vectorTriangle[l]->cVert->x;
			cVertice.y = vectorTriangle[l]->cVert->y;
			cVertice.z = vectorTriangle[l]->cVert->z;
			
					
			Position normalizedN;

			normalizedN.x =  vectorTriangle[l]->normalN->x;
			normalizedN.y =  vectorTriangle[l]->normalN->y;
			normalizedN.z =  vectorTriangle[l]->normalN->z;
					
			Position aminuso;
			
			aminuso = (aVertice - o);
			
			float tTriangle = (dot(aminuso, normalizedN)) / (dot(lh,normalizedN));

					
			Position findedDot;
			Position v1;
			Position v2;
			Position v3;
					
			findedDot.x = o.x + tTriangle * lh.x;  
			findedDot.y = o.y + tTriangle * lh.y;  
			findedDot.z = o.z + tTriangle * lh.z;  
					
			v1 = (bVertice - aVertice) * (findedDot - aVertice);
			v2 = (cVertice - bVertice) * (findedDot - bVertice);
			v3 = (aVertice - cVertice) * (findedDot - cVertice);
			
			if(dot(v1,v2) > -0.000000000000000001 && dot(v2,v3) > -0.000000000000000001 && tTriangle < 1 && tTriangle > 0) 
			{
				flagSh = 1;
				break;
			}
		}
		if (flagSh == 1)
		{
			continue;
		}
		Position normalL;
		normalL = makeNormalized(lh);

		color.r = color.r + ( material.diffuse -> r * vectorPointLight[i] -> intens -> r * maxN(0,dot(normalL,normal)) / (4 * M_PI * dot(lh,lh)));
		color.g = color.g + ( material.diffuse -> g * vectorPointLight[i] -> intens -> g * maxN(0,dot(normalL,normal)) / (4 * M_PI * dot(lh,lh)));
		color.b = color.b + ( material.diffuse -> b * vectorPointLight[i] -> intens -> b * maxN(0,dot(normalL,normal)) / (4 * M_PI * dot(lh,lh)));	
		
		Position h;
		h = addition(normalL,v);
		Position normalh;
		normalh = makeNormalized(h);

		color.r = color.r + (material.specular -> r * vectorPointLight[i] -> intens -> r * pow(maxN(0,dot(normal,normalh)),material.specExp)) / (4 * M_PI * dot(lh,lh));
		color.g = color.g + (material.specular -> g * vectorPointLight[i] -> intens -> g * pow(maxN(0,dot(normal,normalh)),material.specExp)) / (4 * M_PI * dot(lh,lh));
		color.b = color.b + (material.specular -> b * vectorPointLight[i] -> intens -> b * pow(maxN(0,dot(normal,normalh)),material.specExp)) / (4 * M_PI * dot(lh,lh));
		
		
	}//cout << material.reflectance->r << " " << material.reflectance->g << " " << material.reflectance->b << endl;
	if((material.reflectance->r != 0 || material.reflectance->g != 0 || material.reflectance->b != 0) && countRay > 0)
	{
		
		Position r;
		float datVN = dot(v,normal);
		r.x = (2 * datVN) * normal.x - v.x;
		r.y = (2 * datVN) * normal.y - v.y;
		r.z = (2 * datVN) * normal.z - v.z;
		
		r = makeNormalized(r);
		
		Position newHitP;
		newHitP.x = hitPoint.x + r.x * 0.01;
		newHitP.y = hitPoint.y + r.y * 0.01;
		newHitP.z = hitPoint.z + r.z * 0.01;
		
		countRay--;
		return additionInt(color,compoMult(*material.reflectance,rayColor(color, newHitP, r,v)));
	}
	
	return color;
	
}

Intensity rayColor(Intensity color, const Position &o ,const Position &d,const Position &v) 
{
	int flag = 0 ;
	Position aVertice;
	Position bVertice;
	Position cVertice;
	
	Position minCenterSphere;
	Position minNormalTriangle;
				
	Position c;
	
	Triangle *pTriangleObj = NULL;
	Sphere *pSphereObj = NULL;
	
	tmin = std::numeric_limits<float>::max();
	
				
	for(int l = 0 ; l < vectorTriangle.size() ; l++)
	{

		aVertice.x = vectorTriangle[l]->aVert->x;
		aVertice.y = vectorTriangle[l]->aVert->y;
		aVertice.z = vectorTriangle[l]->aVert->z;
		
		bVertice.x = vectorTriangle[l]->bVert->x;
		bVertice.y = vectorTriangle[l]->bVert->y;
		bVertice.z = vectorTriangle[l]->bVert->z;
		
		cVertice.x = vectorTriangle[l]->cVert->x;
		cVertice.y = vectorTriangle[l]->cVert->y;
		cVertice.z = vectorTriangle[l]->cVert->z;
		
		Position normalizedN;

		normalizedN.x =  vectorTriangle[l]->normalN->x;
		normalizedN.y =  vectorTriangle[l]->normalN->y;
		normalizedN.z =  vectorTriangle[l]->normalN->z;
		
		//cout << normalizedN.x << " " << normalizedN.y << " " << normalizedN.z << endl;
		Position aminuso;
		aminuso = (aVertice - o);
			
		float tTriangle = (dot(aminuso, normalizedN)) / (dot(d,normalizedN));
		
		//cout << tTriangle << endl;
		
		Position findedDot;
		Position v1;
		Position v2;
		Position v3;
		
		findedDot.x = o.x + tTriangle * d.x;  
		findedDot.y = o.y + tTriangle * d.y;  
		findedDot.z = o.z + tTriangle * d.z;  
		
		v1 = (bVertice - aVertice) * (findedDot - aVertice);
		v2 = (cVertice - bVertice) * (findedDot - bVertice);
		v3 = (aVertice - cVertice) * (findedDot - cVertice);
		if(dot(v1,v2) > -0.000000000000000001 && dot(v2,v3) > -0.000000000000000001)
		{
			if (tTriangle < tmin && tTriangle > 0)
			{
				pTriangleObj = vectorTriangle[l];
				tmin = tTriangle;
				flag = 1;
				minNormalTriangle = normalizedN;
				
				//cout << minNormalTriangle -> x << " " << minNormalTriangle->y << " " << minNormalTriangle->z << endl;
			}
		}
	}
	
	for(int l = 0 ; l < vectorSphere.size() ; l++)
	{
		
		c.x = vectorSphere[l] -> center -> x;
		c.y = vectorSphere[l] -> center -> y;
		c.z = vectorSphere[l] -> center -> z;
		
		float A = d.x * d.x + d.y * d.y + d.z * d.z;

		Position OminusC;
		OminusC = (o) - c;
		
		float B = 2 * ((d.x) *  (OminusC.x) + (d.y) *  (OminusC.y) + (d.z) *  (OminusC.z));
		float C = (OminusC.x * OminusC.x + OminusC.y * OminusC.y + OminusC.z * OminusC.z) - (vectorSphere[l] -> radious) * (vectorSphere[l] -> radious);
		
		float delta = (B * B - (4 * A * C));
		
		if (delta >= 0)
		{
			// calculate t
			float tSphere1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
			float tSphere2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
			
			if(tSphere1 < tmin && tSphere1 > 0)
			{
				pSphereObj = vectorSphere[l];
				tmin = tSphere1;
				flag = 2;
				minCenterSphere.x = c.x;
				minCenterSphere.y = c.y;
				minCenterSphere.z = c.z;
				
			}
			if (tSphere2 < tmin && tSphere2 > 0)
			{
				pSphereObj = vectorSphere[l];
				tmin = tSphere2;
				flag = 2;
				minCenterSphere.x = c.x;
				minCenterSphere.y = c.y;
				minCenterSphere.z = c.z;
			}
		}

	}

	if (pSphereObj == NULL && pTriangleObj == NULL)
	{
		color.r = backGroundColor->r;
		color.g = backGroundColor->g;
		color.b = backGroundColor->b;
	}
	else
	{
		Position hitPoint; 
		Position tmind;
		
		if(flag == 1) /// triangle
		{
			Material matTriangle;
			for(int c = 0 ; c < vectorMaterial.size() ; c++)
			{
				if (vectorMaterial[c]-> mid == pTriangleObj->mid)
				{
					matTriangle = *vectorMaterial[c];
				}
			}
			hitPoint.x = o.x + tmin * d.x;  
			hitPoint.y = o.y + tmin * d.y;  
			hitPoint.z = o.z + tmin * d.z; 
				
			tmind.x = tmin * d.x;
			tmind.y = tmin * d.y;
			tmind.z = tmin * d.z;
			color = shading(color, hitPoint, minNormalTriangle, matTriangle, o);
			
		}
		else        /// sphere
		{
			Material matSphere;
			tmind.x = tmin * d.x;
			tmind.y = tmin * d.y;
			tmind.z = tmin * d.z;
			for(int c = 0 ; c < vectorMaterial.size() ; c++)
			{
				if (vectorMaterial[c]-> mid == pSphereObj->mid)
				{	
					matSphere = *vectorMaterial[c];
				}
			}
			hitPoint.x = o.x + tmin * d.x;  
			hitPoint.y = o.y + tmin * d.y;  
			hitPoint.z = o.z + tmin * d.z;
			
			
			Position normalSphere;
			normalSphere.x = (hitPoint.x - minCenterSphere.x) / pSphereObj -> radious;
			normalSphere.y = (hitPoint.y - minCenterSphere.y) / pSphereObj -> radious;
			normalSphere.z = (hitPoint.z - minCenterSphere.z) / pSphereObj -> radious;	
			
			Position normalizedSphere;
			normalizedSphere = makeNormalized(normalSphere);
	
			color = shading(color, hitPoint, normalizedSphere, matSphere, o);
		}
		
	}
	return color;
}

vector<vector<Intensity> > rayTracing(int i)
{
	vector<vector<Intensity> > color;
	
	int horRes = vectorCamera[i]->horRes;
	int verRes = vectorCamera[i]->verRes;
	
	color.resize(verRes);
	for(int w = 0 ; w < verRes ; w++)
	{
		color[w].resize(horRes);
	}
	
	float distance = vectorCamera[i]->distance;
	
	Position w;
	Position v;
	Position o;
	o.x = vectorCamera[i]->position->x;
	o.y = vectorCamera[i]->position->y;
	o.z = vectorCamera[i]->position->z;
	
	
	w.x = -(vectorCamera[i]->gaze->x); 				// w 
	w.y = -(vectorCamera[i]->gaze->y); 
	w.z = -(vectorCamera[i]->gaze->z); 
	
	w = makeNormalized(w);
	
	Position u;							// u up -> v
	u = (*vectorCamera[i]->up) * (w);
	
	u = makeNormalized(u);
	v = makeNormalized(*vectorCamera[i]->up);
	
	for (int j = 0 ; j < verRes ; j++)
	{
		for (int k = 0 ; k < horRes; k++)
		{
			float su = vectorCamera[i]->plane->left + (float)(vectorCamera[i]->plane->right - vectorCamera[i]->plane->left) * (k + 0.5) / (float)horRes;
			float sv = vectorCamera[i]->plane->top - (float)(vectorCamera[i]->plane->top - vectorCamera[i]->plane->bottom) * (j + 0.5) / (float) verRes;
			
			Position d;
			d.x = su * u.x + sv * v.x - w.x * distance; 
			d.y = su * u.y + sv * v.y - w.y * distance; 
			d.z = su * u.z + sv * v.z - w.z * distance;  
			
			//cout << d.x << " " << d.y << " " << d.z << endl;
			//////////////////////////////////////////////////////////////////////////////////
			countRay = rayReflectionCount;
			color[j][k] = rayColor(color[j][k] ,o ,d ,v);
			
			//////////////////////////////////////////////////////////////
			
		}
	}
	
	return color;
}

void preProcess()
{
	for(int i = 0 ; i < vectorTriangle.size() ; i++)
	{
		Position *aVert = new Position;
		aVert->x = mapVertex[vectorTriangle[i] -> vertice -> a] -> coordinates -> x;
		aVert->y = mapVertex[vectorTriangle[i] -> vertice -> a] -> coordinates -> y;
		aVert->z = mapVertex[vectorTriangle[i] -> vertice -> a] -> coordinates -> z;
		
		Position *bVert = new Position;
		bVert->x = mapVertex[vectorTriangle[i] -> vertice -> b] -> coordinates -> x;
		bVert->y = mapVertex[vectorTriangle[i] -> vertice -> b] -> coordinates -> y;
		bVert->z = mapVertex[vectorTriangle[i] -> vertice -> b] -> coordinates -> z;
		
		Position *cVert = new Position;
		cVert->x = mapVertex[vectorTriangle[i] -> vertice -> c] -> coordinates -> x;
		cVert->y = mapVertex[vectorTriangle[i] -> vertice -> c] -> coordinates -> y;
		cVert->z = mapVertex[vectorTriangle[i] -> vertice -> c] -> coordinates -> z;
		
		vectorTriangle[i]->aVert = aVert;
		vectorTriangle[i]->bVert = bVert;
		vectorTriangle[i]->cVert = cVert;
		
		Position *normalN = new Position;
		
		Position norm;
		norm = ((*bVert - *aVert) * (*cVert - *aVert));
		norm = makeNormalized(norm);
		
		normalN->x = norm.x;
		normalN->y = norm.y;
		normalN->z = norm.z;
		
		vectorTriangle[i]->normalN = normalN;
	}
	for(int i = 0 ; i < vectorSphere.size() ; i++)
	{
		Position *c = new Position;
		c->x = mapVertex[vectorSphere[i] -> vid] -> coordinates -> x;
		c->y = mapVertex[vectorSphere[i] -> vid] -> coordinates -> y;
		c->z = mapVertex[vectorSphere[i] -> vid] -> coordinates -> z;
		
		vectorSphere[i] -> center = c;
	}
}

int main(int argc, char* argv[])
{

    bool result = parseSceneXML(argv[1]);
    preProcess();
    
    for(int k = 0 ; k < vectorCamera.size() ; k++)
    {
		int horRes = vectorCamera[k]->horRes;
		int verRes = vectorCamera[k]->verRes;
		
		data = new float[horRes * verRes * 3];
		vector<vector<Intensity> > pColor = rayTracing(k);
		int index = 0;
		for(int j = 0 ; j < verRes ; j++)
		{
			for (int i = 0 ; i < horRes ; i++)
			{
				if(pColor[j][i].r > 255)
					data[index] = 255;
				else
					data[index] = pColor[j][i].r;
				if(pColor[j][i].g > 255)
					data[index+1] = 255;
				else
					data[index+1] = pColor[j][i].g;
				if(pColor[j][i].b > 255)
					data[index+2] = 255;
				else
					data[index+2] = pColor[j][i].b;
				index = index + 3; 
			}
		}
		
		writePPM(vectorCamera[k]->outputName.c_str(), horRes, verRes, data);
		delete[] data;
	}
	
	delete ambientLight;
	delete backGroundColor;
	for(int i = 0 ; i < vectorCamera.size() ; i++)
    {
		delete vectorCamera[i]->position;
		delete vectorCamera[i]->gaze;
		delete vectorCamera[i]->up;
		delete vectorCamera[i]->plane;
		delete vectorCamera[i];
	}
	for(int i = 0 ; i < vectorPointLight.size() ; i++)
    {
		delete vectorPointLight[i]->pos;
		delete vectorPointLight[i]->intens;
		delete vectorPointLight[i];
	}
	for(int i = 0 ; i < vectorSphere.size() ; i++)
    {
		delete vectorSphere[i]->center;
		delete vectorSphere[i];
	}
	for(int i = 0 ; i < vectorTriangle.size() ; i++)
    {
		delete vectorTriangle[i]->aVert;
		delete vectorTriangle[i]->bVert;
		delete vectorTriangle[i]->cVert;
		delete vectorTriangle[i]->normalN; 
		delete vectorTriangle[i]->vertice;
		delete vectorTriangle[i];
	}
	for (std::map<int,Vertex *>::iterator it=mapVertex.begin(); it!=mapVertex.end(); ++it)
	{
		delete it->second->coordinates;
		delete it->second;
	}
	for(int i = 0 ; i < vectorMaterial.size() ; i++)
    {
		delete vectorMaterial[i]->ambient;
		delete vectorMaterial[i]->diffuse;
		delete vectorMaterial[i]->specular;
		delete vectorMaterial[i]->reflectance;
		delete vectorMaterial[i];
	}
    return 0;
}

Position operator* (const Position& v1, const Position& v2)
{
	Position crossresult;
	crossresult.x = (v1.y * v2.z) - (v1.z * v2.y);
	crossresult.y = (v1.z * v2.x) - (v1.x * v2.z);
	crossresult.z = (v1.x * v2.y) - (v1.y * v2.x);
	return crossresult;
}

Position operator- (const Position& v1, const Position& v2)
{
	Position crossresult;
	crossresult.x = v1.x - v2.x;
	crossresult.y = v1.y - v2.y;
	crossresult.z = v1.z - v2.z;
	return crossresult;
}










