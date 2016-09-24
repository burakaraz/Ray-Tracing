#include <iostream>
#include <cstdio>
#include "tinyXML/tinyxml.h"
#include "raytracer.h"
#include <map>

extern vector<Material *> vectorMaterial;
//extern vector<Vertex *> vectorVertex;
extern map<int,Vertex * > mapVertex;
extern vector<Triangle *> vectorTriangle;
extern vector<Sphere *> vectorSphere;
extern vector<PointLight *> vectorPointLight;
extern Intensity *ambientLight;
extern Intensity *backGroundColor;
extern int rayReflectionCount;
extern vector<Camera *> vectorCamera;

bool parseSceneXML(const char* filename)
{
	TiXmlDocument doc(filename);
	bool loadOkay = doc.LoadFile();

	if (!loadOkay)
	{
		std::cout << "Could not load file: " << filename << "Error = " << doc.ErrorDesc() << std::endl;
		return false;
	}

	TiXmlNode* pRoot = doc.FirstChild("Scene");
	for (TiXmlNode* pNode = pRoot->FirstChild(); pNode; pNode = pNode->NextSibling())
	{
        if (pNode->Value() == std::string("Material"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // get material index

            //
            // read reflectance coefficients
            //
            float amb[3], dif[3], spe[3], mir[3];
			float phongExp;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Ambient"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &amb[0], &amb[1], &amb[2]);
				}
				else if (pChild->Value() == std::string("Diffuse"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &dif[0], &dif[1], &dif[2]);
				}
				else if (pChild->Value() == std::string("Specular"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f %f", &spe[0], &spe[1], &spe[2], &phongExp);
				}
				else if (pChild->Value() == std::string("Reflectance"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &mir[0], &mir[1], &mir[2]);
				}
			}
			/**************************************************************/
			
			Material *newMaterial = new Material;
			Intensity *ambient = new Intensity;
			Intensity *diffuse = new Intensity;
			Intensity *specular = new Intensity;
			Intensity *reflectance = new Intensity;
			newMaterial -> mid = index;
			ambient -> r = amb[0];
			ambient -> g = amb[1];
			ambient -> b = amb[2];
			newMaterial -> ambient = ambient;
			diffuse -> r = dif[0];
			diffuse -> g = dif[1];
			diffuse -> b = dif[2];
			newMaterial -> diffuse = diffuse;
			specular -> r = spe[0];
			specular -> g = spe[1];
			specular -> b = spe[2];
			newMaterial -> specular = specular;
			newMaterial -> specExp = phongExp;
			reflectance -> r = mir[0];
			reflectance -> g = mir[1];
			reflectance -> b = mir[2];
			newMaterial -> reflectance = reflectance;
			vectorMaterial.push_back(newMaterial);

            /********************************************************************/
            //
            // TODO: Save the scanned values for the current material in your data structures.
            //

        }
		else if (pNode->Value() == std::string("Vertex"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // get vertex index

			float coords[3];
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->FirstChild()->Value(), "%f %f %f", &coords[0], &coords[1], &coords[2]);

			/********************************************************************/
            
            Vertex *newVertex = new Vertex;
            newVertex -> vid = index;
            
            Position *newCoordinates = new Position;
            newCoordinates -> x = coords[0];
            newCoordinates -> y = coords[1];
            newCoordinates -> z = coords[2];
            
            newVertex -> coordinates = newCoordinates;
            //vectorVertex.push_back(newVertex);
            mapVertex[index] = newVertex;
            
            /**********************************************************************/
            //
            // TODO: Save the scanned values for the current vertex in your data structures.
            //

		}
		else if (pNode->Value() == std::string("Triangle"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // get triangle index

			int vIndex[3], mIndex;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Vertices"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d %d %d", &vIndex[0], &vIndex[1], &vIndex[2]);
				}
				else if (pChild->Value() == std::string("MaterialId"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d", &mIndex);
				}
			}


			
			/****************************************************************************/
			Triangle *newTriangle = new Triangle;
			newTriangle -> tid = index;
			newTriangle -> mid = mIndex;
			
			Vertices *newVertices = new Vertices;
			newVertices -> a = vIndex[0];	
			newVertices -> b = vIndex[1];	
			newVertices -> c = vIndex[2];
			
			newTriangle -> vertice = newVertices;
			vectorTriangle.push_back(newTriangle);		
			
			/****************************************************************************/
            //
            // TODO: Save the scanned values for the current triangle in your data structures.
            //


		}
		else if (pNode->Value() == std::string("Sphere"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // Sphere index

			int vIndex, mIndex;
			float rad;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Center"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d", &vIndex);
				}
				else if (pChild->Value() == std::string("MaterialId"))
				{
					sscanf(pChild->FirstChild()->Value(), "%d", &mIndex);
				}
				else if (pChild->Value() == std::string("Radius"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f", &rad);
				}
			}
			
			/*******************************************************************/
			
			Sphere *newSphere = new Sphere;
			newSphere -> sid = index;
			newSphere -> vid = vIndex;
			newSphere -> mid = mIndex;
			newSphere -> radious = rad;
			
			vectorSphere.push_back(newSphere);	
			
			/*******************************************************************/

            //
            // TODO: Save the scanned values for the current sphere in your data structures.
            //

		}
		else if (pNode->Value() == std::string("PointLight"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // Light index

			float pos[3], intensity[3];
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Position"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &pos[0], &pos[1], &pos[2]);
				}
				else if (pChild->Value() == std::string("Intensity"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &intensity[0], &intensity[1], &intensity[2]);
				}
			}
			
			/**********************************************************************/
			PointLight *newPointLight = new PointLight;
			newPointLight -> plid = index;
			
			Position *newPos = new Position;
			Intensity *newIntens = new Intensity;
			
			newPos -> x = pos[0];
			newPos -> y = pos[1];
			newPos -> z = pos[2];
			
			newIntens -> r = intensity[0];
			newIntens -> g = intensity[1];
			newIntens -> b = intensity[2];
			
			newPointLight -> pos = newPos;
			newPointLight -> intens = newIntens;
			
			vectorPointLight.push_back(newPointLight);
			
			/***********************************************************************/

            //
            // TODO: Save the scanned values for the current point light in your data structures.
            //
		}
		else if (pNode->Value() == std::string("AmbientLight"))
        {
			float intensity[3];
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->Value(), "%f %f %f", &intensity[0], &intensity[1], &intensity[2]);


			/***************************************************************/
			
			ambientLight -> r = intensity[0];
			ambientLight -> g = intensity[1];
			ambientLight -> b = intensity[2];
			
			/***************************************************************/
            //
            // TODO: Save the scanned values for the ambient light in your data structures.
            //
		}
		else if (pNode->Value() == std::string("BackgroundColor"))
        {
            float bgColor[3];
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->Value(), "%f %f %f", &bgColor[0], &bgColor[1], &bgColor[2]);
			
			/****************************************************************/
			backGroundColor -> r = bgColor[0];
			backGroundColor -> g = bgColor[1];
			backGroundColor -> b = bgColor[2];
			
			/****************************************************************/

            //
            // TODO: Save the scanned values for the background color in your data structures.
            //
		}
		else if (pNode->Value() == std::string("RayReflectionCount"))
        {
            int rayReflectCount;
			TiXmlNode* pChild = pNode->FirstChild();
			sscanf(pChild->Value(), "%d", &rayReflectCount);

			/*****************************************************************/
			rayReflectionCount = rayReflectCount;
			/*****************************************************************/
            //
            // TODO: Save the scanned values for the ray reflection count in your data structures.
            //
		}
		else if (pNode->Value() == std::string("Camera"))
        {
			TiXmlAttribute* pAtt = pNode->ToElement()->FirstAttribute();
			int index = pAtt->IntValue(); // Camera index

			float gaze[3], up[3], pos[3];
			float left, right, bottom, top, distance;
			int nx, ny;
			std::string imageName;
			for (TiXmlNode* pChild = pNode->FirstChild(); pChild; pChild = pChild->NextSibling())
			{
				if (pChild->Value() == std::string("Position"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &pos[0], &pos[1], &pos[2]);
				}
				else if (pChild->Value() == std::string("Gaze"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &gaze[0], &gaze[1], &gaze[2]);
				}
				else if (pChild->Value() == std::string("Up"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f", &up[0], &up[1], &up[2]);
				}
				else if (pChild->Value() == std::string("ImagePlane"))
				{
					sscanf(pChild->FirstChild()->Value(), "%f %f %f %f %f %d %d", &left, &right, &bottom, &top, &distance, &nx, &ny);
				}
				else if (pChild->Value() == std::string("OutputName"))
				{
					imageName = pChild->FirstChild()->Value();
				}
			}

			/*****************************************************************/
			
			Camera *newCamera = new Camera;
			newCamera -> cid = index;
			
			Position *newPosition = new Position;
			Position *newGaze = new Position;
			Position *newUp = new Position;
			
			newPosition -> x = pos[0];
			newPosition -> y = pos[1];
			newPosition -> z = pos[2];
			
			newGaze -> x = gaze[0];
			newGaze -> y = gaze[1];
			newGaze -> z = gaze[2];
			
			newUp -> x = up[0];
			newUp -> y = up[1];
			newUp -> z = up[2];
			
			Plane *newPlane = new Plane;
			newPlane -> left = left;
			newPlane -> right = right;
			newPlane -> bottom = bottom;
			newPlane -> top = top;
			
			newCamera -> horRes = nx;
			newCamera -> verRes = ny;
			
			newCamera -> outputName = imageName;
			
			newCamera -> position = newPosition;
			newCamera -> gaze = newGaze;
			newCamera -> up = newUp;
			newCamera -> plane = newPlane;
			
			newCamera -> distance = distance;
			
			
			vectorCamera.push_back(newCamera);
			
			
			/*****************************************************************/
            //
            // TODO: Save the scanned values for the current camera in your data structures.
            //

		}
	}

    return true;
}
