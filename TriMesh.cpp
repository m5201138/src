#include "TriMesh.h"
#include "Triangle.h"

#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>


void TriMesh::getFaceVerts(unsigned face, std::vector<unsigned>& verts) const
{
  verts = mFaceToVert[face];
}

unsigned TriMesh::getVertIndexInFace(unsigned face, unsigned vert) const
{
  std::vector<unsigned> fVert;
  getFaceVerts(face, fVert);
  assert(fVert.size() == 3);

  for (unsigned i = 0; i < 3; ++i) 
    if (fVert[i] == vert) return i;

  assert(false);
  return 0;    
}

unsigned TriMesh::getFaceVert(unsigned face, unsigned index) const
{
  return mFaceToVert[face][index];
}

double TriMesh::computeFaceArea(unsigned face) const
{
  std::vector<unsigned> fVert;
  getFaceVerts(face, fVert);
  assert(fVert.size() == 3);

  Vec3 p0 = getVertPos(fVert[0]);
  Vec3 p1 = getVertPos(fVert[1]);
  Vec3 p2 = getVertPos(fVert[2]);
  return Triangle::area(p0, p1, p2);
}

Vec3 TriMesh::computeFaceNormal(unsigned face) const
{
  std::vector<unsigned> fVert;
  getFaceVerts(face, fVert);
  assert(fVert.size() == 3);

  Vec3 p0 = getVertPos(fVert[0]);
  Vec3 p1 = getVertPos(fVert[1]);
  Vec3 p2 = getVertPos(fVert[2]);
  return (p1-p0).cross(p2-p1).normalized();
}

double TriMesh::computeTotalArea() const
{
  double sum = 0.0;
  for (unsigned face = 0; face < numFaces(); ++face) 
    sum += computeFaceArea(face);
  return sum;
}

double TriMesh::computeMaxEdgeLength() const
{
  double maxLen = 0.0;
  for (unsigned face = 0; face < numFaces(); ++face) 
    {
      std::vector<unsigned> fVerts;
      getFaceVerts(face, fVerts);
      assert(fVerts.size() == 3);

      for (unsigned i = 0; i < 3; ++i) 
        {
	  unsigned vert0 = fVerts[(i+1)%3];
	  unsigned vert1 = fVerts[(i+2)%3];
	  double len = computeEdgeLength(vert0, vert1);
	  maxLen = std::max(maxLen, len);
        }
    }
  return maxLen;
}

void TriMesh::computeVertNormals(std::vector<Vec3>& vertNormal) const
{
  vertNormal = std::vector<Vec3>(numVerts(), Vec3::Zero());
    
  for (unsigned face = 0; face < numFaces(); ++face) 
    {
      std::vector<unsigned> fVert;
      getFaceVerts(face, fVert);
      assert(fVert.size() == 3);

      double faceArea = computeFaceArea(face);
      Vec3 faceNormal = computeFaceNormal(face);
      for (unsigned i = 0; i < 3; ++i) 
	vertNormal[fVert[i]] += faceArea * faceNormal;
    }

  for (unsigned vert = 0; vert < numVerts(); ++vert) 
    {
      vertNormal[vert].normalize();
    }
}

double TriMesh::computeEdgeLength(unsigned vert0, unsigned vert1) const
{
  Vec3 p0 = getVertPos(vert0);
  Vec3 p1 = getVertPos(vert1);
  return (p1-p0).norm();
}

void TriMesh::clear()
{
  mPoints.clear();
  mFaceToVert.clear();
}

void TriMesh::read(const char* filename)
{
  clear();
  std::string line;
  std::ifstream in(filename);
  while(getline(in, line))
    {
      std::stringstream ss(line);
      std::string token;
      ss >> token;

      if (token == "v") 
        {
	  double x, y, z;
	  ss >> x >> y >> z;
	  mPoints.push_back(Vec3(x,y,z));
	  continue;
        }

      if (token == "f") 
        {
	  std::vector<unsigned> face;
	  while (ss >> token) 
            {
	      unsigned index;
	      std::string indexstring;
	      std::stringstream tokenstream(token);
	      getline(tokenstream, indexstring, '/');
	      std::stringstream indexstream(indexstring);
	      indexstream >> index;
	      face.push_back(index-1);
            }
	  mFaceToVert.push_back(face);
        }
    }
  in.close();
}

void TriMesh::normalize()
{
  Vec3 c = Vec3::Zero();
  for (unsigned vert = 0; vert < numVerts(); ++vert) {
    c += mPoints[vert];
  }

  c /= numVerts();
  for (unsigned vert = 0; vert < numVerts(); ++vert) {
    mPoints[vert] -= c;
  }

  double scale = std::sqrt(computeTotalArea());
  for (unsigned vert = 0; vert < numVerts(); ++vert) {
    mPoints[vert] /= scale;
  }
}

void TriMesh::setData(POINT3D *m_ppt3dVertices ,int m_nVertices,unsigned int *m_piTriangleIndices,int m_nTriangles){

  mPoints.clear();
  mFaceToVert.clear();

  double x,y,z;
  for(int i=0;i<m_nVertices;i++){
    x=(double)m_ppt3dVertices[i][0];
    y=(double)m_ppt3dVertices[i][1];
    z=(double)m_ppt3dVertices[i][2];
    mPoints.push_back(Vec3(x,y,z));
  }
   
  for(int i=0;i<m_nTriangles;i++){
    std::vector<unsigned> face;
    face.push_back(m_piTriangleIndices[i*3]);
    face.push_back(m_piTriangleIndices[i*3+1]);
    face.push_back(m_piTriangleIndices[i*3+2]);
    mFaceToVert.push_back(face);
    // std::cout<<mFaceToVert[i][0]<<std::endl;
  }
    std::cout<<"number of points: "<<mPoints.size()<<std::endl;
    std::cout<<"number of faces: "<<mFaceToVert.size()<<std::endl;
}


void TriMesh::setData(std::vector<Vec3>& points, 
		      std::vector< std::vector<unsigned> >& faces)
{
  mPoints = points;
  mFaceToVert = faces;
}
