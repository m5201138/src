#ifndef TRIMESH_H
#define TRIMESH_H

#include "types.h"
#include <vector>
#include "Vectors.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <CGAL/Polyhedron_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

class TriMesh
{
public:

  TriMesh() { }

  TriMesh(const TriMesh& mesh)
    : mPoints(mesh.mPoints)
    , mFaceToVert(mesh.mFaceToVert)
  { }

  TriMesh& operator = (const TriMesh& mesh)
  {
    mPoints = mesh.mPoints;
    mFaceToVert = mesh.mFaceToVert;
    return *this;
  }
    std::vector<double> colors;
  void clear();
  void normalize();
  void read(const char* filename);

  inline unsigned numVerts() const { return mPoints.size(); }
  inline unsigned numFaces() const { return mFaceToVert.size(); }
  inline const Vec3& getVertPos(unsigned vert) const { return mPoints[vert]; }
 std::pair<std::multimap<int,Point>::iterator, std::multimap<int,Point>::iterator> getEqual_range(int i){return segmentedPointMultimap.equal_range(i);}

  void getFaceVerts(unsigned face, std::vector<unsigned>& verts) const;
  unsigned getVertIndexInFace(unsigned face, unsigned vert) const;
  unsigned getFaceVert(unsigned face, unsigned index) const;
  double computeFaceArea(unsigned face) const;
  Vec3 computeFaceNormal(unsigned face) const;
  double computeTotalArea() const;
  double computeMaxEdgeLength() const;
  void computeVertNormals(std::vector<Vec3>& vertNormal) const;
  double computeEdgeLength(unsigned vert0, unsigned vert1) const;

  void setData(POINT3D *m_ppt3dVertices ,int m_nVertices,unsigned int *m_piTriangleIndices,int m_nTriangles);

  void setData(std::vector<Vec3>& points,
               std::vector< std::vector<unsigned> >& faces);
        void createOFFFile(const std::string& outFileName);
    void segmentation(void);
    int returnNumber_of_segments(){return number_of_segments;}
    int Numcolors(){return colors.size();}
    int getSegmentNumberForPoint(Point p){return segmentedPointMap[p];}
    void replacePoints(std::vector<Vec3>& points){mPoints=points;}
    void makeMap();

private:
    
  std::vector<Vec3> mPoints;
  std::vector< std::vector<unsigned> > mFaceToVert;
       std::vector<int> f;
   std::multimap<int,Point> segmentedPointMultimap;
    std::map<Point,int> segmentedPointMap;
    std::size_t number_of_segments;
    std::vector<int> segmentNumbers;

};

#endif // TRIMESH_H

