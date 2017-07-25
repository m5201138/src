#include "TriMesh.h"
#include "Triangle.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>

#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::Point_3 Point;
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
            int i=0;
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

void TriMesh::createOFFFile(const std::string& outFileName)
{

    std::ofstream out(outFileName.c_str());
        out<<"OFF"<<std::endl;
    out<<mPoints.size()<<" "<<mFaceToVert.size()<<" "<<0<<std::endl;
    for (int i = 0; i<mPoints.size(); ++i) {
        out << mPoints[i].x << " " << mPoints[i].y << " " << mPoints[i].z <<std::endl;
    }
    for(int i=0;i<mFaceToVert.size();i++){
        out<<"3 "<<mFaceToVert[i][0]<<" "<<mFaceToVert[i][1]<<" "<<mFaceToVert[i][2]<<std::endl;
    }
    
    out.close();
}
    

void TriMesh::segmentation(void){
    normalize();
    //Polyhedron mesh;
    std::ifstream input("mesh.off");
    if ( !input || !(input >> mesh) || mesh.empty() ) {
        std::cerr << "Not a valid off file." << std::endl;
        //  return EXIT_FAILURE;
    }
    // create a property-map for segment-ids
    typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);
    // calculate SDF values and segment the mesh using default parameters.
    number_of_segments = CGAL::segmentation_via_sdf_values(mesh, segment_property_map);
    std::cout << "Number of segments: " << number_of_segments << std::endl;
    // print segment-ids
    for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
        std::cout << segment_property_map[facet_it] << " ";
    }
    
    //visualize segmentation
    std::cout << std::endl;
    std::ofstream out("result.off");
    out<<"OFF"<<std::endl;
    out<<mPoints.size()<<" "<<mFaceToVert.size()<<" "<<0<<std::endl;
    for(Polyhedron::Vertex_const_iterator vit=mesh.vertices_begin();vit!=mesh.vertices_end();vit++){
        out<<CGAL::to_double(vit->point().x())<<" "<<CGAL::to_double(vit->point().y())<<" "<<CGAL::to_double(vit->point().z())<<std::endl;
        
    }
    double color=1.0/(double)number_of_segments;
    int count=0;
    std::cout<<"color="<<number_of_segments<<std::endl;
    typedef typename Polyhedron::Vertex_const_iterator VCI;
    typedef typename Polyhedron::Facet_const_iterator FCI;
    typedef typename Polyhedron::Halfedge_around_facet_const_circulator HFCC;
    typedef CGAL::Inverse_index<VCI> Index;
    Index index(mesh.vertices_begin(), mesh.vertices_end());
    for (FCI fi = mesh.facets_begin();fi != mesh.facets_end();++fi){
        HFCC hc = fi->facet_begin();
        HFCC hc_end = hc;
        out<<"3";
        segmentNumbers.push_back(segment_property_map[fi]);
        do {
            out<<" ";
            out<<index[VCI(hc->vertex())];
                        /*segmentedPointMultimap.insert(std::make_pair(segment_property_map[fi],Point(hc->vertex()->point().x(),hc->vertex()->point().y(),hc->vertex()->point().z())));
             segmentedPointMap.insert(std::make_pair(Point(hc->vertex()->point().x(),hc->vertex()->point().y(),hc->vertex()->point().z()),segment_property_map[fi]));*/
            ++hc;
        } while(hc != hc_end);
        if(segment_property_map[fi]==0){
            out<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<1.0f;
            colors.push_back(1.0);
            colors.push_back(1.0);
            colors.push_back(1.0);
        }
        else if(segment_property_map[fi]==1){
            out<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<0.0f<<" "<<std::fixed<<0.0f;
            colors.push_back(1.0);
            colors.push_back(0.0);
            colors.push_back(0.0);
        }
        else if(segment_property_map[fi]==2){
            out<<" "<<std::fixed<<0.0f<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<0.0f;
            colors.push_back(0.0);
            colors.push_back(1.0);
            colors.push_back(0.0);
        }
        else if(segment_property_map[fi]==3){
            out<<" "<<std::fixed<<0.0f<<" "<<std::fixed<<0.0f<<" "<<std::fixed<<1.0f;
            colors.push_back(0.0);
            colors.push_back(0.0);
            colors.push_back(1.0);
        }
        else if(segment_property_map[fi]==4){
            out<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<0.0f<<" "<<std::fixed<<1.0f;
            colors.push_back(1.0);
            colors.push_back(0.0);
            colors.push_back(1.0);
        }
        else if(segment_property_map[fi]==5){
            out<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<0.0f;
            colors.push_back(1.0);
            colors.push_back(1.0);
            colors.push_back(0.0);
        }
        else if(segment_property_map[fi]==6){
            out<<" "<<std::fixed<<0.0f<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<1.0f;
            colors.push_back(0.0);
            colors.push_back(1.0);
            colors.push_back(1.0);
        }
        else if(segment_property_map[fi]==7){
            out<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<0.5f<<" "<<std::fixed<<0.0f;
            colors.push_back(1.0);
            colors.push_back(0.5);
            colors.push_back(0.0);
        }
        else if(segment_property_map[fi]==8){
            out<<" "<<std::fixed<<0.0f<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<0.5f;
            colors.push_back(0.0);
            colors.push_back(1.0);
            colors.push_back(0.5);
        }
        else if(segment_property_map[fi]==9){
            out<<" "<<std::fixed<<0.5f<<" "<<std::fixed<<1.0f<<" "<<std::fixed<<0.0f;
            colors.push_back(0.5);
            colors.push_back(1.0);
            colors.push_back(0.0);
        }
        out<<std::endl;
    }
    
}


void TriMesh::makeMap(){
    for (int i=0;i<numFaces();i++){
        for(int j=0;j<3;j++){
             segmentedPointMap.insert(std::make_pair(Point(mPoints[mFaceToVert[i][j]].x,mPoints[mFaceToVert[i][j]].y,mPoints[mFaceToVert[i][j]].z),segmentNumbers[i]));
            segmentedPointMultimap.insert(std::make_pair(segmentNumbers[i],Point(mPoints[mFaceToVert[i][j]].x,mPoints[mFaceToVert[i][j]].y,mPoints[mFaceToVert[i][j]].z)));

        }
    }
                  replacePolyhedron();
}

void TriMesh::replacePolyhedron(){
    mesh.clear();
    createOFFFile("meshForReplacePolyhedron.off");
    std::ifstream input("meshForReplacePolyhedron.off");
    if ( !input || !(input >> mesh) || mesh.empty() ) {
        std::cerr << "Not a valid off file." << std::endl;
    }
    
}




