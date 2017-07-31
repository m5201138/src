#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <eigen/Eigen/Core>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
void Segmentation::init(){
int segmentNumber;
std::vector<Vec3> vec3Points;

std::vector<PointVectorPair> points;
PointList vertices;
//    std::vector<std::vector<Vector3> > pointsSets;
//  std::vector<std::vector<Vector3> > normalsSets;
//std::vector<Poisson_reconstruction_function> functionSets;
for(unsigned i=0;i<meshPtr->numVerts();i++){
    Vec3 p_neighbor = meshPtr->getVertPos(i);
    vec3Points.push_back(p_neighbor);
    
}
meshPtr->replacePoints(vec3Points);
meshPtr->makeMap();

std::set<Point> settmp;
for(int i=0;i<meshPtr->numSeg();i++){
    std::cout<<"meshPtr->numSeg()"<<i<<std::endl;
    std::pair<std::multimap<int,Point>::iterator, std::multimap<int,Point>::iterator> p = meshPtr->getEqual_range(i);
    for(std::multimap<int,Point>::iterator it = p.first;it!=p.second;it++){
        std::cout<<"segpoints[i]="<<it->second<<std::endl;
        
        settmp.insert(it->second);
    }
    // segPoints.push_back(settmp);
    
    
    for(auto itr=settmp.begin();itr!=settmp.end();itr++){
        points.push_back(std::make_pair(*itr, Vector(0,0,0)));
    }
    
    CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(), CGAL::First_of_pair_property_map<PointVectorPair>(), CGAL::Second_of_pair_property_map<PointVectorPair>(), nb_neighbors);
    
    std::vector<PointVectorPair>::iterator unoriented_points_begin =
    CGAL::mst_orient_normals(points.begin(), points.end(),
                             CGAL::First_of_pair_property_map<PointVectorPair>(),
                             CGAL::Second_of_pair_property_map<PointVectorPair>(),
                             nb_neighbors);
    std::vector<Vector3> points2;
    std::vector<Vector3> normals2;
    for(int j=0;j<points.size();j++){
        points2.push_back(Vector3(points[j].first.x(),points[j].first.y(),points[j].first.z()));
        normals2.push_back(Vector3(points[j].second.x(),points[j].second.y(),points[j].second.z()));
    }
    pointsSets.push_back(points2);
    normalsSets.push_back(normals2);
    points.erase(points.begin(),points.end());
    points2.erase(points2.begin(),points2.end());
    normals2.erase(normals2.begin(),normals2.end());
    settmp.erase(settmp.begin(),settmp.end());
    }
    
    std::vector<Vector3> pointsVector3;
    std::vector<Vector3> normalsVector3;
    for(int i=0;i<meshPtr->numSeg();i++){
        for(int j=0;j<pointsSets[i].size();j++){
            pointsVector3.push_back(pointsSets[i][j]);
            normalsVector3.push_back(normalsSets[i][j]);
        }
        std::cout<<"i="<<i<<std::endl;
        functionSets.push_back(Poisson_reconstruction(pointsVector3,normalsVector3,meshPtr));
        pointsVector3.erase(pointsVector3.begin(),pointsVector3.end());
        normalsVector3.erase(normalsVector3.begin(),normalsVector3.end());
    }
    // std::vector<Vector3> maxPoints;
    // std::vector<Vector3> maxNormals;
    for(int i=0;i<meshPtr->numSeg();i++){
        for(int j=0;j<pointsSets[i].size();j++){
            std::cout<<"func="<<functionSets[i](Point(pointsSets[i][j].x(),pointsSets[i][j].y(),pointsSets[i][j].z()))<<std::endl;
            if(functionSets[i](Point(pointsSets[i][j].x(),pointsSets[i][j].y(),pointsSets[i][j].z()))>=0){
                
                maxPoints.push_back(pointsSets[i][j]);
                maxNormals.push_back(normalsSets[i][j]);
            }
        }
    }
    std::cout<<"maxpoints="<<maxPoints.size()<<std::endl;
    
    std::vector<Point> points3;
    for(unsigned i=0; i<maxPoints.size(); i++){
        points3.push_back(Point(maxPoints[i].x(), maxPoints[i].y(),
                                maxPoints[i].z()));
    }
    
    originalFunction.push_back(Poisson_reconstruction(maxPoints,maxNormals,meshPtr));
}