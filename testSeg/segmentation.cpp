#include <iostream>
#include <fstream>
#include <eigen/Eigen/Core>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <utility>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef boost::tuple<Point_3,int>                           Point_and_int;
typedef CGAL::Random_points_in_cube_3<Point_3>              Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
Traits_base>                                              Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef Neighbor_search::Tree Tree;
typedef K_neighbor_search::Tree                             Tree;
typedef K_neighbor_search::Distance                         Distance;

void build_graph(std::vector<Point> points, int k){
            std::vector<int>     indices;
    for(int i=0;i<points.size();i++){
        indices.push_back(i);
    }
    Tree tree(
              boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
              boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) )
              );
    Distance tr_dist;
    MatrixXd idx;
    idx.resize(points.size(),)
    for(int i=0;i<points.size();i++;){
            K_neighbor_search search(tree, points[i], k);
        for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
            
            << tr_dist.inverse_of_transformed_distance(it->second) << " "
            << boost::get<0>(it->first)<< " " << boost::get<1>(it->first) << std::endl;
        }

    }
    Neighbor_search search(tree, points[i], N);
    for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it){
            std::cout << it->first << " "<< std::sqrt(it->second) << std::endl;
        }
    }
    
}


void conv_clustering(std::vector<Point> points, std::vector<Vector> normals, int k, int pnum){
    //build_graph(points, k);
}
void initseg(std::vector<Point> points, std::vector<Vector> normals){
    conv_clustering(points,normals,10,75);
}
int main(){
    std::string line;
    std::ifstream in("Figure1_Deer.ply");
    std::vector<Point> points;
    std::vector<Vector> normals;
    while(getline(in, line))
    {
        std::stringstream ss(line);
        std::string token;
        double p1,p2,p3,n1,n2,n3;
        ss >> p1 >> p2 >> p3 >> n1 >> n2 >> n3;
        points.push_back(Point(p1,p2,p3));
        normals.push_back(Vector(n1,n2,n3));
    }
    initseg(points,normals);
    
    
}