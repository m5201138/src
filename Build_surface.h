#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <map>
#include <algorithm>
#include <iterator>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::FT FT;
//typedef CGAL::Triangulation_vertex_base_with_info_3<Point ,Kernel> Vb;
//typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
typedef Delaunay::Vertex_handle    Vertex_handle;
template <class HDS>
class Build_surface : public CGAL::Modifier_base<HDS> {
public:
    Build_surface() {
    }
    void addVertex(Point p){
        point.push_back(p);
        
    }
    void addFace(std::set<std::vector<int> > s){
        st=s;
        //copy(s.begin(), s.end(), back_inserter(face) );
    }
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        std::cout << point.size() << std::endl;
        std::cout << st.size() << std::endl;
        
        B.begin_surface( point.size(), st.size());
        for(int i= 0; i<point.size(); ++i) {
            std::cout << "adding: " << point[i] << std::endl;
            B.add_vertex(point[i]);
        }
        int i=0;
        std::set< std::vector<int> >::iterator itr;

        for(itr=st.begin(); itr!=st.end(); ++itr){
            std::cout << i << std::endl;
            i++;
            
            if (B.test_facet(itr->begin(), itr->end())) {
                B.add_facet(itr->begin(), itr->end());
            }
        }
        B.end_surface();
        
    }
private:
    std::vector<Point> point;
    std::vector<int> face;
    std::map<Point, int> mp;
    std::set<std::vector<int> > st;
};