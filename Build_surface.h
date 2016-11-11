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
typedef CGAL::Triangulation_vertex_base_with_info_3<Point ,Kernel> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds> Delaunay;
typedef Delaunay::Vertex_handle    Vertex_handle;
template <class HDS>
class Build_surface : public CGAL::Modifier_base<HDS> {
public:
    Build_surface() {
    }
    void addVertex(Point p){
        point.push_back(p);
        
    }
    void addFace(std::set<int> s){
        st=s;
        //copy(s.begin(), s.end(), back_inserter(face) );
    }
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( point.size(), st.size(), point.size()*2);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        for(int i= 0; i<point.size()/3; ++i) {
            B.add_vertex(point[i*3]);
            B.add_vertex(point[i*3+1]);
            B.add_vertex(point[i*3+2]);
        }
        for(auto itr=st.begin(); itr!=st.end(); ++itr){
            B.begin_facet();
            B.add_vertex_to_facet(*itr);
            itr++;
            if(itr==st.end())break;
            B.add_vertex_to_facet(*itr);
            itr++;
            if(itr==st.end())break;
            B.add_vertex_to_facet(*itr);
            B.end_facet();
        }
        B.end_surface();
        
    }
private:
    std::vector<Point> point;
    std::vector<int> face;
    std::map<Point, int> mp;
    std::set<int> st;
};