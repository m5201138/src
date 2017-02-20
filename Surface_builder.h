//
// https://svn.zib.de/lenne3d/!svn/bc/1964/lib/qpl_cgal/3.5.1/include/CGAL/IO/Complex_2_in_triangulation_3_polyhedron_builder.h
//

#ifndef SURFACE_BUILDER_H
#define SURFACE_BUILDER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <map>
#include <set>
#include <stack>
#include <utility>


template <class Triangulation, class Cell_inside, class Polyhedron>
class Surface_builder 
  : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> 
{
  // Types
public:
  typedef typename Polyhedron::HalfedgeDS HDS;

private:
  typedef CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> Base;
  typedef typename Triangulation::Finite_facets_iterator Finite_facets_iterator;

  typedef typename Triangulation::Cell_handle Cell_handle;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Point Point;
  typedef typename Triangulation::Geom_traits::Vector_3 Vector;
  typedef typename Triangulation::Edge Edge;
  typedef typename Triangulation::Facet Facet;
  typedef typename Triangulation::Finite_vertices_iterator Finite_vertices_iterator;

    
  // Fields
private:
  Triangulation& tr;
  Cell_inside cell_inside;

  
  // Methods
public:
  Surface_builder(Triangulation& t, Cell_inside& cell_in)
    : Base(), tr(t), cell_inside(cell_in) {}

  void operator() (HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
    int number_of_faces = compute_number_of_faces();
    builder.begin_surface(tr.number_of_vertices(), number_of_faces);
    {
      // coordinates of vertices
      std::map<Vertex_handle, int> V;
      int inum = 0;
      for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
	   vit != tr.finite_vertices_end();
	   ++vit)
      {
        V[vit] = inum++;
        Point p = static_cast<Point>(vit->point());
        builder.add_vertex(p);
      }
      
      // oriented faces
      
      Finite_facets_iterator fit = tr.finite_facets_begin();
      int number_of_faces = compute_number_of_faces();
      std::set<Facet> oriented_set;
      while (oriented_set.size() != number_of_faces) {
	while (!select_facet(*fit) || 
	       oriented_set.find(*fit) != oriented_set.end()) 
        {
          ++fit;
        }
        
	oriented_set.insert(*fit);
      }

      typename std::set<Facet>::const_iterator sfit;
      for (sfit = oriented_set.begin();
	   sfit != oriented_set.end();
	   ++sfit)
      {

        Cell_handle cell_handle = sfit->first;
        int vert_index = sfit->second;

        int vert_indices[3];
        if (cell_inside(cell_handle)) {
          for (int i = 0; i < 3; ++i) {
            vert_indices[i] = cw_order(vert_index, i);
          }
        } else {
          for (int i = 0; i < 3; ++i) {
            vert_indices[i] = ccw_order(vert_index, i);
          }
        }
        
        int indices[3];
        int index = 0;
        for (int i = 0; i < 3; ++i) {
          /*
            indices[index++] = 
	      V[sfit->first->vertex(tr.vertex_triple_index(sfit->second, i))];
          */
          indices[index++] = V[sfit->first->vertex(vert_indices[i])];
        }

        if (builder.test_facet(indices+0, indices+3))
          builder.add_facet(indices+0, indices+3);
        else {
          /*
          builder.begin_facet();
          builder.add_vertex_to_facet(indices[2]);
          builder.add_vertex_to_facet(indices[1]);
          builder.add_vertex_to_facet(indices[0]);
          builder.end_facet();
          */
        }
      }
    }
    builder.end_surface();
    
  }
  
private:
  // compute number of faces in tr that satisfies the predicate fs()
  int compute_number_of_faces() {
    int num_faces = 0;
    Finite_facets_iterator fit;
    for (fit = tr.finite_facets_begin();
	 fit != tr.finite_facets_end();
	 ++fit) 

      {
	if (select_facet(*fit)) {
	  num_faces++;
	}
      }

    return num_faces;
  }

  // decide if a facet belongs to the surface based on
  // its 2 adjacent cells
  template<typename Facet>
  bool select_facet(Facet& facet) {
    Cell_handle cell_handle = facet.first;
    int vertex_index = facet.second;

    Cell_handle opposite_cell_handle = cell_handle->neighbor(vertex_index);

    bool is_cell_in = cell_inside(cell_handle);
    bool is_opposite_cell_in = cell_inside(opposite_cell_handle);

    if (is_cell_in && is_opposite_cell_in) return false;
    if (!is_cell_in && !is_opposite_cell_in) return false;
    return true;
  }

  int ccw_order(int vert, int position) {
    int ccw[4][3] = {
      {1,3,2},
      {2,3,0},
      {0,3,1},
      {0,1,2}
    };

    return ccw[vert][position];
  }

  int cw_order(int vert, int position) {
    int cw[4][3] = {
      {1,2,3},
      {2,0,3},
      {0,1,3},
      {0,2,1}
    };

    return cw[vert][position];
  }
};

#endif
