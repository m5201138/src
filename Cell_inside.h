// Predicate:
// operator() determines if a cell is inside the solid defined by a given implicit function

#ifndef CELL_INSIDE_H
#define CELL_INSIDE_H

#include <CGAL/centroid.h>

template <class Triangulation, class Function>
class Cell_inside {
 public:
  typedef typename Triangulation::Cell_handle Cell_handle;
  typedef typename Triangulation::Tetrahedron Tetrahedron;
  typedef typename Triangulation::Point Point;

 private:
  Function fun;
  Triangulation& tr;
  
 public:
  Cell_inside(Triangulation& t, Function& f)
      : tr(t), fun(f) {}

  bool operator()(Cell_handle& cell_handle) {
    if (tr.is_infinite(cell_handle))
      return false;
    
    Tetrahedron tet = tr.tetrahedron(cell_handle);
    Point centroid = CGAL::centroid(tet);
    double f = fun(centroid);
    return f < -1.0e-5;//0.0;
  }
};

#endif
