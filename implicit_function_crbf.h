#ifndef IMPLICIT
#define IMPLICIT


#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <eigen/Eigen/LU>
#include <eigen/Eigen/Cholesky>
#include <CGAL/Implicit_surface_3.h>
#include "HRBF_closed.h"
 template<typename FT,typename Point>class implicit_function_crbf{
 public:
typedef Eigen::Matrix<double,3,1> Vector3;
     HRBF_closed crbf;
     implicit_function_crbf(HRBF_closed c){
         crbf=c;
     }
     FT operator()(Point p) const{
     Vector3 pt(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
         return crbf.eval(pt);
                }
 };


#endif