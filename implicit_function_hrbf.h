#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <eigen/Eigen/LU>
#include <eigen/Eigen/Cholesky>
#include <vector>
#include <iostream>
#include "hrbf_core.h"
#include "hrbf_phi_funcs.h"
#include <CGAL/Implicit_surface_3.h>
 template<typename FT,typename Point>class implicit_function_hrbf
{
    public:
    typedef HRBF_fit<double, 3, Rbf_pow3<double> > HRBF_fit;
    HRBF_fit hrbf;
    //typedef  Scalar;
    typedef Eigen::Matrix<double,3,1> Vector;
    //typedef GT Geom_traits;
    //typedef Function_ Function;
    //typedef typename Geom_traits::Sphere_3 Sphere_3;
    //typedef typename Geom_traits::FT FT;
    //typedef typename Geom_traits::Point_3 Point;
    
    //_Scalar HRBF_fit<_Scalar, _Dim, Rbf>::(*fpEval)(Vector&)
    //fpEval= _Scalar HRBF_fit<_Scalar,_Dim,Rbf>::eval
    implicit_function_hrbf(HRBF_fit h ){
        hrbf=h;
        
    }
    FT operator()(Point p) const { return hrbf.eval(Vector(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()))); }
};
