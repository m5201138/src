#ifndef HRBF_CLOSED
#define HRBF_CLOSED


#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <eigen/Eigen/LU>
#include <eigen/Eigen/Cholesky>
#include <CGAL/Implicit_surface_3.h>
class HRBF_closed{
typedef Eigen::Matrix<double,3,1> Vector3;
 private:
     std::vector<Vector3> pts;
     std::vector<Vector3> normals;
     std::vector<Vector3> centers;
    //double d;
     double rho, eta;
 public:
    HRBF_closed() {}
     HRBF_closed(std::vector<Vector3>& pt,std::vector<Vector3>& normal,std::vector<Vector3>& cent,double r, double e){
        copy(pt.begin(), pt.end(), back_inserter(pts));
         copy(normal.begin(), normal.end(), back_inserter(normals) );
         copy(cent.begin(), cent.end(), back_inserter(centers) );
         //copy(d_tmp.begin(),d_tmp.end(),back_inserter(d));
         rho=r;
         eta=e;

         //for(unsigned int k=0;k<pts.size();k++){
           //  d=0;
         //}
         
     }
     

     double eval(Vector3 p) const{
         std::vector<double> g(centers.size()*3);
         for(unsigned int i=0;i<centers.size();i++){
             g[i*3] = 0.0;
             g[i*3+1] = 0.0;
             g[i*3+2] = 0.0;
             double r=sqrt(pow((p(0)-centers[i](0)),2)+pow((p(1)-centers[i](1)),2)+pow((p(2)-centers[i](2)),2));
             
             if (r <= rho&&(p(0)!=centers[i](0)||p(1)!=centers[i](1)||p(2)!=centers[i](2))!=0){
                 g[i*3]=-20.0/pow(rho,2)*(p(0)-centers[i](0))*pow((1.0-r/rho),3);
                 g[i*3+1]=-20.0/pow(rho,2)*(p(1)-centers[i](1))*pow((1.0-r/rho),3);
                 g[i*3+2]=-20.0/pow(rho,2)*(p(2)-centers[i](2))*pow((1.0-r/rho),3);
             }
         }

         double d=0;
         for(unsigned int i=0;i<centers.size();i++){
             double rho2 = rho*rho;
             Vector3 nr;
             
             nr(0)=rho2/(20.0+eta*rho2)*normals[i](0);
             nr(1)=rho2/(20.0+eta*rho2)*normals[i](1);
             nr(2)=rho2/(20.0+eta*rho2)*normals[i](2);
             
             double dt=nr(0)*g[i*3]+nr(1)*g[i*3+1]+nr(2)*g[i*3+2];
             
             d-=dt;

         }
         //std::cout<<p.x()<<" "<<p.y()<<" "<<p.z()<<" "<<d<<std::endl;
         return d;
     }
 };


#endif