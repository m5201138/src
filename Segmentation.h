#include "TriMesh.h"
#include <set>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <eigen/Eigen/Core>

typedef Eigen::Matrix<double,3,1> Vector3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
class Segmentation{
public:
void init();
TriMesh*   meshPtr;
private:
    
std::vector<std::vector<Vector3> > pointsSets;
std::vector<std::vector<Vector3> > normalsSets;
std::vector<Poisson_reconstruction_function> functionSets;
std::vector<Vector3> maxPoints;
std::vector<Vector3> maxNormals;
std::vector<Poisson_reconstruction_function> originalFunction;
};
