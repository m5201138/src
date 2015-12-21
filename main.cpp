// is getopt() in stdlib?
#include <cstdlib>
#include <iostream>
#include "Viewer.h"
#include "TriMesh.h"

bool parseCmdLine(int argc, char** argv, 
                  char* & meshFile, 
                  double& scale,
                  int   & verbose)
{
    int arg;
    while ((arg = getopt(argc, argv, "m:v:s:")) != -1)
    {
        switch (arg)
        {
            // Mesh file
            case 'm':
                meshFile = optarg;
                break;
            // Extra parameters
            case 's':
                scale = atof(optarg);
                break;
            case 'v':
                verbose = atoi(optarg);
                break;
        }
    }

    if (meshFile == 0) 
        std::cout << "Usage: " << argv[0] 
        << " -m mesh.obj"
        << " [-s scale=1]"
        << " [-v verbose=0]"
        << std::endl;

    return meshFile;
}

int main(int argc, char** argv)
{
    int    verbose  = 0;
    char*  meshFile = 0;
    double scale    = 1.0;
    if (!parseCmdLine(argc, argv, meshFile, scale, verbose)) return 0;

    TriMesh* mesh=new TriMesh;
    mesh->read(meshFile);
    mesh->normalize();
        
    Viewer viewer;
    viewer.meshPtr = mesh;
    viewer.verbose = verbose;
    viewer.clearData();
    viewer.init(argc,argv);
    return 1;
}
