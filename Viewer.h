#ifndef VIEWER_H
#define VIEWER_H

#include "Camera.h"
#include "Shader.h"
#include "types.h"
#include "TriMesh.h"
#include <set>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;


class Viewer
{
public:
    // Parameters
   static int  windowSize[2];
   static bool renderWireframe;
   static bool renderSelected;

   // Data
   static std::set<unsigned> selectedVert;
   static TriMesh*        meshPtr;
    static TriMesh marching;
   static int    verbose;
    //static TriMesh points;
   
   // GL Data
   static Camera camera;
   static Shader shader;
   static GLuint surfaceDL;  

   // init
   static void init(int argc, char** argv);
   static void initGLUT(int argc, char** argv);
   static void initGLSL();

   // GLUT callbacks
   static void mouse(int button, int state, int x, int y);
   static void keyboard(unsigned char c, int x, int y);
   static void special(int i, int x, int y);
   static void motion(int x, int y);
   static void display();
   static void idle();

   // GL setup   
   static void setGL();
   static void setLighting();
   static void callDisplayList();
   static void updateDisplayList();

   // draw routines
   static void drawScene();
   static void drawVerts();
   static void drawPolygons();
   static void drawWireframe();
   static void drawSelectedVerts();
   static void pickVertex(int x, int y);
   static void takeScreenshot();

   // Application
   static void clearData();
   static void pickCenter();
    //added
    static void selectedVertDeformation(Vec3& selected_point,Vec3& selected_normal);
    static void read(const char*);
    void initSeg();
};

#endif // VIEWER_H
