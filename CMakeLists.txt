project(pick)

cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Options (to edit if necessary)
option(USE_DEBUG "Debug mode?"   OFF)
option(USE_64BIT_ARCH "64 arch?" OFF)

if(USE_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else(USE_DEBUG)
    set(CMAKE_BUILD_TYPE Release)
endif(USE_DEBUG)

if(USE_64BIT_ARCH)
    set(LINK_FLAGS -m64)
    add_definitions(-m64)
endif(USE_64BIT_ARCH)


find_package(GLUT)
find_package(OpenGL)
# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
     include( ${EIGEN3_USE_FILE} )
  else()
     message(STATUS "NOTICE: some of the executables in this directory need Eigen 3.1 (or greater)")
  endif()
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
 
 # include helper file
endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()
set(LIB_SRC  TriMesh.cpp types.cpp Camera.cpp Image.cpp Shader.cpp Viewer.cpp CIsoSurface.cpp Vectors.cpp )
set(SRCS ${LIB_SRC} main.cpp)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-long-long -DGL_GLEXT_PROTOTYPES")

include_directories(${GLUT_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})


if(APPLE)
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -framework Accelerate")
endif(APPLE)

add_executable(${PROJECT_NAME} ${SRCS})
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
