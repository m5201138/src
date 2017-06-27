project(test)

cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

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


# CGAL and its components
find_package( CGAL REQUIRED COMPONENTS Core )


if ( CGAL_FOUND )
include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

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
set(SRCS ${LIB_SRC} segmentation.cpp)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-long-long")

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})


if(APPLE)
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -framework Accelerate")
endif(APPLE)

add_executable(${PROJECT_NAME} ${SRCS})
target_link_libraries(${PROJECT_NAME})
