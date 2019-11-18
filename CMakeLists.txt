cmake_minimum_required(VERSION 3.0.0)
project(OFRunner VERSION 0.1.0)

include(CTest)
enable_testing()

FIND_PACKAGE(OpenGL REQUIRED)
IF(OPENGL_FOUND)
	MESSAGE("===== support OpenGL =====")
	MESSAGE(STATUS "OpenGL library status:")
	MESSAGE(STATUS "	include path: ${OPENGL_INCLUDE_DIR}")
	MESSAGE(STATUS "	libraries: ${OPENGL_LIBRARIES}")
	INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
ELSE()
	MESSAGE("##### not support OpenGL #####")
ENDIF()
FIND_PACKAGE(GLUT REQUIRED)
IF(GLUT_FOUND)
	MESSAGE("===== support GLUT =====")
	MESSAGE(STATUS "	include path: ${GLUT_INCLUDE_DIR}")
	MESSAGE(STATUS "	libraries: ${GLUT_LIBRARIES}")
	INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
ELSE()
	MESSAGE("##### not support GLUT #####")
ENDIF()
find_package( X11 REQUIRED )
  list( APPEND GLSLCOOKBOOK_LIBS ${X11_Xrandr_LIB} ${X11_Xi_LIB} ${X11_Xrender_LIB} ${X11_Xdamage_LIB} ${X11_Xxf86vm_LIB} ${X11_Xext_LIB}  ${X11_X11_LIB} ${X11_Xinerama_LIB}  ${X11_Xcursor_LIB})


SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O2 -std=c11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -g -Wall -O2 -std=c++11 -lX11 -ldl")

link_directories(/home/tlh/Dev/OpenMesh-8.0/BUILD/Build/lib)
link_libraries("libOpenMeshCore.a" "libOpenMeshTools.a")

link_libraries(${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${X11_LIBRARIES} pthread "libglfw3.a" "libGLEW.so")


include_directories(./Sources)
include_directories(/usr/local/include/eigen3)
include_directories(/home/tlh/Dev/OpenMesh-8.0/src)


aux_source_directory(./Sources SourceFiles)
aux_source_directory(./Sources/Input SourceFiles)
aux_source_directory(./Sources/OfGenerator SourceFiles)
aux_source_directory(./Sources/OfInterpolator SourceFiles)
aux_source_directory(./Sources/QuickSolver SourceFiles)
aux_source_directory(./Sources/Viewer SourceFiles)
aux_source_directory(./Sources/Viewer/common SourceFiles)
add_executable(OFRunner ${SourceFiles})

set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_SOURCE_DIR}/Bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
