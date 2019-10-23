#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H
#include "stdlib.h"
#include"time.h"
#include "math.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <dirent.h>
#include "Eigen/Eigenvalues"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <Viewer/common/shader.hpp>
#include <Viewer/common/texture.hpp>
#include <Viewer/common/controls.hpp>
#include <Viewer/common/objloader.hpp>

struct PolyTraits : public OpenMesh::DefaultTraits
{	
    typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	typedef OpenMesh::Vec4f Color;
};
typedef OpenMesh::PolyMesh_ArrayKernelT<PolyTraits>  Mesh;
typedef Mesh::VertexHandle VertexHandle;
typedef Mesh::FaceHandle FaceHandle;
typedef Mesh::VertexIter VertexIter;
typedef Mesh::FaceIter FaceIter;
typedef Mesh::HalfedgeIter HalfedgeIter;
typedef Mesh::HalfedgeHandle HalfedgeHandle;
typedef Mesh::EdgeHandle EdgeHandle;
typedef Mesh::EdgeIter EdgeIter;
typedef Mesh::Point Point;
typedef Mesh::FaceHalfedgeIter FaceHalfedgeIter;
typedef Mesh::VertexVertexIter VertexVertexIter;
typedef Mesh::VertexEdgeIter VertexEdgeIter;
typedef Mesh::VertexOHalfedgeIter VertexOHalfedgeIter;
typedef Mesh::FaceFaceIter FaceFaceIter;
typedef Mesh::FaceEdgeIter FaceEdgeIter;
typedef Mesh::VertexFaceIter VertexFaceIter;
typedef Mesh::FaceVertexIter FaceVertexIter;

//////////FILES parameters
extern std::string         g_project_path;
extern std::string         g_model_folder;
extern std::string         g_model_name;
extern std::string         g_OFtemplate_folder;
extern std::string         g_output_folder;
//////////GEOMETRY parameters
extern Eigen::Vector3d     g_rotation_axis;
extern double              g_rotation_angle;
extern double              g_model_length;
extern double              g_ref_length;
extern double              g_ref_area;
extern int				   g_render_resolution;
//////////CFD parameters
extern double              g_velocity;
extern double              g_viscosity;
//////////////////////////////////////////////
extern Eigen::Vector3d     g_avg_normal;
extern double              g_projection_area;
extern Mesh                g_model;
extern Mesh                g_rotated_model;
extern std::string         g_case_path;
extern Eigen::Matrix<double,2,3> g_coefficients;
extern double			   g_Re_number;
extern double			   g_time_step;
extern int				   g_total_step;
extern double              g_mesh_min_length;
extern bool				   g_delete_result;

//BATCH GENERATING/////////////////////////////
extern double g_min_Re;
extern double g_max_Re;
extern double g_min_length;
extern double g_max_length;
extern double g_min_rotation_y;
extern double g_max_rotation_y;
extern double g_min_rotation_z;
extern double g_max_rotation_z;
extern int g_Re_count;
extern int g_length_count;
extern int g_rotation_count_y;
extern int g_rotation_count_z;




bool LoadConfig(std::string filename);


#endif