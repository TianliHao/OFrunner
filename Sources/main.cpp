#include "OfGenerator/of_processor.h"

//define global variables
//////////FILES parameters
std::string         g_project_path;
std::string         g_model_folder;
std::string         g_model_name;
std::string         g_OFtemplate_folder;
std::string         g_output_folder;
//////////GEOMETRY parameters
Eigen::Vector3d     g_rotation_axis;
double              g_rotation_angle;
double              g_model_length;
double              g_ref_length;
double              g_ref_area;
int                 g_render_resolution;
//////////CFD parameters
double              g_velocity;
double              g_viscosity;

///////////////not input global variables
Eigen::Vector3d     g_avg_normal;
double              g_projection_area;
//These should change when other g_ variables change
Mesh                g_model;
Mesh                g_rotated_model;
std::string         g_case_path;
Eigen::Matrix<double,2,3> g_coefficients;
double              g_Re_number;
double              g_time_step;
int                 g_total_step;
double              g_mesh_min_length;
bool                g_delete_result;

//////////BATCH GENERATION
double g_min_Re;
double g_max_Re;
double g_min_length;
double g_max_length;
double g_min_rotation_y;
double g_max_rotation_y;
double g_min_rotation_z;
double g_max_rotation_z;
int g_Re_count;
int g_length_count;
int g_rotation_count_y;
int g_rotation_count_z;


int main() {
    OFProcessor of_processor;
    return 9999;
    }