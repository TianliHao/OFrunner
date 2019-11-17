#include "OfGenerator/of_processor.h"
#include "OfInterpolator/of_interpolator.h"
#include "QuickSolver/quick_solver.h"

//define global variables
//////////FILES parameters
std::string         g_project_path;
std::string         g_model_folder;
std::string         g_model_name;
std::string         g_OFtemplate_folder;
std::string         g_output_folder;
std::string         g_dataset_folder;
//////////GEOMETRY parameters
Eigen::Vector3d     g_rotation_axis;
double              g_rotation_angle;
double              g_model_length;
double              g_diagonal_length;
double              g_ref_length;
double              g_ref_area;
int                 g_render_resolution;
//////////CFD parameters
double              g_velocity;
double              g_viscosity;
double              g_flow_density;
//////////manual factors
double              g_timestep_manual_factor;
double              g_coeff_avg_start_rate;

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


/*************************************/
int g_program_mode;
std::vector<std::vector<double>> g_dataset;
QuickSolver g_quick_solver;
OFInterpolator g_of_interpolator;
OFProcessor g_of_processor;
double g_render_deltaT;

int main() {
    if(!LoadConfig("Config"))
        {std::cout<< "Cannot Find The Config File. Please put it at the same directory. Exit. " <<std::endl; return 9;}



    g_of_processor.Init();
    if(g_program_mode==2)
    {
        g_of_interpolator.Init();
        g_quick_solver.Init();

        g_of_processor.RunViewer();
    }
    
    if(g_program_mode==3)
    {
        g_of_interpolator.Init();
    }
    
    

    // std::vector<double> test_lowd_input={31.3,0,1,0};
    // std::vector<double> test_output_vec=of_interpolator.MinDistInterpolator(test_lowd_input);

    // std::cout<<"/*/*/*/*/*/*/* interpolation output: "<<std::endl;
    // for(int i=0;i<test_output_vec.size();i++)
    // { 
    //     std::cout<<test_output_vec[i]<<"\t";
    // }std::cout<<std::endl;
    


    return 9999;
    }