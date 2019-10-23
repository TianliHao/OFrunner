#include "Input/global_variables.h"
using namespace std;
bool LoadConfig(string filename)
{
    ifstream fin(filename);
	if (fin.fail() == true)
	{
		cout << "Cannot load config file." << endl;
		return false;
	}
	string param;
	while (fin >> param)
	{
        if (param == "g_project_path")
            fin >> g_project_path;
		else if (param == "g_model_folder")
            fin >> g_model_folder;
        else if (param == "g_model_name")
            fin >> g_model_name;
        else if (param == "g_OFtemplate_folder")
            fin >> g_OFtemplate_folder;
        else if (param == "g_output_folder")
            fin >> g_output_folder;
        ////    
        else if (param == "g_rotation_axis")
            fin >> g_rotation_axis[0] >> g_rotation_axis[1] >> g_rotation_axis[2];
        else if (param == "g_rotation_angle")
            fin >> g_rotation_angle;
        else if (param == "g_model_length")
            fin >> g_model_length;
        else if (param == "g_ref_length")
            fin >> g_ref_length;
        else if (param == "g_ref_area")
            fin >> g_ref_area;
        else if (param == "g_render_resolution")
            fin >> g_render_resolution;
        ////
        else if (param == "g_velocity")
            fin >> g_velocity;
        else if (param == "g_Re_number")
            fin >> g_Re_number;
        else if (param == "g_viscosity")
            fin >> g_viscosity;
        else if (param == "g_time_step")
            fin >> g_time_step;
        else if (param == "g_total_step")
            fin >> g_total_step;
        else if (param == "g_delete_result")
        {
            int flag;
            fin>>flag;
            if(flag!=0) g_delete_result=true;
            else g_delete_result=false;
        }
//////////////////////
        else if (param == "g_min_Re")
            fin >> g_min_Re;
        else if (param == "g_max_Re")
            fin >> g_max_Re;
        else if (param == "g_min_length")
            fin >> g_min_length;
        else if (param == "g_max_length")
            fin >> g_max_length;
        else if (param == "g_min_rotation_y")
            fin >> g_min_rotation_y;
        else if (param == "g_max_rotation_y")
            fin >> g_max_rotation_y;
        else if (param == "g_min_rotation_z")
            fin >> g_min_rotation_z;
        else if (param == "g_max_rotation_z")
            fin >> g_max_rotation_z;
        else if (param == "g_Re_count")
            fin >> g_Re_count;
        else if (param == "g_length_count")
            fin >> g_length_count;
        else if (param == "g_rotation_count_y")
            fin >> g_rotation_count_y;
        else if (param == "g_rotation_count_z")
            fin >> g_rotation_count_z;
	}
	cout << "Config File Loaded Successfully." << endl;
	return true;
}

