#include "OfGenerator/of_processor.h"
using namespace std;

void OFProcessor::Init()
{
    OpenMesh::IO::_OBJReader_();
    OpenMesh::IO::_OBJWriter_();
    ChangeGlobalVariables();
    if(g_program_mode!=1)
        return;

    RunViewer();
    ComputeCasesForOneModel(); 
}

void OFProcessor::RunViewer()
{
    Viewer viewer;
}

void ComputeVNormal(Mesh &mesh)
{
    mesh.release_face_normals();
    mesh.release_vertex_normals();
    mesh.request_face_normals();
    mesh.update_face_normals();
    mesh.request_vertex_normals();
    mesh.update_vertex_normals();
};
void OFProcessor::LoadObject()
{
    OpenMesh::IO::read_mesh(g_model, g_project_path+"/"+g_model_folder+"/"+g_model_name);
    NormalizeModel();
    g_model_length=g_model_length/g_diagonal_length;
    NormalizeModel();
    RotateModel();
    g_rotated_model=rotated_model_;
    //rotated_model_=g_model;
    ComputeVNormal(g_rotated_model);
    if(g_program_mode==1)
        RunViewer();
}

void OFProcessor::GenerateOneCase()
{
    successful_run_=false;
    int total_step=g_total_step;
    double timestep_manual_factor=g_timestep_manual_factor;
    ChangeGlobalVariables();
    for(int i=0;i<1;i++)
    {
        double rerun_rate=0.7;
        GenerateOneInput();
        AdjustTemplate();
        RunOneCase();
        ReadOneOFOutput();
        if(successful_run_) break;
        g_timestep_manual_factor=g_timestep_manual_factor*rerun_rate;
        g_total_step=(int)(g_total_step/rerun_rate);
    }
    g_total_step= total_step;
    g_timestep_manual_factor=timestep_manual_factor;
}

void OFProcessor::GenerateOneInput()
{//This function is for generating one OpenFoam input case
    string command = "cd " + g_project_path + ";";
    command += "mkdir " + g_output_folder+"/"+g_model_name + ";";
    command += "cp -r " + g_OFtemplate_folder + " " + g_output_folder+"/"+g_model_name + ";";
    command += "rm -r " + g_case_path + ";";
    command += "mv " + g_output_folder+"/"+g_model_name+"/"+g_OFtemplate_folder + " "
        + g_case_path + ";";
    std::cout<<command<<std::endl;
    system(command.c_str());

    OpenMesh::IO::write_mesh(rotated_model_, g_project_path+"/"+g_case_path+"/constant/triSurface/model.obj");
}

void OFProcessor::NormalizeModel()
{//This function is for normalizing the model coordinates to [-0.5,0.5]
    double normalized_length=g_model_length;
    double max_dist=0;
    double min_coord[3]={DBL_MAX}; double max_coord[3]={DBL_MIN};
    for(int i=0;i<g_model.n_vertices();i++)
    {
        OpenMesh::Vec3d v_coord = g_model.point(g_model.vertex_handle(i));
        for(int a=0;a<3;a++)
        {
            if(v_coord[a]<min_coord[a]) min_coord[a]=v_coord[a];
            if(v_coord[a]>max_coord[a]) max_coord[a]=v_coord[a];
            if(max_coord[a]-min_coord[a]>max_dist) max_dist=max_coord[a]-min_coord[a];
        }
    }
    for(int i=0;i<g_model.n_vertices();i++)
    {
        OpenMesh::Vec3d v_coord = g_model.point(g_model.vertex_handle(i));
        for(int a=0;a<3;a++)
        {
            v_coord[a]-=(max_coord[a]+min_coord[a])/2;
            v_coord[a]*=normalized_length/max_dist;
        }
        g_model.set_point(g_model.vertex_handle(i),v_coord);
    }

    g_diagonal_length=sqrt(pow(max_coord[0]-min_coord[0],2)+pow(max_coord[1]-min_coord[1],2)+pow(max_coord[2]-min_coord[2],2))
                        *normalized_length/max_dist;
    cout<<"g_diagnoal_length: "<<g_diagonal_length<<endl;
}

void OFProcessor::RotateModel()
{
    rotated_model_.clean();
    rotated_model_=g_model;
    Eigen::AngleAxisd t_V(g_rotation_angle/180.0*M_PI, g_rotation_axis);
    for(int i=0;i<rotated_model_.n_vertices();i++)
    {
        OpenMesh::Vec3d v_coord = rotated_model_.point(rotated_model_.vertex_handle(i));
        Eigen::Vector3d point(v_coord[0],v_coord[1],v_coord[2]);
        point=t_V*point;
        v_coord=OpenMesh::Vec3d(point[0],point[1],point[2]);
        rotated_model_.set_point(rotated_model_.vertex_handle(i),v_coord);
    }
}

void OFProcessor::RunOneCase()
{
    double actual_total_step=g_total_step;
    double Co_test_total_step=30;

    //get mesh
    string command = "cd " + g_project_path + "/" + g_case_path + "; ./Allclean; ./GenerateMesh";
    system(command.c_str());
    ReadMinEdgeLength();

    for(int i=0;i<1;i++)
    {
        //run short time CFD to get Co number, then adjust time step
        g_total_step=Co_test_total_step;
        AdjustTemplatecontrolDict();
        command = "cd " + g_project_path + "/" + g_case_path + "; ./RunCFD";
        system(command.c_str());
        ReadCourantNumber();
        command = "cd " + g_project_path + "/" + g_case_path + "; rm -r 0.* [1-9]* post* proce*/0.* proce*/[1-9]* log.patchSummary log.potentialFoam; mv log.pisoFoam log.pisoFoam"+to_string(i);
        system(command.c_str());
    }

    //run CFD
    g_total_step=actual_total_step;
    AdjustTemplatecontrolDict();
    command = "cd " + g_project_path + "/" + g_case_path + "; ./RunCFD";
    system(command.c_str());

    //save disk space
    if(g_delete_result==true)
    {
        command = "cd " + g_project_path + "/" + g_case_path + "; rm -r 0.* [1-9]* proce*";
        system(command.c_str());
    }

}

void OFProcessor::ReadMinEdgeLength()
{
    string log_name=g_project_path + "/" + g_case_path +"/log.checkMesh";
    ifstream login(log_name);
    if(login.fail())
    {
        cout<<"Cannot load checkMesh logfile."<<endl;
        return;
    }
    string word;
    while(login>>word)
    {
        if(word!="Min/max") continue;
        login>>word;
        if(word!="edge") continue;
        login>>word>>word>>g_mesh_min_length;
        break;
    }
    login.close();
    cout<<"3d model minimal length is: "<<g_mesh_min_length<<endl;

    g_time_step=g_mesh_min_length/g_velocity*0.001;
    cout<<"time step is: "<<g_time_step<<endl;

}

void OFProcessor::ReadCourantNumber()
{
    int start_count=15;

    string log_name=g_project_path + "/" + g_case_path +"/log.pisoFoam";
    ifstream login(log_name);
    if(login.fail())
    {
        cout<<"Cannot load pisoFoam logfile."<<endl;
        return;
    }
    string word;
    double mean_Co=0, max_Co=0;
    int count=0;
    int devide_count=0;
    while(login>>word)
    {
        double num1, num2;
        if(word!="Courant") continue;
        login>>word;
        if(word!="Number") continue;
        login>>word>>num1>>word>>num2;
        if(count>start_count)
        {
            mean_Co+=num1; max_Co+=num2;
            devide_count++;
        }
        count++;
    }
    login.close();
    mean_Co/=(double)devide_count; max_Co/=(double)devide_count;
    cout<<"Count"<<devide_count<<" Courant Number mean: "<<mean_Co<<"||||| max: "<<max_Co<<endl;

    //Change the time step to make sure Co number is about 1
    double rate=1.0/max_Co;
    g_time_step*=rate*g_timestep_manual_factor;
    cout<<"|||||||||||||New time step is: "<<g_time_step<<endl;
}


void OFProcessor::ReadOneOFOutput()
{//Start index of this file is 1
    string log_name=g_project_path+"/"+g_case_path+"/../CFDlog.TLH";//TLH means Transient Line Holding data
    string dataset_name=g_project_path+"/Dataset/"+g_model_name+".TLH";

    int average_start_index=g_total_step*g_coeff_avg_start_rate;//values in former interations are not stable enough
    ifstream fin[3];
    fin[0]=ifstream(g_project_path + "/" + g_case_path +"/postProcessing/forceCoeffs1/0/forceCoeffs.dat");
    fin[1]=ifstream(g_project_path + "/" + g_case_path +"/postProcessing/forceCoeffs2/0/forceCoeffs.dat");
    fin[2]=ifstream(g_project_path + "/" + g_case_path +"/postProcessing/forceCoeffs3/0/forceCoeffs.dat");
    ifstream forcein(g_project_path + "/" + g_case_path +"/postProcessing/forces/0/forces.dat");
    if (forcein.fail() == true)
    {
        cout << "Cannot load forces file." << endl;
        return;
    }
    string file_head;
    int output_index=0;
    string output_str[18];
    double output_value[18];
    Eigen::Vector3d c_force_average(0,0,0);
    Eigen::Vector3d c_torque_average(0,0,0);
    int file_head_line_count=0;
    while(getline(forcein,file_head)&&file_head_line_count<3) file_head_line_count++;
    while(forcein>>file_head)
    {
        output_index++;
        
        char tmpchar;
        forcein>>tmpchar;
        for(int i=0;i<18;i++)
        {
            if(i%3==0) forcein>>tmpchar>>output_value[i];
            if(i%3==1) forcein>>output_value[i];
            if(i%3==2) forcein>>output_value[i]>>tmpchar;
            if(i==8) forcein>>tmpchar>>tmpchar;
        }
        forcein>>tmpchar;
        if(output_index>=average_start_index)
        {
            c_force_average[0]+=output_value[0]+output_value[3];
            c_force_average[1]+=output_value[1]+output_value[4];
            c_force_average[2]+=output_value[2]+output_value[5];
            c_torque_average[0]+=output_value[9]+output_value[12];
            c_torque_average[1]+=output_value[10]+output_value[13];
            c_torque_average[2]+=output_value[11]+output_value[14];
        }
    }
    if(output_index>=average_start_index)
    {
        c_force_average/=(double)(output_index-average_start_index+1);
        c_torque_average/=(double)(output_index-average_start_index+1);
        cout<<"force: "<<c_force_average<<endl;
        cout<<"moment: "<<c_torque_average<<endl;
    }
    
    for(int a=0;a<3;a++)
    {
    if (fin[a].fail() == true)
        {
            cout << "Cannot load forceCoeffs file." << endl;
            return;
        }
        string file_head;
        double output_value[5];
        output_index=0;
        double c_force_value_average=0;
        double c_torque_value_average=0;
        int file_head_line_count=0;
        while(getline(fin[a],file_head)&&file_head_line_count<9) file_head_line_count++;
        while(fin[a]>>file_head)
        {
            output_index++;
            fin[a]>>output_value[0]>>output_value[1]>>output_value[2]>>output_value[3]>>output_value[4];
            if(output_index>=average_start_index)
            {

                c_force_value_average+=output_value[1];
                c_torque_value_average+=output_value[0];
            }
        }
        if(output_index>=average_start_index)
        {
            c_force_value_average/=(double)(output_index-average_start_index+1);
            c_torque_value_average/=(double)(output_index-average_start_index+1);
            g_coefficients(0,a)=c_force_value_average;
            g_coefficients(1,a)=c_torque_value_average;
            cout<<"coord: "<<a<<endl;
            cout<<"force coefficient: "<<g_coefficients(0,a)<<endl;
            cout<<"moment coefficient: "<<g_coefficients(1,a)<<endl;
        }
    }
    //do not output datapair for one-case-test
    if(g_one_case_CFD_test!=0)
        return;

    if(output_index>=g_total_step-1)
    {
        successful_run_=true;
        ofstream of(dataset_name,ios::app);
        of<<g_rotation_angle<<" "<<g_rotation_axis(0)<<" "<<g_rotation_axis(1)<<" "<<g_rotation_axis(2)<<"    "
            <<g_diagonal_length<<" "<<g_projection_area<<" "<<g_velocity<<" "
            <<g_Re_number<<"    "<<g_coefficients(0,0)<<" "
            <<g_coefficients(0,1)<<" "<<g_coefficients(0,2)<<"    "
            <<g_coefficients(1,0)<<" "<<g_coefficients(1,1)<<" "
            <<g_coefficients(1,2)<<"    "
            <<c_force_average(0)<<" "<<c_force_average(1)<<" "<<c_force_average(2)<<"    "
            <<c_torque_average(0)<<" "<<c_torque_average(1)<<" "<<c_torque_average(2)<<"    "
            <<19960420<<" "<<g_alpha_degree<<" "<<g_beta_degree<<"\n";//temp change: save alpha beta, replacing the normal
            //<<g_avg_normal(0)<<" "<<g_avg_normal(1)<<" "<<g_avg_normal(2)<<"\n";
        of.close();
    }
}

void OFProcessor::ChangeGlobalVariables()
{

    LoadObject();
    
    g_case_path = g_output_folder+"/"+g_model_name+"/"
        +to_string(g_rotation_axis[0]).substr(0,to_string(g_rotation_axis[0]).size()-4)+"_"
        +to_string(g_rotation_axis[1]).substr(0,to_string(g_rotation_axis[1]).size()-4)+"_"
        +to_string(g_rotation_axis[2]).substr(0,to_string(g_rotation_axis[2]).size()-4)+"_"
        +to_string(g_rotation_angle).substr(0,to_string(g_rotation_angle).size()-4)+"_"
        +to_string((int)g_Re_number)+"_"
        +to_string(g_model_length).substr(0,to_string(g_model_length).size()-4);

    g_case_path = g_output_folder+"/"+g_model_name+"/"
        +to_string(g_alpha_degree).substr(0,to_string(g_alpha_degree).size()-4)+"_"
        +to_string(g_beta_degree).substr(0,to_string(g_beta_degree).size()-4)+"_"
        +to_string((int)g_Re_number)+"_"
        +to_string(g_model_length).substr(0,to_string(g_model_length).size()-4);
    
    //add flag of testing one case
    if(g_one_case_CFD_test!=0)
        g_case_path+="_OneCaseTest";


    g_velocity=g_Re_number/1.225/g_diagonal_length*1.5e-5;
    cout<<"\n\n********Re Number= "<<g_Re_number<<"\n"<<endl;
    cout<<"*******Velocity= "<<g_velocity<<"\n\n"<<endl;

}

void OFProcessor::ComputeAllCFD()
{//outer loop is different models, inner loops are changing different rotation
    string input_models_root=g_project_path+"/"+g_model_folder+"/";
    struct dirent *ptr;    
    DIR *dir;
    dir=opendir(input_models_root.c_str());
    printf("all models:\n");
    while((ptr=readdir(dir))!=NULL)
    {
        if(ptr->d_name[0] == '.')
            continue;
        string model_name=ptr->d_name;
        cout<<model_name<<"|||||||||||||||||||||||||||||||||||"<<endl;

        //for one model test
        if(model_name!="cube.obj") continue;
        g_model_name=model_name;
        //change rotaion axis
        for(int i=0;i<1;i++)
        {
            //change rotation angle
            double rotation=0, init_rotation=0, max_rotation=180;
            int rotation_num=10;
            for(int a=0;a<rotation_num;a++)
            {
                g_rotation_angle=rotation;
                GenerateOneCase();
                rotation+=(max_rotation-init_rotation)/(double)(rotation_num-1);
            }
        }
    }
    closedir(dir);
}

void OFProcessor::ComputeCasesForOneModel()
{
    double min_Re=g_min_Re;
    double max_Re=g_max_Re;
    double min_length=g_min_length;
    double max_length=g_max_length;
    double min_rotation_y=g_min_rotation_y;
    double max_rotation_y=g_max_rotation_y;
    double min_rotation_z=g_min_rotation_z;
    double max_rotation_z=g_max_rotation_z;
    int Re_count=g_Re_count;
    int length_count=g_length_count;
    int rotation_count_y=g_rotation_count_y;
    int rotation_count_z=g_rotation_count_z;

    if(g_one_case_CFD_test!=0)
    {
        //if g_one_case_CFD_test!=0
        //all the changes:
        //here, ignore sampling number and delete_result_flag
        //ChangeGlobalVariables(): g_case_path, in case of deleting dataset generating files
        //ReadOneOFOutput(): do not generate datapair
        g_rotation_count_y=1;
        g_rotation_count_z=1;
        g_Re_count=1;
        g_length_count=1;
        g_delete_result=0;
    }
    
    //rotation along y-axis
    for(int y=0;y<rotation_count_y;y++)
    {
        Eigen::Vector3d rotation_axis_y(0,1,0);
        double rotation_angle_y=min_rotation_y;
        if(rotation_count_y>1) rotation_angle_y+=(max_rotation_y-min_rotation_y)/(double)(rotation_count_y-1)*y;
        g_alpha_degree=rotation_angle_y;
        rotation_angle_y*=M_PI/180;
        Eigen::AngleAxisd rot_y(rotation_angle_y, rotation_axis_y);

        //rotation along z-axis
        for(int z=0;z<rotation_count_z;z++)
        {
            Eigen::Vector3d rotation_axis_z(0,0,1);
            double rotation_angle_z=min_rotation_z;
            if(rotation_count_z>1) rotation_angle_z+=(max_rotation_z-min_rotation_z)/(double)(rotation_count_z-1)*z;
            g_beta_degree=rotation_angle_z;
            rotation_angle_z*=M_PI/180;
            Eigen::AngleAxisd rot_z(rotation_angle_z, rotation_axis_z);
            //ORDER CHANGED!!
            Eigen::AngleAxisd rot_combine(rot_z.matrix()*rot_y.matrix());//first_mat*second_mat means use local coord
            g_rotation_angle=rot_combine.angle()*180/M_PI;
            g_rotation_axis=rot_combine.axis();
            for(int re=0;re<Re_count;re++)
            {
                double Re_number=min_Re;
                if(Re_count>1) Re_number+=(max_Re-min_Re)/(double)(Re_count-1)*re;
                g_Re_number=Re_number;
                for(int l=0;l<length_count;l++)
                {
                    //g_delete_result=true;////////////////////delete result for saving diskspace

                    double length=min_length;
                    if(length_count>1) length+=(max_length-min_length)/(double)(length_count-1)*l;
                    g_model_length=length;
                
                    GenerateOneCase();
                }
            }            
        }
    }
}

void OFProcessor::GenerateNormalList()
{
    int longitudinal_div=36;
    int latitudinal_div=2;
    double delta_alpha_degree=360/(double)longitudinal_div;
    double delta_beta_degree=180/(double)latitudinal_div;

    std::vector<std::vector<double>> normal_list;
    string listname="../Render/NormalList/"+g_model_name+".NL";
    std::fstream fout(listname.c_str(),std::ios::out);


    int small_longitudinal_div=1;
    for(int lon=0;lon<small_longitudinal_div;lon++)
    {
        for(int lat=0;lat<latitudinal_div+1;lat++)
        {
            g_alpha_degree=delta_alpha_degree*lon;
            g_beta_degree=delta_beta_degree*lat-90;
            Eigen::Vector3d rotation_axis_y(0,1,0);
            Eigen::Vector3d rotation_axis_z(0,0,1);
            Eigen::AngleAxisd rot_y(g_alpha_degree*M_PI/180.0, rotation_axis_y);
            Eigen::AngleAxisd rot_z(g_beta_degree*M_PI/180.0, rotation_axis_z);
            Eigen::AngleAxisd rot_combine(rot_z.matrix()*rot_y.matrix());
            g_rotation_angle=rot_combine.angle()*180/M_PI;
            g_rotation_axis=rot_combine.axis();
            
            RotateModel();
            g_rotated_model=rotated_model_;
            ComputeVNormal(g_rotated_model);
            Viewer viewer;

            std::vector<double> normal_data;
            normal_data.push_back(g_alpha_degree);
            normal_data.push_back(g_beta_degree);
            normal_data.push_back(g_projection_area);
            normal_data.push_back(g_avg_normal[0]);
            normal_data.push_back(g_avg_normal[1]);
            normal_data.push_back(g_avg_normal[2]);
            normal_list.push_back(normal_data);
        }
    }
    std::vector<std::vector<double>> new_normal_list;
    for(int i=0;i<normal_list.size();i++)
    {
        if(g_model_name=="foamdisk.obj"||g_model_name=="foamsemisphere.obj")
        {
            for(int lon=0;lon<longitudinal_div;lon++)
            {
                std::vector<double> new_data;
                g_alpha_degree=delta_alpha_degree*lon;
                new_data=normal_list[i];
                new_data[0]=g_alpha_degree;
                new_normal_list.push_back(new_data);
            }
        }
        if(g_model_name=="realbox.obj")
        {
            for(int a=0;a<4;a++)
            {
                std::vector<double> new_data;
                new_data=normal_list[i];
                new_data[0]=new_data[0]+a*90;
                new_normal_list.push_back(new_data);
            }
        }
    }
    normal_list=new_normal_list;
    g_normal_list=normal_list;

    //write file
    for(int i=0;i<normal_list.size();i++)
    {
        for(int a=0;a<6;a++)
        {
            fout<<normal_list[i][a]<<" ";
        }fout<<endl;
    }
}