#include "OfInterpolator/of_interpolator.h"
#include "OfGenerator/of_processor.h"
using namespace std;
using namespace Eigen;

//data format-index
//rotation_angle, axis, length, area, velo, Re, Cd, Cl, Cs, coeffM, force, torque, avg_normal
//0               123   4       5     6     7   8   9   10  111213  141516 171819  202122
class OFProcessor;

void AlphaBetaRadiansFromRotMat(Matrix3d rot, double &alpha_radians, double &beta_radians)
{
    Vector3d sphere_point(1,0,0);
    sphere_point=rot*sphere_point;
    alpha_radians=-asin(sphere_point[2]/sqrt(sphere_point[0]*sphere_point[0]+sphere_point[2]*sphere_point[2]));
    if(sphere_point[0]<0) alpha_radians=M_PI-alpha_radians;
    beta_radians=asin(sphere_point[1]/1.0);
}

void RotMatFromAlphaBeta(Matrix3d &rot, double alpha_radians, double beta_radians)
{
    Eigen::AngleAxisd rot_alpha(alpha_radians,Vector3d(0,1,0));
    Eigen::AngleAxisd rot_beta(beta_radians,Vector3d(0,0,1));
    rot = rot_alpha.matrix()*rot_beta.matrix();         
}

void ViewXYZFromAlphaBeta(double &x, double &y, double &z, double alpha_radians, double beta_radians)
{
    //get viewpoint xyz from alpha_point
    //note: minus vector go to the upstream of airflow
    Matrix3d alpha_beta_rotation;
    RotMatFromAlphaBeta(alpha_beta_rotation,alpha_radians,beta_radians);
    Vector3d translate_vector(1,0,0);
    translate_vector=alpha_beta_rotation*translate_vector;
    translate_vector=-translate_vector;
    x=translate_vector[0]; y=translate_vector[1]; z=translate_vector[2];
}

double SmallSphereAreaPartFromBeta(double beta_radians, double delta_alpha_radians, double delta_beta_radians)
{
    //use nearby alpha-beta get a trapezoid region
    double lower_line_length=abs(2*sin(0.5*delta_alpha_radians)*cos(beta_radians));
    double smaller_beta_radians=beta_radians;
    if(beta_radians>0) smaller_beta_radians=beta_radians-delta_beta_radians;
    if(beta_radians<0) smaller_beta_radians=beta_radians+delta_beta_radians;
    double upper_line_length=abs(2*sin(0.5*delta_alpha_radians)*cos(smaller_beta_radians));
    double height=abs(2*sin(0.5*delta_beta_radians));
    double small_area=(lower_line_length+upper_line_length)*0.5*height;
    return small_area;
}

void OFInterpolator::GetDatasetForOneModel()
{
    dataset_.clear();
    input_size_=4;//viewpoint_alpha-from_x, viewpoint_beta_from_xz, length, velocity, force, torque
    output_size_=6;

    //keep diagonal length = 1
    //save x, y, z on the sphere
    //save Cd*S, Cm*S*diagonal_l
    //*** the FORCE is inversed, while keeping the torque 

    max_force_=0;max_torque_=0;
    for(int i=0;i<g_dataset.size();i++)
    {
        vector<double> datapair;
        double rotation_angle, diagonal_length, velocity;
        Vector3d rotation_axis, force_in_OF_coord, torque_in_OF_coord;

        rotation_angle=g_dataset[i][0];
        rotation_axis=Vector3d(g_dataset[i][1],g_dataset[i][2],g_dataset[i][3]);
        diagonal_length=g_dataset[i][4];
        velocity=g_dataset[i][6];

        bool output_real_force_torque=true;
        if(output_real_force_torque)
        {
            force_in_OF_coord=Vector3d(g_dataset[i][14],g_dataset[i][15],g_dataset[i][16]);
            torque_in_OF_coord=Vector3d(g_dataset[i][17],g_dataset[i][18],g_dataset[i][19]);
            //use coeff rather than force/torque
            force_in_OF_coord=Vector3d(g_dataset[i][8],g_dataset[i][9],g_dataset[i][10]);
            torque_in_OF_coord=Vector3d(g_dataset[i][11],g_dataset[i][12],g_dataset[i][13]);
            //get rid of area problem
            force_in_OF_coord*=g_dataset[i][5];
            force_in_OF_coord=force_in_OF_coord/g_dataset[i][4]/g_dataset[i][4];
            torque_in_OF_coord*=g_dataset[i][5]*g_dataset[i][4];
            torque_in_OF_coord=torque_in_OF_coord/pow(g_dataset[i][4],3);

            force_in_OF_coord=-force_in_OF_coord;
            torque_in_OF_coord*=-10;
        }
        else
        {
            force_in_OF_coord=Vector3d(g_dataset[i][8],g_dataset[i][9],g_dataset[i][10]);
            torque_in_OF_coord=Vector3d(g_dataset[i][11],g_dataset[i][12],g_dataset[i][13]);
            force_in_OF_coord*=g_dataset[i][5];
            torque_in_OF_coord*=g_dataset[i][5]*g_dataset[i][4];
            //caution! scale by the diagonal normalization!
            force_in_OF_coord/=pow(diagonal_length,2);
            torque_in_OF_coord/=pow(diagonal_length,3);

            // //make force&torque become the Coeff Cd*A and Cm*A*l: 
            force_in_OF_coord=Vector3d(g_dataset[i][14],g_dataset[i][15],g_dataset[i][16]);
            torque_in_OF_coord=Vector3d(g_dataset[i][17],g_dataset[i][18],g_dataset[i][19]);
            force_in_OF_coord/=0.5*velocity*velocity*g_flow_density;
            torque_in_OF_coord/=0.5*velocity*velocity*g_flow_density;
                        //caution! scale by the diagonal normalization!
            force_in_OF_coord/=pow(diagonal_length,2);
            torque_in_OF_coord/=pow(diagonal_length,3);
        }
        max_force_=max(max_force_,sqrt(force_in_OF_coord.dot(force_in_OF_coord)));
        max_torque_=max(max_torque_,sqrt(torque_in_OF_coord.dot(torque_in_OF_coord)));

        Matrix3d model_rotation_in_world=AngleAxisd(rotation_angle/180*M_PI, rotation_axis).toRotationMatrix();
        Matrix3d coord_rotation_fixed_model=model_rotation_in_world.transpose();
        Vector3d force_in_world, torque_in_world;
        force_in_world=coord_rotation_fixed_model*force_in_OF_coord;
        torque_in_world=coord_rotation_fixed_model*torque_in_OF_coord;
        double alpha_radians,beta_radians;
        AlphaBetaRadiansFromRotMat(coord_rotation_fixed_model,alpha_radians,beta_radians);

        double alpha_degree=alpha_radians*180/M_PI;
        double beta_degree=beta_radians*180/M_PI;
        datapair.push_back(alpha_degree);
        datapair.push_back(beta_degree);
        datapair.push_back(diagonal_length);
        datapair.push_back(velocity);
        datapair.push_back(force_in_world[0]);
        datapair.push_back(force_in_world[1]);
        datapair.push_back(force_in_world[2]);
        datapair.push_back(torque_in_world[0]);
        datapair.push_back(torque_in_world[1]);
        datapair.push_back(torque_in_world[2]);
        dataset_.push_back(datapair);
        for(int j=0;j<datapair.size();j++)
        {
            cout<<datapair[j]<<" ";
        }cout<<endl;
    }
}

void OFInterpolator::GenerateOneModelList()
{
    //data augmentation: duplicate data using symmetry
    vector<vector<double>> new_data;


    if(g_model_name=="disk.obj"||g_model_name=="foamdisk.obj"
        ||g_model_name=="semisphereshell.obj"||g_model_name=="foamsemisphere.obj")
    {
        if(g_model_name=="disk.obj"||g_model_name=="foamdisk.obj")
        {
            int origin_datasize=dataset_.size();
            for(int i=0;i<origin_datasize;i++)
            {
                vector<double> new_datapair=dataset_[i];
                if(abs(new_datapair[1])<5) continue;
                new_datapair[1]=-new_datapair[1];
                new_datapair[input_size_+1]=-new_datapair[input_size_+1];
                new_datapair[input_size_+2]=-new_datapair[input_size_+2];
                new_datapair[input_size_+4]=-new_datapair[input_size_+4];
                new_datapair[input_size_+5]=-new_datapair[input_size_+5];
                dataset_.push_back(new_datapair);
            }
        }

        double augment_angle_gap=10;
        for(int i=0;i<dataset_.size();i++)
        {
            //for disk, only beta works for different values
            for(int j=-180;j<180;j+=augment_angle_gap)
            {
                vector<double> new_datapair=dataset_[i];
                //if(new_datapair[1]>87/*||new_datapair[1]<-85*/) continue;
                //rotate force and torque alpha the second time (because it's already in world coord)
                double relative_angle_radians=(j-new_datapair[0])*M_PI/180.0;
                Matrix3d force_torque_rotation;
                RotMatFromAlphaBeta(force_torque_rotation,relative_angle_radians,0);
                Vector3d new_force=force_torque_rotation*Vector3d(new_datapair[input_size_],new_datapair[input_size_+1],new_datapair[input_size_+2]);
                Vector3d new_torque=force_torque_rotation*Vector3d(new_datapair[input_size_+3],new_datapair[input_size_+4],new_datapair[input_size_+5]);
                
                new_datapair[0]=j;
                new_datapair[input_size_+0]=new_force[0];
                new_datapair[input_size_+1]=new_force[1];
                new_datapair[input_size_+2]=new_force[2];
                new_datapair[input_size_+3]=new_torque[0];
                new_datapair[input_size_+4]=new_torque[1];
                new_datapair[input_size_+5]=new_torque[2];
                new_data.push_back(new_datapair);
            }
        }
        dataset_.clear();
        for(int i=0;i<new_data.size();i++)
            dataset_.push_back(new_data[i]);
    }
    if(g_model_name=="openbox.obj"||g_model_name=="realbox.obj"||g_model_name=="leaflike.obj")
    {
        //alpha 0-90, beta 0-90
        //swap the value of axis when it's the symmetry axis
        for(int i=0;i<dataset_.size();i++)
        {//rotate 90 degree for 3 times along y-axis
            vector<double> new_datapair=dataset_[i];
            for(int a=0;a<4;a++)
            {
                double rotation_degree=a*90;
                Matrix3d rotation=AngleAxisd(rotation_degree*M_PI/180.0,Vector3d(0,1,0)).toRotationMatrix();
                new_datapair[0]+=rotation_degree;
                Vector3d new_force=rotation*Vector3d(new_datapair[input_size_],new_datapair[input_size_+1],new_datapair[input_size_+2]);
                Vector3d new_torque=rotation*Vector3d(new_datapair[input_size_+3],new_datapair[input_size_+4],new_datapair[input_size_+5]);
                new_datapair[input_size_+0]=new_force[0];
                new_datapair[input_size_+1]=new_force[1];
                new_datapair[input_size_+2]=new_force[2];
                new_datapair[input_size_+3]=new_torque[0];
                new_datapair[input_size_+4]=new_torque[1];
                new_datapair[input_size_+5]=new_torque[2];
                new_data.push_back(new_datapair);
            }
        }
        dataset_.clear();
        for(int i=0;i<new_data.size();i++)
            dataset_.push_back(new_data[i]);
    }
    OutputDataset();
    
}

void OFInterpolator::OutputArrowObj(int mode)
{
    double arrow_radius=0.01*g_model_length;
    double force_arrow_length_rate=1/max_force_*g_model_length*0.01;
    double torque_arrow_length_rate=1/max_torque_*g_model_length*0.01;

    //input arrow Model
    Mesh arrow_mesh, force_mesh, torque_mesh;
    OpenMesh::IO::read_mesh(arrow_mesh, g_project_path+"/"+g_model_folder+"/arrow.obj");
    force_mesh=g_model; torque_mesh=g_model;
    //add arrow
    for(int i=0;i<dataset_.size();i++)
    {
        double alpha_radians=dataset_[i][0]*M_PI/180;
        double beta_radians=dataset_[i][1]*M_PI/180;
        double alpha_degree=dataset_[i][0];
        double beta_degree=dataset_[i][1];

        double x,y,z;
        ViewXYZFromAlphaBeta(x,y,z,alpha_radians,beta_radians);
        Vector3d translate_vector(x,y,z);

        Vector3d force(dataset_[i][input_size_],dataset_[i][input_size_+1],dataset_[i][input_size_+2]);
        Vector3d force_scale_vector(sqrt(force.dot(force))*force_arrow_length_rate,arrow_radius, arrow_radius);//origin arrow along x-axis
        Vector3d torque(dataset_[i][input_size_+3],dataset_[i][input_size_+4],dataset_[i][input_size_+5]);
        Vector3d torque_scale_vector(sqrt(torque.dot(torque))*torque_arrow_length_rate,arrow_radius, arrow_radius);//origin arrow along x-axis
        
        //rotate arrow to force direction, using angle&axis
        Vector3d x_axis(1,0,0);
        x_axis.normalize();
        force.normalize();
        double force_rot_angle_radians=acos(force.dot(x_axis));//always be positive, 0~M_PI
        Vector3d force_rot_axis=x_axis.cross(force);
        force_rot_axis.normalize();
        AngleAxisd force_t_V(force_rot_angle_radians,force_rot_axis);
        x_axis=Vector3d(1,0,0);
        x_axis.normalize();
        torque.normalize();
        double torque_rot_angle_radians=acos(torque.dot(x_axis));
        Vector3d torque_rot_axis=x_axis.cross(torque);
        torque_rot_axis.normalize();
        AngleAxisd torque_t_V(torque_rot_angle_radians,torque_rot_axis);

        Mesh this_arrow=arrow_mesh;
        VertexIter v_it=this_arrow.vertices_begin();
        VertexIter v_end=this_arrow.vertices_end();
        for(;v_it!=v_end;v_it++)
        {//index start from 0
            //arrow point movement
            OpenMesh::Vec3d new_point=this_arrow.point(v_it);
            Vector3d force_arrow_point
                (new_point[0]*force_scale_vector[0],new_point[1]*force_scale_vector[1],new_point[2]*force_scale_vector[2]);
            Vector3d torque_arrow_point
                (new_point[0]*torque_scale_vector[0],new_point[1]*torque_scale_vector[1],new_point[2]*torque_scale_vector[2]);
            

            force_arrow_point=force_t_V*force_arrow_point+translate_vector;
            torque_arrow_point=torque_t_V*torque_arrow_point+translate_vector;
            OpenMesh::Vec3d force_p(force_arrow_point[0],force_arrow_point[1],force_arrow_point[2]);
            OpenMesh::Vec3d torque_p(torque_arrow_point[0],torque_arrow_point[1],torque_arrow_point[2]);
            force_mesh.add_vertex(force_p);
            torque_mesh.add_vertex(torque_p);
        }

        FaceIter f_it=this_arrow.faces_begin();
        FaceIter f_end=this_arrow.faces_end();
        for(;f_it!=f_end;f_it++)
        {
            FaceVertexIter fv_it=this_arrow.fv_begin(f_it);
            FaceVertexIter fv_end=this_arrow.fv_end(f_it);
            std::vector<Mesh::VertexHandle> face_vhs_force;
            std::vector<Mesh::VertexHandle> face_vhs_torque;
            for(;fv_it!=fv_end;fv_it++)
            {
                int arrow_v_index=fv_it->idx();
                VertexHandle vh_force=force_mesh.vertex_handle(arrow_v_index+force_mesh.n_vertices()-this_arrow.n_vertices());
                VertexHandle vh_torque=torque_mesh.vertex_handle(arrow_v_index+torque_mesh.n_vertices()-this_arrow.n_vertices());
                face_vhs_force.push_back(vh_force);
                face_vhs_torque.push_back(vh_torque);
            }
            force_mesh.add_face(face_vhs_force);
            torque_mesh.add_face(face_vhs_torque);
        }

    }
    if(mode==0)
    {
        OpenMesh::IO::write_mesh(force_mesh,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name+"_f_part.obj");
        OpenMesh::IO::write_mesh(torque_mesh,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name+"_t_part.obj");
    }
    else
    {
        OpenMesh::IO::write_mesh(force_mesh,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name+"_force_vis.obj");
        OpenMesh::IO::write_mesh(torque_mesh,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name+"_torque_vis.obj");
    }
    
}

void OFInterpolator::OutputCsvFile(int mode)//0 for origin data, 1 for duplicated data by symmetry
{

    string new_data_case_name=g_model_name.substr(0,g_model_name.length()-4);
    string command = "cd " + g_project_path+"/"+g_dataset_folder+"/ModelsData;";
    command += "mkdir " + new_data_case_name +";";
    system(command.c_str());

    g_model_length=g_model_length/g_diagonal_length;
    g_of_processor.NormalizeModel();
    OpenMesh::IO::write_mesh(g_model,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+new_data_case_name+"/"+g_model_name);

    string out_filename;
    if(mode==0)
        out_filename=g_project_path+"/"+g_dataset_folder+"/ModelsData/"+new_data_case_name+"/"+g_model_name+".csv";
    if(mode==1)
        out_filename=g_project_path+"/"+g_dataset_folder+"/ModelsData/"+new_data_case_name+"/Aug.csv";
    ofstream fout(out_filename);
    for(int i=0;i<dataset_.size();i++)
    {
        double alpha_radians=dataset_[i][0]*M_PI/180;
        double beta_radians=dataset_[i][1]*M_PI/180;
        double x,y,z;
        ViewXYZFromAlphaBeta(x,y,z,alpha_radians,beta_radians);
        fout<<x<<","<<y<<","<<z<<",";

        for(int j=4;j<dataset_[i].size();j++)
        {
            fout<<dataset_[i][j];
            if(j!=dataset_[i].size()-1) fout<<",";
            else fout<<endl;
        }

        if(mode==0) continue;
        Vector3d normal(x,y,z);
        Vector3d force(dataset_[i][input_size_+0],dataset_[i][input_size_+1],dataset_[i][input_size_+2]);
        double force_dir=force.dot(normal);
        if(force_dir<0)
            cout<<"wrong force dirc!!"<<force_dir<<"  "<<dataset_[i][0]<<"||"<<dataset_[i][1]<<endl;
    }
    fout.close();
}

void OFInterpolator::OutputSHParameters()
{
    int longitudinal_div=36;
    int latitudinal_div=36;
    if(g_model_name=="foamdisk.obj")
        latitudinal_div=34;
    Eigen::Vector3d up(0,1,0);
    int SH_order=9;
    
    int SH_num=(SH_order+1)*(SH_order+1);
    vector<double> SH_coeff_F(SH_num*3,0);
    vector<double> SH_coeff_T(SH_num*3,0);
    double delta_alpha_radians=2*M_PI/(double)longitudinal_div;
    double delta_beta_radians=M_PI/(double)latitudinal_div;

    //for each point instead of polar point
        //compute the small surface area with one value in the corner
    for(int sh=0;sh<SH_num;sh++)
    {
        //start intergration for one parameter
        for(int i=0;i<dataset_.size();i++)
        {
            double alpha_radians=dataset_[i][0]*M_PI/180;
            double beta_radians=dataset_[i][1]*M_PI/180;
            double x,y,z;
            ViewXYZFromAlphaBeta(x,y,z,alpha_radians,beta_radians);

            if(abs(dataset_[i][1])<2) continue;
            double small_area=SmallSphereAreaPartFromBeta(beta_radians,delta_alpha_radians,delta_beta_radians);
            
            double SH_a[100]; makeArray_SphericalHarmonics(SH_a, SH_order, x, y, z);
            Vector3d dir_vec(x,y,z);
            Vector3d force_coeff_vec(dataset_[i][input_size_+0],dataset_[i][input_size_+1],dataset_[i][input_size_+2]);
            Vector3d torque_coeff_vec(dataset_[i][input_size_+3],dataset_[i][input_size_+4],dataset_[i][input_size_+5]);
            for(int dim=0;dim<=2;dim++)
            {
                // if(dim==0) dir_vec=dir_vec;             //drag value
                // if(dim==1) dir_vec=dir_vec.cross(up);   //lift value
                // if(dim==2)                              //side value
                // {
                //     Vector3d dir_1=dir_vec.cross(up);
                //     dir_vec=dir_vec.cross(dir_1);
                // }
                // dir_vec.normalize();
                // double force_coeff=force_coeff_vec.dot(dir_vec);
                // double torque_coeff=torque_coeff_vec.dot(dir_vec);
                double force_coeff=force_coeff_vec[dim];
                double torque_coeff=torque_coeff_vec[dim];

                SH_coeff_F[sh*3+dim]+=force_coeff*SH_a[sh]*small_area;
                SH_coeff_T[sh*3+dim]+=torque_coeff*SH_a[sh]*small_area;
            }
        }
    }    
    //save txt file for SH coeff
    string obj_name=g_model_name.substr(0,g_model_name.length()-4);
    string fpath=g_project_path+"/"+g_dataset_folder+"/ModelsData/"+obj_name+"/"+obj_name+".txt";
    std::fstream fout(fpath,std::ios::out);
    fout << "gravity: " << 0 << " " << -10 << " " << 0 << std::endl;
    fout << "rho_rigid: " << 5 << std::endl;
    fout << "rho_fluid: " << 1 << std::endl;
    fout << "replen: " << 0.1 << std::endl;
    fout << "Iscale: " << 0.05 << std::endl;
    fout << "Dscale: " << 0.01 << std::endl;
    fout << "shorderF: "  << SH_order << std::endl;
    for(int ish=0;ish<SH_num;++ish){
      fout << "   " << ish << " ";
      fout << SH_coeff_F[ish*3+0] << " ";
      fout << SH_coeff_F[ish*3+1] << " ";
      fout << SH_coeff_F[ish*3+2] << std::endl;
    }
    fout << "shorderT: "  << SH_order << std::endl;
    for(int ish=0;ish<SH_num;++ish){
      fout << "   " << ish << " ";
      fout << SH_coeff_T[ish*3+0] << " ";
      fout << SH_coeff_T[ish*3+1] << " ";
      fout << SH_coeff_T[ish*3+2] << std::endl;
    }

    SH_coeff_F_=SH_coeff_F;
    SH_coeff_T_=SH_coeff_T;

    //test stability
    for(int a=0;a<longitudinal_div;a++)
    {
        for(int b=0;b<latitudinal_div;b++)
        {
            std::vector<double> datapair(input_size_+6,0);
            
            double new_delta_alpha_radians=2*M_PI/(double)longitudinal_div;
            double new_delta_beta_radians=M_PI/(double)latitudinal_div;
            double new_alpha_radians=a*new_delta_alpha_radians;
            double new_beta_radians=b*new_delta_beta_radians-0.5*M_PI;
            
            Vector3d new_sample_direction;
            new_sample_direction[0]=cos(new_beta_radians)*cos(new_alpha_radians);
            new_sample_direction[1]=sin(new_beta_radians);
            new_sample_direction[2]=-cos(new_beta_radians)*sin(new_alpha_radians);
            
            //get force from SH
            double Y[100]; makeArray_SphericalHarmonics(Y,SH_order,
                new_sample_direction[0],new_sample_direction[1],new_sample_direction[2]);
            Vector3d new_sample_force(0,0,0);
            for(int sh=0;sh<SH_num;sh++)
            {
                new_sample_force[0]+=Y[sh]*SH_coeff_F[sh*3+0];
                new_sample_force[1]+=Y[sh]*SH_coeff_F[sh*3+1];
                new_sample_force[2]+=Y[sh]*SH_coeff_F[sh*3+2];
            }
            double force_direction=new_sample_force.dot(new_sample_direction);
            if(force_direction<0)
            {
                cout<<"insane accelerating point! alpha: "<<new_alpha_radians*180/M_PI
                <<"||| beta: "<<new_beta_radians*180/M_PI<<endl;
            }
        }
    }    
}

void OFInterpolator::OutputSHParameters_interpolated()
{
    int longitudinal_div=36;
    int latitudinal_div=36;
    Eigen::Vector3d up(0,1,0);
    int SH_order=9;
    
    int SH_num=(SH_order+1)*(SH_order+1);
    vector<double> SH_coeff_F(SH_num*3,0);
    vector<double> SH_coeff_T(SH_num*3,0);
    double delta_alpha_radians=2*M_PI/(double)longitudinal_div;
    double delta_beta_radians=M_PI/(double)latitudinal_div;

    for(int sh=0;sh<SH_num;sh++)
    {
        for(int lon=0;lon<longitudinal_div;lon++)
        {
            for(int lat=0;lat<latitudinal_div+1;lat++)
            {
                double new_alpha_radians=delta_alpha_radians*lon;
                double new_beta_radians=delta_beta_radians*lat-M_PI*0.5;
                int near_index=0;
                double near_dist=DBL_MAX;
                for(int i=0;i<dataset_.size();i++)
                {
                    double alpha_radians=dataset_[i][0]*M_PI/180;
                    double beta_radians=dataset_[i][1]*M_PI/180;
                    double this_dist=sqrt(pow(alpha_radians-new_alpha_radians,2)+pow(alpha_radians-new_alpha_radians,2));
                    if(this_dist<near_dist)
                    {
                        near_dist=this_dist; near_index=i;
                    }
                }

                //just copy the nearest value
                int i=near_index;
                double alpha_radians=new_alpha_radians;
                double beta_radians=new_beta_radians;
                
                double x=cos(beta_radians)*cos(alpha_radians), y=sin(beta_radians), z=-cos(beta_radians)*sin(alpha_radians);
                x=-x;y=-y;z=-z;
                
                if(abs(dataset_[i][1])<2) continue;
                //use nearby alpha-beta get a trapezoid region
                double lower_line_length=abs(2*sin(0.5*delta_alpha_radians)*cos(beta_radians));
                double smaller_beta_radians=beta_radians;
                if(beta_radians>0) smaller_beta_radians=beta_radians-delta_beta_radians;
                if(beta_radians<0) smaller_beta_radians=beta_radians+delta_beta_radians;
                double upper_line_length=abs(2*sin(0.5*delta_alpha_radians)*cos(smaller_beta_radians));
                double height=abs(2*sin(0.5*delta_beta_radians));
                double small_area=(lower_line_length+upper_line_length)*0.5*height;
                
                double SH_a[100]; makeArray_SphericalHarmonics(SH_a, SH_order, x, y, z);
                Vector3d dir_vec(x,y,z);
                Vector3d force_coeff_vec(dataset_[i][input_size_+0],dataset_[i][input_size_+1],dataset_[i][input_size_+2]);
                Vector3d torque_coeff_vec(dataset_[i][input_size_+3],dataset_[i][input_size_+4],dataset_[i][input_size_+5]);
                for(int dim=0;dim<=2;dim++)
                {
                    if(dim==0) dir_vec=dir_vec;             //drag value
                    if(dim==1) dir_vec=dir_vec.cross(up);   //lift value
                    if(dim==2)                              //side value
                    {
                        Vector3d dir_1=dir_vec.cross(up);
                        dir_vec=dir_vec.cross(dir_1);
                    }
                    dir_vec.normalize();
                    double force_coeff=force_coeff_vec.dot(dir_vec);
                    double torque_coeff=torque_coeff_vec.dot(dir_vec);
                    force_coeff=force_coeff_vec[dim];
                    torque_coeff=torque_coeff_vec[dim];
                    //sanity check
                    //force_coeff=up[dim];
                    //torque_coeff=up[dim];

                    SH_coeff_F[sh*3+dim]+=force_coeff*SH_a[sh]*small_area;
                    SH_coeff_T[sh*3+dim]+=torque_coeff*SH_a[sh]*small_area;
                }
            }
        }
    }    
    //save txt file for SH coeff
    string obj_name=g_model_name.substr(0,g_model_name.length()-4);
    string fpath=g_project_path+"/"+g_dataset_folder+"/ModelsData/"+obj_name+"/"+obj_name+".txt";
    std::fstream fout(fpath,std::ios::out);
    fout << "gravity: " << 0 << " " << -10 << " " << 0 << std::endl;
    fout << "rho_rigid: " << 5 << std::endl;
    fout << "rho_fluid: " << 1 << std::endl;
    fout << "replen: " << 0.1 << std::endl;
    fout << "Iscale: " << 0.05 << std::endl;
    fout << "Dscale: " << 0.01 << std::endl;
    fout << "shorderF: "  << SH_order << std::endl;
    for(int ish=0;ish<SH_num;++ish){
      fout << "   " << ish << " ";
      fout << SH_coeff_F[ish*3+0] << " ";
      fout << SH_coeff_F[ish*3+1] << " ";
      fout << SH_coeff_F[ish*3+2] << std::endl;
    }
    fout << "shorderT: "  << SH_order << std::endl;
    for(int ish=0;ish<SH_num;++ish){
      fout << "   " << ish << " ";
      fout << SH_coeff_T[ish*3+0] << " ";
      fout << SH_coeff_T[ish*3+1] << " ";
      fout << SH_coeff_T[ish*3+2] << std::endl;
    }

    SH_coeff_F_=SH_coeff_F;
    SH_coeff_T_=SH_coeff_T;

    //test stability
    for(int a=0;a<longitudinal_div;a++)
    {
        for(int b=0;b<latitudinal_div;b++)
        {
            std::vector<double> datapair(input_size_+6,0);
            
            double new_delta_alpha_radians=2*M_PI/(double)longitudinal_div;
            double new_delta_beta_radians=M_PI/(double)latitudinal_div;
            double new_alpha_radians=a*new_delta_alpha_radians;
            double new_beta_radians=b*new_delta_beta_radians-0.5*M_PI;
            
            Vector3d new_sample_direction;
            new_sample_direction[0]=cos(new_beta_radians)*cos(new_alpha_radians);
            new_sample_direction[1]=sin(new_beta_radians);
            new_sample_direction[2]=-cos(new_beta_radians)*sin(new_alpha_radians);
            
            //get force from SH
            double Y[100]; makeArray_SphericalHarmonics(Y,SH_order,
                new_sample_direction[0],new_sample_direction[1],new_sample_direction[2]);
            Vector3d new_sample_force(0,0,0);
            for(int sh=0;sh<SH_num;sh++)
            {
                new_sample_force[0]+=Y[sh]*SH_coeff_F[sh*3+0];
                new_sample_force[1]+=Y[sh]*SH_coeff_F[sh*3+1];
                new_sample_force[2]+=Y[sh]*SH_coeff_F[sh*3+2];
            }
            double force_direction=new_sample_force.dot(new_sample_direction);
            
                cout<<force_direction<<endl;
            if(force_direction<0)
            {
                cout<<"insane accelerating point! alpha: "<<new_alpha_radians*180/M_PI
                <<"||| beta: "<<new_beta_radians*180/M_PI<<endl;
            }
        }
    }    
}