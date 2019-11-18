#include "OfInterpolator/of_interpolator.h"
using namespace std;
using namespace Eigen;

//data format-index
//rotation_angle, axis, length, area, velo, Re, Cd, Cl, Cs, coeffM, force, torque, avg_normal
//0               123   4       5     6     7   8   9   10  111213  141516 171819  202122


void OFInterpolator::GetDatasetForOneModel()
{
    dataset_.clear();
    input_size_=4;//alpha-from_x, beta_from_xz, length, velocity, force, torque

    max_force_=0;max_torque_=0;
    for(int i=0;i<g_dataset.size();i++)
    {
        vector<double> datapair;
        double rotation_angle, length, velocity;
        Vector3d rotation_axis, force_in_OF_coord, torque_in_OF_coord;

        rotation_angle=g_dataset[i][0];
        rotation_axis=Vector3d(g_dataset[i][1],g_dataset[i][2],g_dataset[i][3]);
        length=g_dataset[i][4];
        velocity=g_dataset[i][6];

        bool output_real_force_torque=false;
        if(output_real_force_torque)
        {
            force_in_OF_coord=Vector3d(g_dataset[i][14],g_dataset[i][15],g_dataset[i][16]);
            torque_in_OF_coord=Vector3d(g_dataset[i][17],g_dataset[i][18],g_dataset[i][19]);
        }
        else
        {
            force_in_OF_coord=Vector3d(g_dataset[i][8],g_dataset[i][9],g_dataset[i][10]);
            torque_in_OF_coord=Vector3d(g_dataset[i][11],g_dataset[i][12],g_dataset[i][13]);
            force_in_OF_coord*=g_dataset[i][5];
            torque_in_OF_coord*=g_dataset[i][5]*g_dataset[i][4];        
            // //make force&torque become the Coeff Cd*A and Cm*A*l: 
            force_in_OF_coord=Vector3d(g_dataset[i][14],g_dataset[i][15],g_dataset[i][16]);
            torque_in_OF_coord=Vector3d(g_dataset[i][17],g_dataset[i][18],g_dataset[i][19]);
            force_in_OF_coord/=0.5*velocity*velocity*g_flow_density;
            torque_in_OF_coord/=0.5*velocity*velocity*g_flow_density;

        }
        max_force_=max(max_force_,sqrt(force_in_OF_coord.dot(force_in_OF_coord)));
        max_torque_=max(max_torque_,sqrt(torque_in_OF_coord.dot(torque_in_OF_coord)));

        //compute the rotation between world & flow coord, flow_coord means: flow direction=x axis
        Vector3d force_in_world, torque_in_world;
        Matrix3d model_rotation_in_world=AngleAxisd(rotation_angle/180*M_PI, rotation_axis).toRotationMatrix();
        Matrix3d coord_from_world_to_flow=model_rotation_in_world.transpose();
        Matrix3d coord_from_flow_to_world=coord_from_world_to_flow;//.transpose();????
        force_in_world=coord_from_flow_to_world*force_in_OF_coord;
        cout<<"force in of:"<<force_in_OF_coord<<endl;
        cout<<"force in world:"<<force_in_world<<endl;
        torque_in_world=coord_from_flow_to_world*torque_in_OF_coord;
        Vector3d sphere_point(1,0,0);
        sphere_point=coord_from_world_to_flow*sphere_point;//????

        double alpha, beta;
        alpha=acos(sphere_point[0]/sqrt(sphere_point[0]*sphere_point[0]+sphere_point[2]*sphere_point[2]));
        if(sphere_point[2]<0)alpha=-alpha;
        beta=asin(sphere_point[1]/1.0);
        alpha*=180/M_PI;beta*=180/M_PI;
        datapair.push_back(alpha);
        datapair.push_back(beta);
        datapair.push_back(length);
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


    if(g_model_name=="disk.obj")
    {
    double augment_angle_gap=20;
    for(int i=0;i<dataset_.size();i++)
    {
        //for disk, only beta works for different values
        for(int j=-180;j<180;j+=augment_angle_gap)
        {
            vector<double> new_datapair=dataset_[i];
            if(new_datapair[1]>89) continue;
            //rotate force and torque
            double relative_angle=new_datapair[0]-j;
            AngleAxisd t_V(relative_angle/180*M_PI,Vector3d(0,1,0));
            Vector3d new_force=t_V*Vector3d(new_datapair[input_size_],new_datapair[input_size_+1],new_datapair[input_size_+2]);
            Vector3d new_torque=t_V*Vector3d(new_datapair[input_size_+3],new_datapair[input_size_+4],new_datapair[input_size_+5]);
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
    for(int i=0;i<new_data.size();i++)
        dataset_.push_back(new_data[i]);
    }

    else
    {
        //alpha 0-90, beta 0-90
        //swap the value of axis when it's the symmetry axis
    for(int i=0;i<dataset_.size();i++)
    {//x
        vector<double> new_datapair=dataset_[i];
        if(new_datapair[0]==0||new_datapair[0]==-180||new_datapair[0]==180) continue;
        new_datapair[0]=-new_datapair[0];
        new_datapair[input_size_+0]=-new_datapair[input_size_+0];
        new_datapair[input_size_+3]=-new_datapair[input_size_+3];
        new_data.push_back(new_datapair);
    }
    for(int i=0;i<new_data.size();i++)
        dataset_.push_back(new_data[i]);
    for(int i=0;i<dataset_.size();i++)
    {//z
        vector<double> new_datapair=dataset_[i];
        if(new_datapair[0]==90||new_datapair[0]==-90) continue;
        if(new_datapair[0]>0)
            new_datapair[0]=180-new_datapair[0];
        else
            new_datapair[0]=-180-new_datapair[0];
        new_datapair[input_size_+2]=-new_datapair[input_size_+2];
        new_datapair[input_size_+5]=-new_datapair[input_size_+5];
        new_data.push_back(new_datapair);
    }
    for(int i=0;i<new_data.size();i++)
        dataset_.push_back(new_data[i]);
    // for(int i=0;i<dataset_.size();i++)
    // {//y
    //     vector<double> new_datapair=dataset_[i];
    //     if(new_datapair[1]==0) continue;
    //         new_datapair[1]=-new_datapair[1];
    //     new_datapair[input_size_+1]=-new_datapair[input_size_+1];
    //     new_datapair[input_size_+4]=-new_datapair[input_size_+4];
    //     new_data.push_back(new_datapair);
    // }
    // for(int i=0;i<new_data.size();i++)
    //     dataset_.push_back(new_data[i]);  

    }
    OutputDataset();
    
}

void OFInterpolator::OutputArrowObj()
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
        double alpha=dataset_[i][0]*M_PI/180;
        double beta=dataset_[i][1]*M_PI/180;
        Vector3d translate_vector(cos(alpha)*cos(beta)*g_model_length,sin(beta)*g_model_length,sin(alpha)*cos(beta)*g_model_length);
        Vector3d force(dataset_[i][input_size_],dataset_[i][input_size_+1],dataset_[i][input_size_+2]);
        Vector3d force_scale_vector(sqrt(force.dot(force))*force_arrow_length_rate,arrow_radius, arrow_radius);//origin arrow along x-axis
        Vector3d torque(dataset_[i][input_size_+3],dataset_[i][input_size_+4],dataset_[i][input_size_+5]);
        Vector3d torque_scale_vector(sqrt(torque.dot(torque))*torque_arrow_length_rate,arrow_radius, arrow_radius);//origin arrow along x-axis
        Vector3d x_axis(1,0,0);
        force.normalize();
        double force_rot_angle=force.dot(x_axis);
        Vector3d force_rot_axis=x_axis.cross(force);
        AngleAxisd force_t_V(acos(force_rot_angle),force_rot_axis);
        torque.normalize();
        double torque_rot_angle=torque.dot(x_axis);
        Vector3d torque_rot_axis=x_axis.cross(torque);
        AngleAxisd torque_t_V(acos(torque_rot_angle),torque_rot_axis);
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
    OpenMesh::IO::write_mesh(force_mesh,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name+"_force_vis.obj");
    OpenMesh::IO::write_mesh(torque_mesh,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name+"_torque_vis.obj");
}

void OFInterpolator::OutputCsvFile()
{
    string out_filename=g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name+".txt";
    ofstream fout(out_filename);
    for(int i=0;i<dataset_.size();i++)
    {
        double alpha=dataset_[i][0]*M_PI/180;
        double beta=dataset_[i][1]*M_PI/180;
        double x=cos(beta)*cos(alpha), y=sin(beta), z=cos(beta)*sin(alpha);
        fout<<x<<","<<y<<","<<z<<",";
        for(int j=2;j<dataset_[i].size();j++)
        {
            fout<<dataset_[i][j];
            if(j!=dataset_[i].size()-1) fout<<",";
            else fout<<endl;
        }
    }
    fout.close();

    OpenMesh::IO::write_mesh(g_model,g_project_path+"/"+g_dataset_folder+"/ModelsData/"+g_model_name);
    // string command = "cd " + g_project_path + ";";
    // command += "cp " + g_model_folder+"/"+g_model_name + " " + g_dataset_folder+"/ModelsData;";
    // system(command.c_str());
}