#include "OfInterpolator/of_interpolator.h"
using namespace std;
using namespace OpenMesh;
using namespace Eigen;
//data format-index
//rotation_angle, axis, length, area, velo, Re, Cd, Cl, Cs, coeffM, force, torque, avg_normal
//0               123   4       5     6     7   8   9   10  111213  141516 171819  202122
void OFInterpolator::Init()
{
    data_element_count_=23;
    ReadInputPair();
    
    if(g_program_mode==2)
    {
        //judge if need compute physical quantities
        ifstream fin(g_project_path+"/"+g_model_folder+"/Quantities/"+g_model_name+".qu");
        if(!fin)
            Voxelize();

        ChangeDataFormat();

        //flag of interpolators
        linear_computed_=false;
    }
    if(g_program_mode==3)
    {
        GetDatasetForOneModel();
        OutputCsvFile(0);
        OutputArrowObj(0);
        
        GenerateOneModelList();
        OutputArrowObj(1);
        OutputCsvFile(1);

        OutputSHParameters();
    }
    if(g_program_mode==5)
    {
        ResetDatasetNormalArea();
        ChangeDataFormat();
        AugmentNormalData();
        InterpolateForNormalMap();
        OutputArrowObj(1);
        OutputSHParameters();
    }
}

void OFInterpolator::ReadInputPair()
{
    g_dataset.clear();
    string filename=g_project_path+"/"+g_dataset_folder+"/"+g_model_name+".TLH";
    ifstream dataset_in(filename);
    if (dataset_in.fail() == true)
    {
        cout << "Cannot load dataset files." << endl;
        return;
    }
    int data_element_count=0;
    double current_data_element;
    vector<double> current_data_line;
    while(dataset_in>>current_data_element)
    {
        current_data_line.push_back(current_data_element);
        if(data_element_count==data_element_count_-1)
        {
            g_dataset.push_back(current_data_line);
            current_data_line.clear();
            data_element_count=0;
            data_line_count_++;
        }
        else data_element_count++;
    }
    cout<<"Load dataset successfully. data_line_num: "<<g_dataset.size()<<endl;
    cout<<"                           data_element_num: "<<g_dataset[0].size()<<endl;
}

double ComputeRadiansFrom2DVec(double x, double y)
{
    x/=sqrt(x*x+y*y);
    double theta=acos(x);
    if(y<0) theta=-theta;
    return theta;
}

void OFInterpolator::ReadNormalList()
{
    g_normal_list.clear();
    std::string NL_filename="../Render/NormalList/"+g_model_name+".NL";
    std::fstream fin(NL_filename,std::ios::in);
    double param;
    while(fin>>param)
    {
        std::vector<double> normal_data(6,0);
        normal_data[0]=param;
        for(int i=1;i<6;i++)
            fin>>normal_data[i];
        g_normal_list.push_back(normal_data);
    }
}

void OFInterpolator::ResetDatasetNormalArea()
{
    ReadNormalList();
    for(int i=0;i<g_dataset.size();i++)
    {
        Vector3d rotation_axis(g_dataset[i][1],g_dataset[i][2],g_dataset[i][3]);
        Matrix3d model_rotation_in_world=AngleAxisd(g_dataset[i][0]/180*M_PI, rotation_axis).toRotationMatrix();
        double alpha_radians, beta_radians;
        AlphaBetaRadiansFromRotMat(model_rotation_in_world.transpose(),alpha_radians,beta_radians);
        double alpha_degree=alpha_radians*180/M_PI;
        double beta_degree=beta_radians*180/M_PI;
        
        double min_dist=DBL_MAX;
        int min_index=0;
        for(int n=0;n<g_normal_list.size();n++)
        {
            double this_dist=pow(g_normal_list[n][0]-alpha_degree,2)+pow(g_normal_list[n][1]-beta_degree,2);
            if(this_dist<min_dist)
            {
                min_dist=this_dist;
                min_index=n;
            }
        }
        g_dataset[i][20]=g_normal_list[min_index][3];
        g_dataset[i][21]=g_normal_list[min_index][4];
        g_dataset[i][22]=g_normal_list[min_index][5];
    }
}

void OFInterpolator::ChangeDataFormat()
{
    dataset_.clear();
    for(int i=0;i<g_dataset.size();i++)
    {   
        Eigen::Vector3d local_x;
        Eigen::Vector3d local_y;
        Eigen::Vector3d local_z;
        Eigen::Vector3d f_coeff(g_dataset[i][8],g_dataset[i][9],g_dataset[i][10]);
        f_coeff*=g_dataset[i][5];//save the Cd*A value
        Eigen::Vector3d t_coeff(g_dataset[i][11],g_dataset[i][12],g_dataset[i][13]);
        t_coeff*=g_dataset[i][5]*g_dataset[i][4];//save the Cm*A*l value
        Eigen::Vector3d normal_area(g_dataset[i][20],g_dataset[i][21],g_dataset[i][22]);
        normal_area*=g_dataset[i][5];//save the normal*area value
        Vector3d rotation_axis(g_dataset[i][1],g_dataset[i][2],g_dataset[i][3]);
        Eigen::Matrix3d model_rotation_in_world=AngleAxisd(g_dataset[i][0]/180.0*M_PI, rotation_axis).toRotationMatrix();
        Eigen::Matrix3d viewpoint_rotation=model_rotation_in_world.transpose();
        
        local_z=-viewpoint_rotation*Vector3d(1,0,0);
        local_y=viewpoint_rotation*Vector3d(0,1,0);
        local_x=viewpoint_rotation*Vector3d(0,0,1);


        std::vector<double> datapair;
        input_size_=12;
        datapair.push_back(local_x[0]); datapair.push_back(local_x[1]); datapair.push_back(local_x[2]);
        datapair.push_back(local_y[0]); datapair.push_back(local_y[1]); datapair.push_back(local_y[2]);
        datapair.push_back(local_z[0]); datapair.push_back(local_z[1]); datapair.push_back(local_z[2]);
        datapair.push_back(normal_area[0]); datapair.push_back(normal_area[1]); datapair.push_back(normal_area[2]);
        datapair.push_back(f_coeff[0]); datapair.push_back(f_coeff[1]); datapair.push_back(f_coeff[2]);
        datapair.push_back(t_coeff[0]); datapair.push_back(t_coeff[1]); datapair.push_back(t_coeff[2]);
        dataset_.push_back(datapair);
    }
}

void OFInterpolator::AugmentNormalData()//rotation by the flow-axis
{
    vector<vector<double>> new_dataset;
    for(int i=0;i<dataset_.size();i++)
    {
        Vector3d local_x(dataset_[i][0],dataset_[i][1],dataset_[i][2]);
        Vector3d local_y(dataset_[i][3],dataset_[i][4],dataset_[i][5]);
        Vector3d local_z(dataset_[i][6],dataset_[i][7],dataset_[i][8]);
        Vector3d normal_area(dataset_[i][9],dataset_[i][10],dataset_[i][11]);
        Vector3d f_coeff(dataset_[i][12],dataset_[i][13],dataset_[i][14]);
        Vector3d t_coeff(dataset_[i][15],dataset_[i][16],dataset_[i][17]);

        int rot_num=36;
        for(int j=0;j<rot_num;j++)
        {
            vector<double> datapair;
            double rot_radians=2*M_PI/(double)rot_num*j;
            Matrix3d world_rot=AngleAxisd(rot_radians,local_z).toRotationMatrix();
            local_x=world_rot*local_x;
            local_y=world_rot*local_y;
            Vector3d local_rot_axis(0,0,1);
            Matrix3d local_rot=AngleAxisd(rot_radians,local_rot_axis).toRotationMatrix();
            normal_area=local_rot*normal_area;
            f_coeff=local_rot*f_coeff;
            t_coeff=local_rot*t_coeff;
            datapair.push_back(local_x[0]); datapair.push_back(local_x[1]); datapair.push_back(local_x[2]);
            datapair.push_back(local_y[0]); datapair.push_back(local_y[1]); datapair.push_back(local_y[2]);
            datapair.push_back(local_z[0]); datapair.push_back(local_z[1]); datapair.push_back(local_z[2]);
            datapair.push_back(normal_area[0]); datapair.push_back(normal_area[1]); datapair.push_back(normal_area[2]);
            datapair.push_back(f_coeff[0]); datapair.push_back(f_coeff[1]); datapair.push_back(f_coeff[2]);
            datapair.push_back(t_coeff[0]); datapair.push_back(t_coeff[1]); datapair.push_back(t_coeff[2]);
            new_dataset.push_back(datapair);
        }
    }
    dataset_=new_dataset;
    
    // for(int i=0;i<dataset_.size();i++)
    // {
    //     for(int j=0;j<dataset_[0].size();j++)
    //     {
    //         cout<<dataset_[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
}


void OFInterpolator::InterpolateForNormalMap()
{

    //delete x,y,z lines in dataset, only keep normal-f/t
    input_size_=3;
    int origin_dataset_size=dataset_[0].size();
    for(int i=0;i<dataset_.size();i++)
    {
        vector<double> new_data;
        for(int j=9;j<origin_dataset_size;j++)
            new_data.push_back(dataset_[i][j]);
        dataset_[i]=new_data;
    }
    
    g_model_name=g_model_name;
    ReadNormalList();
    
    ComputeLeastSquareParams();

    vector<vector<double>> new_dataset;
    for(int i=0;i<g_normal_list.size();i++)
    {
        vector<double> new_data(8,0);//alpha,beta, normal*area, force, torque
        new_data[0]=g_normal_list[i][0];
        new_data[1]=g_normal_list[i][1];
        vector<double> interpolate_data(9,0);
        Vector3d normal_area(g_normal_list[i][0],g_normal_list[i][1],g_normal_list[i][2]);
        for(int a=0;a<3;a++)
            interpolate_data[a]=normal_area[a];
        LeastSquareFitting(interpolate_data);

        for(int a=0;a<6;a++)
            new_data[a+2]=interpolate_data[a+3];
        new_dataset.push_back(new_data);
    }
    input_size_=2;
    dataset_=new_dataset;

    for(int i=0;i<dataset_.size();i++)
    {
        for(int j=0;j<dataset_[0].size();j++)
        {
            cout<<dataset_[i][j]<<" ";
        }
        cout<<endl;
    }
}









// OD: ray, p0p1p2: triangle
bool ray_triangle_intersect_mt(const Vec3f O, const Vec3f D, const Vec3f p0, const Vec3f p1, const Vec3f p2, float &t, float &u, float &v){
  Vec3f e1 = p1-p0;
  Vec3f e2 = p2-p0;
  Vec3f pvec = cross(D, e2);
  float det = dot(e1,pvec);
  if (det < 1.0e-9) return false;

  Vec3f tvec = O - p0;
  u = dot(tvec,pvec*(1./det));
  if (u < 0 || u > 1) return false;

  Vec3f qvec = cross(tvec, e1);
  v = dot(D,qvec*(1./det));
  if (v < 0 || u + v > 1) return false;

  t = dot(e2,qvec*(1./det));
  return t > 1.0e-9;
}

void OFInterpolator::Voxelize()
{
    int voxelize_res=20;
    int voxel_count=0;
    voxel_res_=voxelize_res;
    gap_=g_model_length/(float)voxelize_res;
    voxel_model_=vector<vector<vector<int>>>(voxelize_res,vector<vector<int>>(voxelize_res,vector<int>(voxelize_res,0)));
    for(int i=0;i<voxelize_res;i++)
    {
        cout<<"voxelization: "<<i<<"/"<<voxelize_res<<endl;
        for(int j=0;j<voxelize_res;j++)
        {
            for(int k=0;k<voxelize_res;k++)
            {
                Vec3f voxel;
                voxel[0]=g_model_length/(float)voxelize_res*(i+0.5)-0.5*g_model_length;
                voxel[1]=g_model_length/(float)voxelize_res*(j+0.5)-0.5*g_model_length;
                voxel[2]=g_model_length/(float)voxelize_res*(k+0.5)-0.5*g_model_length;
                int intersect_count=0;
                for(int f=0;f<g_model.n_faces();f++)
                {
                    FaceHandle fh=g_model.face_handle(f);
                    FaceVertexIter fv_it=g_model.fv_begin(fh);
                    FaceVertexIter fv_it_end=g_model.fv_end(fh);
                    vector<Vec3f> triangle;
                    for(;fv_it!=fv_it_end;fv_it++)
                    {
                        VertexHandle vh=fv_it.handle();
                        Point point= g_model.point(vh);
                        triangle.push_back(Vec3f(point[0],point[1],point[2]));
                    }
                    float t, u, v;
                    if(ray_triangle_intersect_mt(voxel,Vec3f(0,2*g_model_length,3*g_model_length),triangle[1],triangle[0],triangle[2],t,u,v)
                        ||ray_triangle_intersect_mt(voxel,Vec3f(0,2*g_model_length,3*g_model_length),triangle[0],triangle[1],triangle[2],t,u,v))
                        intersect_count++;
                }
                if(intersect_count%2==1)
                {
                    voxel_model_[i][j][k]=1;
                    voxel_count++;
                }
            }
        }
    }
    voxel_to_obj();
    ComputePhysicalQuantities();
}

void OFInterpolator::ComputePhysicalQuantities()
{
    //compute centroid
    //compute inertia tensor    *******this inertia should multiply the density
    //compute volume

    Vec3f centroid(0,0,0);
    glm::mat3 inertia_tensor_wo_density(0);
    float volume=0;
    int voxel_count=0;
    int voxelize_res=voxel_res_;
    for(int i=0;i<voxelize_res;i++)
    {
        for(int j=0;j<voxelize_res;j++)
        {
            for(int k=0;k<voxelize_res;k++)
            {
                if(voxel_model_[i][j][k]==1)
                {
                    Vec3f voxel;
                    voxel[0]=g_model_length/(float)voxelize_res*(i+0.5)-0.5*g_model_length;
                    voxel[1]=g_model_length/(float)voxelize_res*(j+0.5)-0.5*g_model_length;
                    voxel[2]=g_model_length/(float)voxelize_res*(k+0.5)-0.5*g_model_length;
                    voxel_count++;
                    volume+=pow(gap_,3);
                    centroid+=voxel;
                    //compute inertia tensor
                    inertia_tensor_wo_density[0][0]+=(voxel[1]*voxel[1]+voxel[2]*voxel[2])*pow(gap_,3);
                    inertia_tensor_wo_density[1][1]+=(voxel[0]*voxel[0]+voxel[2]*voxel[2])*pow(gap_,3);
                    inertia_tensor_wo_density[2][2]+=(voxel[1]*voxel[1]+voxel[0]*voxel[0])*pow(gap_,3);
                    inertia_tensor_wo_density[0][1]+=-voxel[0]*voxel[1]*pow(gap_,3);
                    inertia_tensor_wo_density[1][2]+=-voxel[1]*voxel[2]*pow(gap_,3);
                    inertia_tensor_wo_density[2][0]+=-voxel[2]*voxel[0]*pow(gap_,3);
                }
            }
        }
    }

    centroid/=(float)voxel_count;
    inertia_tensor_wo_density[1][0]=inertia_tensor_wo_density[0][1];
    inertia_tensor_wo_density[2][1]=inertia_tensor_wo_density[1][2];
    inertia_tensor_wo_density[0][2]=inertia_tensor_wo_density[2][0];

    cout<<"volume: "<<volume<<endl;
    cout<<"centroid: "<<centroid<<endl;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            cout<<inertia_tensor_wo_density[i][j]<<"\t";
        }cout<<endl;
    }

    //save quantities
    string filename = g_project_path+"/"+g_model_folder+"/Quantities/"+g_model_name+".qu";
    FILE *fp = fopen(filename.c_str(),"w");
	if (fp == (FILE *)NULL) {
		printf("File I/O Error:  Cannot create file \n");
		return;
	}
    fprintf(fp,"volume %f\n",volume);
    fprintf(fp,"centroid %f %f %f\n",centroid[0],centroid[1],centroid[2]);
    fprintf(fp,"inertia_tensor_wo_density %f %f %f %f %f %f %f %f %f\n",
        inertia_tensor_wo_density[0][0],inertia_tensor_wo_density[0][1],inertia_tensor_wo_density[0][2],
        inertia_tensor_wo_density[1][0],inertia_tensor_wo_density[1][1],inertia_tensor_wo_density[1][2],
        inertia_tensor_wo_density[2][0],inertia_tensor_wo_density[2][1],inertia_tensor_wo_density[2][2]);
    

    voxel_model_.clear();
}