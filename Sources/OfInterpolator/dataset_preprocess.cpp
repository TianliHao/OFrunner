#include "OfInterpolator/of_interpolator.h"
using namespace std;
using namespace OpenMesh;

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
        //GenerateOneModelList();
        OutputArrowObj();
        OutputCsvFile();
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

void OFInterpolator::ChangeDataFormat()
{
    dataset_.clear();
    for(int i=0;i<g_dataset.size();i++)
    {
        float normal_para, normal_ortho, alpha, f_para, f_ortho, beta, t_para, t_ortho, length, velocity;
        //all angle here is in radians
        normal_para=sqrt(g_dataset[i][20]*g_dataset[i][20]+g_dataset[i][21]*g_dataset[i][21]);//from camera coord
        normal_ortho=g_dataset[i][22];
        f_para=sqrt(g_dataset[i][15]*g_dataset[i][15]+g_dataset[i][16]*g_dataset[i][16]);//from world coord, into camera coord
        f_ortho=g_dataset[i][14];
        t_para=sqrt(g_dataset[i][18]*g_dataset[i][18]+g_dataset[i][19]*g_dataset[i][19]);
        t_ortho=g_dataset[i][17];
        double angle_normal=ComputeRadiansFrom2DVec(g_dataset[i][20],g_dataset[i][21]);
        double angle_force=ComputeRadiansFrom2DVec(g_dataset[i][16],g_dataset[i][15]);
        double angle_torque=ComputeRadiansFrom2DVec(g_dataset[i][19],g_dataset[i][18]);
        alpha=angle_force-angle_normal;//angle from n to force
        beta=angle_torque-angle_normal;//angle from n to torque
        length=g_dataset[i][4];
        velocity=g_dataset[i][6];
        std::vector<double> datapair;
        input_size_=4;
        datapair.push_back(length);//input
        datapair.push_back(velocity);
        datapair.push_back(normal_para);
        datapair.push_back(normal_ortho);
        datapair.push_back(f_para);//output
        datapair.push_back(f_ortho);
        datapair.push_back(alpha);
        datapair.push_back(t_para);
        datapair.push_back(t_ortho);
        datapair.push_back(beta);
        
        dataset_.push_back(datapair);
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