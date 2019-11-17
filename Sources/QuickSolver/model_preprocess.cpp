#include <QuickSolver/quick_solver.h>
#include <OfInterpolator/of_interpolator.h>
using namespace OpenMesh;
using namespace std;
using namespace glm;

void QuickSolver::Init()
{
    density_=0.0001;
    mass_=0.01;      //kg    
    gravitional_acc_=9.8e-2;     //m/s^2
    velocity_=glm::vec3(0.1,0.6,0);    //m/s


    //initialize
    //rotation
    initial_rotation_=glm::mat4(1.0);
    initial_rotation_=glm::rotate(initial_rotation_,
        (float)g_rotation_angle,glm::vec3(g_rotation_axis[0],g_rotation_axis[1],g_rotation_axis[2]));
    angular_velocity_=glm::vec3(0.05,0,0);
    camera_right_=glm::vec3(1,0,0);
    force_=glm::vec3(0,0,0);
    torque_=glm::vec3(0,0,0);
    coord_=glm::vec3(0,0,0);

    ReadPhysicalQuantities();
}

void QuickSolver::GetOFOutput()
{
    double normal_para=sqrt(g_avg_normal[0]*g_avg_normal[0]+g_avg_normal[1]*g_avg_normal[1]);
    cout<<"g0"<<g_avg_normal[0]<<"g1"<<g_avg_normal[1]<<"normalpara"<<normal_para<<endl;
    double normal_ortho=g_avg_normal[2];
    //use interpolator to get f, alpha(angle between f and n), t, beta(angle between t and n)
    float alpha=0, beta=0, f_para=0, f_ortho=0, t_para=0, t_ortho=0;
    float theta=0;//theta is between (1,0) and (n0,n1) //all angle use radians

    vector<double> testdatapair={g_model_length, sqrt(dot(velocity_,velocity_)), normal_para, normal_ortho,0,0,0,0,0,0};
    //g_of_interpolator.MinDistInterpolator(testdatapair);
    g_of_interpolator.LinearInterpolator(testdatapair);
    f_para=(float)testdatapair[g_of_interpolator.input_size_+0];
    f_ortho=(float)testdatapair[g_of_interpolator.input_size_+1];
    alpha=(float)testdatapair[g_of_interpolator.input_size_+2];
    t_para=(float)testdatapair[g_of_interpolator.input_size_+3];
    t_ortho=(float)testdatapair[g_of_interpolator.input_size_+4];
    beta=(float)testdatapair[g_of_interpolator.input_size_+5];
    
    theta=(float)ComputeRadiansFrom2DVec(g_avg_normal[0],g_avg_normal[1]);
    alpha=theta+alpha;
    beta=theta+beta;
    glm::vec4 force_in_velocity_coord(f_para*cos(alpha),f_para*sin(alpha),f_ortho,0);
    glm::vec4 torque_in_velocity_coord(t_para*cos(alpha),t_para*sin(alpha),t_ortho,0);
    glm::vec4 force_in_world_coord=transpose(velocity_coord_)*force_in_velocity_coord;
    glm::vec4 torque_in_world_coord=transpose(velocity_coord_)*torque_in_velocity_coord;

    force_=glm::vec3(force_in_world_coord[0],force_in_world_coord[1],force_in_world_coord[2]);
    torque_=glm::vec3(torque_in_world_coord[0],torque_in_world_coord[1],torque_in_world_coord[2]);
}

void outputmat4(glm::mat4& input)
{
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout<<input[i][j]<<" ";
        }cout<<endl;
    }
}

void QuickSolver::ComputeMotion(float deltaT)
{
    glm::vec3 acc=-force_/mass_+glm::vec3(0,-gravitional_acc_,0);
    cout<<"acc: "<<acc[0]<<" "<<acc[1]<<" "<<acc[2]<<endl;
    coord_+=velocity_*deltaT+0.5f*acc*deltaT*deltaT;
    velocity_+=acc*deltaT;
    cout<<"velocity: "<<velocity_[0]<<" "<<velocity_[1]<<" "<<velocity_[2]<<endl;
    cout<<"coord: "<<coord_[0]<<" "<<coord_[1]<<" "<<coord_[2]<<endl;

    //rotation
    glm::vec3 up=cross(camera_right_,-velocity_);
    velocity_coord_=glm::lookAt(glm::vec3(0,0,0),-velocity_,cross(up,-velocity_));
    rotate_coord_=transpose(initial_rotation_);
    glm::mat4 rotation_in_velocity_coord=velocity_coord_*initial_rotation_;
    //compute inertia
    glm::mat4 velo_rotate_relative_rotation=rotate_coord_*transpose(velocity_coord_);
    //cout<<"OUTPUT velo rotate reletave:"<<endl;outputmat4(velo_rotate_relative_rotation);
    glm::mat4 inertia_tensor_4(inertia_tensor_);
    inertia_tensor_4=(transpose(velo_rotate_relative_rotation))*inertia_tensor_4*velo_rotate_relative_rotation;    
    //cout<<"OUTPUT inertia:"<<endl;outputmat4(inertia_tensor_4);
    glm::vec3 inertia(abs(inertia_tensor_4[0][0]),abs(inertia_tensor_4[1][1]),abs(inertia_tensor_4[2][2]));//to be updated, actually computed from tensor

    //compute omega
    glm::vec3 angular_acc(torque_[0]/inertia[0],torque_[1]/inertia[1],torque_[2]/inertia[2]);
    angular_velocity_+=angular_acc*deltaT;
    glm::vec3 omega_axis=normalize(angular_velocity_);
    float omega_value=length(angular_velocity_);
    //simulate roation only use theta=new_omega*t
    glm::mat4 local_rotate_in_velocity_coord(1);
    local_rotate_in_velocity_coord=glm::rotate(local_rotate_in_velocity_coord,omega_value*deltaT,omega_axis);
    glm::mat4 new_rotation_in_velocity_coord=local_rotate_in_velocity_coord*rotation_in_velocity_coord;
    glm::mat4 new_rotation_in_model_coord=transpose(velocity_coord_)*new_rotation_in_velocity_coord;
    initial_rotation_=new_rotation_in_model_coord;
    //cout<<"OUTPUT initial-rotation:"<<endl;outputmat4(initial_rotation_);
}

void QuickSolver::ReadPhysicalQuantities()
{
    string filename=g_project_path+"/"+g_model_folder+"/Quantities/"+g_model_name+".qu";
    ifstream fin(filename);
    	if (fin.fail() == true)
	{
		cout << "Cannot load PhysicalQuantities file." << endl;
		return;
	}
    string param;
    while(fin>>param)
    {
        if(param=="volume")
            fin>>volume_;
        else if(param=="centroid")
            fin>>centroid_[0]>>centroid_[1]>>centroid_[2];
        else if(param=="inertia_tensor_wo_density")
        {
            fin>>inertia_tensor_wo_density_[0][0]>>inertia_tensor_wo_density_[0][1]>>inertia_tensor_wo_density_[0][2]
                >>inertia_tensor_wo_density_[1][0]>>inertia_tensor_wo_density_[1][1]>>inertia_tensor_wo_density_[1][2]
                >>inertia_tensor_wo_density_[2][0]>>inertia_tensor_wo_density_[2][1]>>inertia_tensor_wo_density_[2][2];
            inertia_tensor_=density_*inertia_tensor_wo_density_;
        }
    }
}