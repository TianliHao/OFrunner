#ifndef OF_INTERPOLATOR_H
#define OF_INTERPOLATOR_H

#include <Input/global_variables.h>
#include <OfInterpolator/mathfuncs.h>
double ComputeRadiansFrom2DVec(double x, double y);

class OFProcessor;
void AlphaBetaRadiansFromRotMat(Eigen::Matrix3d rot, double &alpha_radians, double &beta_radians);
void RotMatFromAlphaBeta(Eigen::Matrix3d &rot, double alpha_radians, double beta_radians);
void ViewXYZFromAlphaBeta(double &x, double &y, double &z, double alpha_radians, double beta_radians);
double SmallSphereAreaPartFromBeta(double beta_radians, double delta_alpha_radians, double delta_beta_radians);

class OFInterpolator
{
public:
    OFInterpolator(){};
    void Init();

//read data pairs
    void ReadInputPair();
    void ReadNormalList();
    void NormalMapToVector();
    
//get normal-force dataset
    void ChangeDataFormat();
    void ResetDatasetNormalArea();
    void AugmentNormalData();//by rotation of flow-direction-axis
    void InterpolateForNormalMap();
    
//compute the arrow visualization for force and torque
    void GetDatasetForOneModel();
    void GenerateOneModelList();
    void OutputArrowObj(int mode);
    void OutputCsvFile(int mode);
    double max_force_;
    double max_torque_;

    void OutputSHParameters();
    void OutputSHParameters_interpolated();//interpolated sample around the sphere
    void TestStability();
    std::vector<double> SH_coeff_F_;
    std::vector<double> SH_coeff_T_;

//interpolators
    void MinDistInterpolator(std::vector<double>& test_lowd_data);
    void LinearInterpolator(std::vector<double>& test_lowd_data);
    void ComputeLeastSquareParams();
    void LeastSquareFitting(std::vector<double>& input_data);
    
//get physical quantities
    void Voxelize();
    void voxel_to_obj();
    void save_node_info(FILE* meshfp,int a,int b,int c,int new_index,float x,float y,float z);
    void ComputePhysicalQuantities();

//output dataset to file
    void WriteVectorPair();
    int data_line_count_;
    int data_element_count_;


//datapair
    std::vector<std::vector<double>> dataset_;
    int input_size_;
    int output_size_;

//interpolator model parameters
    bool linear_computed_;
    std::vector<Eigen::VectorXf> linear_weights_;

//voxel data
    std::vector<std::vector<std::vector<int>>> voxel_model_;
    int voxel_res_;
    float gap_;//size of voxel

    inline void OutputDataset()
    {
        for(int i=0;i<dataset_.size();i++)
        {
            for(int j=0;j<dataset_[i].size();j++)
            {
                std::cout<<dataset_[i][j]<<" ";
            }std::cout<<std::endl;
        }
    }

};
#endif