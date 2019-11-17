#ifndef OF_INTERPOLATOR_H
#define OF_INTERPOLATOR_H

#include <Input/global_variables.h>
double ComputeRadiansFrom2DVec(double x, double y);
class OFInterpolator
{
public:
    OFInterpolator(){};
    void Init();

//read data pairs
    void ReadInputPair();
    void ReadNormalMap();
    void NormalMapToVector();

//get real dataset //lower dimension than input-23-numbers
    void ChangeDataFormat();

//compute the arrow visualization for force and torque
    void GetDatasetForOneModel();
    void GenerateOneModelList();
    void OutputArrowObj();
    void OutputCsvFile();
    double max_force_;
    double max_torque_;

//interpolators
    void MinDistInterpolator(std::vector<double>& test_lowd_data);
    void LinearInterpolator(std::vector<double>& test_lowd_data);
    void DNNInterpolator();

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