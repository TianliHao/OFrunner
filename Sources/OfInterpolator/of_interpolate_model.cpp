#include "OfInterpolator/of_interpolator.h"
using namespace std;

//data format-index
//rotation_angle, axis, length, area, velo, Re, Cd, Cl, Cs, coeffM, force, torque, avg_normal
//0               123   4       5     6     7   8   9   10  111213  141516 171819  202122

//returned value is the output part of test_datapair
void OFInterpolator::MinDistInterpolator(std::vector<double>& test_datapair)
{
    vector<double> min_value(input_size_,DBL_MAX);
    vector<double> max_value(input_size_,DBL_MIN);
    //get min-max range for every dimensions
    for(int i=0;i<dataset_.size();i++)
    {   for(int a=0;a<input_size_;a++)
        {
            if(min_value[a]>dataset_[i][a]) min_value[a]=dataset_[i][a];
            if(max_value[a]<dataset_[i][a]) max_value[a]=dataset_[i][a];
        }
    }
    double min_dist=DBL_MAX;
    int nearest_index=0;
    for(int i=0;i<dataset_.size();i++)
    {   
        double current_dist=0;
        for(int a=0;a<input_size_;a++)
        {
            current_dist+=pow((dataset_[i][a]-test_datapair[a])/(max_value[a]-min_value[a]),2);
        }
        current_dist=sqrt(current_dist);
        if(current_dist<min_dist)
        {
            min_dist=current_dist;
            nearest_index=i;
        }
    }
    //just copy the nearest output value
    cout<<"interpolated output: ";
    for(int a=input_size_;a<test_datapair.size();a++)
    {
        test_datapair[a]=dataset_[nearest_index][a];
        cout<<test_datapair[a]<<" ";
    }cout<<endl;

    cout<<"input range: ";
    for(int a=0;a<input_size_;a++)
    {
        cout<<"min:"<<min_value[a]<<"-max:"<<max_value[a]<<"-input:"<<test_datapair[a]<<endl;;
    }

}

void OFInterpolator::LinearInterpolator(vector<double>& test_lowd_data)
{
    if(!linear_computed_)
    {
        linear_computed_=true;
        linear_weights_.clear();    
        for(int i=input_size_;i<dataset_[0].size();i++)
        {
            Eigen::VectorXf weights_for_one_output(input_size_+1);

            Eigen::VectorXf n_outputs(dataset_.size());
            Eigen::MatrixXf interpolate_mat(dataset_.size(),input_size_+1);
            for(int j=0;j<dataset_.size();j++)
            {
                n_outputs[j]=dataset_[j][i];
                for(int k=0;k<input_size_;k++)
                {
                    interpolate_mat(j,k)=(float)dataset_[j][k];
                }
                interpolate_mat(j,input_size_)=1;
            }
            Eigen::MatrixXf mat_square_inverse=(interpolate_mat.transpose()*interpolate_mat).inverse();
            
            weights_for_one_output=mat_square_inverse*interpolate_mat.transpose()*n_outputs;
            linear_weights_.push_back(weights_for_one_output);
        }

    }
    for(int i=input_size_;i<dataset_[0].size();i++)
    {
        test_lowd_data[i]=0;
        for(int j=0;j<input_size_;j++)
        {
            cout<<test_lowd_data[i]<<" "<<linear_weights_[i-input_size_](j)<<" "<<test_lowd_data[j]<<endl;
            test_lowd_data[i]+=linear_weights_[i-input_size_][j]*test_lowd_data[j];
        }
        test_lowd_data[i]+=linear_weights_[i-input_size_][input_size_]*1.0;
    }
    
    vector<double> min_value(input_size_,DBL_MAX);
    vector<double> max_value(input_size_,DBL_MIN);
    //get min-max range for every dimensions
    for(int i=0;i<dataset_.size();i++)
    {   for(int a=0;a<input_size_;a++)
        {
            if(min_value[a]>dataset_[i][a]) min_value[a]=dataset_[i][a];
            if(max_value[a]<dataset_[i][a]) max_value[a]=dataset_[i][a];
        }
    }

    cout<<"interpolated output: ";
    for(int a=input_size_;a<test_lowd_data.size();a++)
    {
        cout<<test_lowd_data[a]<<" ";
    }cout<<endl;
    cout<<"input range: ";
    for(int a=0;a<input_size_;a++)
    {
        cout<<"min:"<<min_value[a]<<"-max:"<<max_value[a]<<"-input:"<<test_lowd_data[a]<<endl;;
    }
}