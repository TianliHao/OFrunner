#include "OfGenerator/of_processor.h"
using namespace std;

void OFProcessor::AdjustTemplate()
{
    AdjustTemplatecontrolDict();
    AdjustTemplateInitialConditions();
    AdjustTemplateblockMeshDict();
    AdjustTemplatesnappyHexMeshDict();
    AdjustTemplateForces();
    AdjustTemplateForcesCoeffs();
}

void OFProcessor::AdjustTemplatecontrolDict()
{
    string filename=g_project_path+"/"+g_case_path+"/system/controlDict";
    string tempfilename=filename+"_tmp";
    string line;
    ifstream fin(filename);
    ofstream fout(tempfilename);
    
    while(!fin.eof())
    {
        getline(fin,line);
        string headword;
        stringstream line_stream;
        line_stream<<line;
        line_stream>>headword;
        if(headword=="deltaT")
            fout<<"deltaT          "<<g_time_step<<";"<<endl;
        else if(headword=="endTime")
            fout<<"endTime         "<<g_time_step*g_total_step<<";"<<endl;
        else
            fout<<line<<endl;
    }
    rename(tempfilename.c_str(),filename.c_str());
}

void OFProcessor::AdjustTemplateInitialConditions()
{
    string filename=g_project_path+"/"+g_case_path+"/0/include/initialConditions";
    string tempfilename=filename+"_tmp";
    string line;
    string headword;
    ifstream fin(filename);
    ofstream fout(tempfilename);
    double ke=1.5*pow(g_velocity*0.16*pow(g_Re_number,-0.125),2);
    double omega=pow(ke,0.5)/(double)(0.038*g_model_length);
    cout<<"\n\nTurbulence Model Coeffs:\n"<<endl;
    cout<<"K=\t\t"<<ke<<endl;
    cout<<"omega=\t"<<omega<<endl;
    while(!fin.eof())
    {
        getline(fin,line);
        stringstream line_stream;
        line_stream<<line;
        line_stream>>headword;
        if(headword=="flowVelocity")
            fout<<"flowVelocity         ("<<g_velocity<<" 0 0);"<<endl;
        else if(headword=="turbulentKE")
            fout<<"turbulentKE          "<<ke<<";"<<endl;
        else if(headword=="turbulentOmega")
            fout<<"turbulentOmega       "<<omega<<";"<<endl;
        else
            fout<<line<<endl;
    }
    rename(tempfilename.c_str(),filename.c_str());
}

void OFProcessor::AdjustTemplateblockMeshDict()
{
    string filename=g_project_path+"/"+g_case_path+"/system/blockMeshDict";
    string tempfilename=filename+"_tmp";
    string line;
    string headword;
    ifstream fin(filename);
    ofstream fout(tempfilename);
    while(!fin.eof())
    {
        getline(fin,line);
        stringstream line_stream;
        line_stream<<line;
        line_stream>>headword;
        if(headword=="convertToMeters")
            fout<<"convertToMeters "<<g_model_length<<";"<<endl;
        else
            fout<<line<<endl;
    }
    rename(tempfilename.c_str(),filename.c_str());
}

void OFProcessor::AdjustTemplatesnappyHexMeshDict()
{
    double normalized_refine_box_min[3]={-1,-1.2,-1.2};
    double normalized_refine_box_max[3]={4,1.2,1.2};
    double normalized_location_in_mesh[3]={0,2,2};

    string filename=g_project_path+"/"+g_case_path+"/system/snappyHexMeshDict";
    string tempfilename=filename+"_tmp";
    ifstream fin(filename);
    ofstream fout(tempfilename);
    while(!fin.eof())
    {
        string line;
        string headword;
        getline(fin,line);
        stringstream line_stream;
        line_stream<<line;
        line_stream>>headword;
        if(headword=="min")
            fout<<"        min ("
            <<normalized_refine_box_min[0]*g_model_length<<" "
            <<normalized_refine_box_min[1]*g_model_length<<" "
            <<normalized_refine_box_min[2]*g_model_length<<");"<<endl;
        else if(headword=="max")
            fout<<"        max ("
            <<normalized_refine_box_max[0]*g_model_length<<" "
            <<normalized_refine_box_max[1]*g_model_length<<" "
            <<normalized_refine_box_max[2]*g_model_length<<");"<<endl;
        else if(headword=="locationInMesh")
            fout<<"    locationInMesh ("
            <<normalized_location_in_mesh[0]*g_model_length<<" "
            <<normalized_location_in_mesh[1]*g_model_length<<" "
            <<normalized_location_in_mesh[2]*g_model_length<<");"<<endl;
        else
            fout<<line<<endl;
    }
    rename(tempfilename.c_str(),filename.c_str());
}


void OFProcessor::AdjustTemplateForces()
{
    string filename=g_project_path+"/"+g_case_path+"/system/forces";
    string tempfilename=filename+"_tmp";
    string line;
    string headword;
    ifstream fin(filename);
    ofstream fout(tempfilename);
    while(!fin.eof())
    {
        getline(fin,line);
        stringstream line_stream;
        line_stream<<line;
        line_stream>>headword;
        if(headword=="CofR")
            fout<<"    CofR          (0 0 0);"<<endl;
        else
            fout<<line<<endl;
    }
    rename(tempfilename.c_str(),filename.c_str());
}

void OFProcessor::AdjustTemplateForcesCoeffs()
{
    string filename=g_project_path+"/"+g_case_path+"/system/forceCoeffs";
    string tempfilename=filename+"_tmp";
    string line;
    string headword;
    ifstream fin(filename);
    ofstream fout(tempfilename);
    while(!fin.eof())
    {
        getline(fin,line);
        stringstream line_stream;
        line_stream<<line;
        line_stream>>headword;
        if(headword=="CofR")
            fout<<"    CofR            (0 0 0);"<<endl;
        else if(headword=="magUInf")
            fout<<"    magUInf         "<<g_velocity<<";"<<endl;
        else if(headword=="lRef")
            fout<<"    lRef            "<<g_diagonal_length<<";"<<endl;
        else if(headword=="Aref")
            fout<<"    Aref            "<<g_projection_area<<";"<<endl;
        else
            fout<<line<<endl;
    }
    rename(tempfilename.c_str(),filename.c_str());
}