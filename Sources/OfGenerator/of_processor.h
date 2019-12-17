#ifndef OF_INPUT_GENERATOR_H
#define OF_INPUT_GENERATOR_H
#include "Input/global_variables.h"
#include "Viewer/Viewer.h"

class OFProcessor
{
    //this class is for generate correct OpenFoam input folder
    //needed parameters are provided in global variables
public:
    OFProcessor(){};
    void Init();
    void RunViewer();

    void LoadObject();
    void GenerateOneCase();
    void GenerateOneInput();

    //get normal-area list for all projections of a model
    void GenerateNormalList();

    void RunOneCase();

    void ReadMinEdgeLength();//for computing Courant Number
    void ReadCourantNumber();//for reading Courant Number in short time
    void ReadOneOFOutput();

    void AdjustTemplate();
    void AdjustTemplatecontrolDict();
    void AdjustTemplateInitialConditions();//contains velocity, k-omega parameters
    void AdjustTemplateblockMeshDict();//model size
    void AdjustTemplatesnappyHexMeshDict();//model size
    void AdjustTemplateForces();//density, axis, rotation center
    void AdjustTemplateForcesCoeffs();//velocity, axis, rotation center, area, length

    void ChangeGlobalVariables();
    void ComputeAllCFD();
    void ComputeCasesForOneModel();//loops of rotation, size, reynolds

    void NormalizeModel();
    void RotateModel();
    void GetMinLength();

private:
    Mesh rotated_model_;
    bool successful_run_;
};

#endif