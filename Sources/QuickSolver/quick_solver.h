#ifndef OF_QUICKSOLVER_H
#define OF_QUICKSOLVER_H
#include "Input/global_variables.h"

void outputmat4(glm::mat4& input);
class QuickSolver
{
public:
    QuickSolver(){};
    void Init();
    void ModelPreprocess();
    void GetOFOutput();
    void ComputeMotion(float deltaT);

    void ReadPhysicalQuantities();

    float refresh_timestep_;

    //model
    float density_;
    float volume_;
    glm::vec3 centroid_;

    float mass_;
    glm::vec3 velocity_;
    glm::vec3 coord_;
    glm::vec3 force_;
    float gravitional_acc_;//perhaps ignore first

    //rotation
    glm::vec3 torque_;
    glm::mat3 inertia_tensor_wo_density_;
    glm::mat3 inertia_tensor_;
    glm::vec3 angular_velocity_;

    glm::vec3 camera_right_;
    glm::mat4 initial_rotation_;
    glm::mat4 rotate_coord_;//only compute inertia, compute angle between velocity_coord and this
    glm::mat4 velocity_coord_;
};


#endif