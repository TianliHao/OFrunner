#include "Input/global_variables.h"
#include <QuickSolver/quick_solver.h>
#include "controls.hpp"
using namespace glm;
using namespace std;


glm::mat4 ViewMatrix;
glm::mat4 ProjectionMatrix;

glm::mat4 getViewMatrix(){
	return ViewMatrix;
}
glm::mat4 getProjectionMatrix(){
	return ProjectionMatrix;
}

glm::mat4 computeModelMotion()
{
	g_quick_solver.ComputeMotion(g_render_deltaT);

	//rotation first
	glm::mat4 model_motion_matrix=g_quick_solver.initial_rotation_;
	
	//translate
	glm::mat4 translate_matrix(1);
	translate_matrix=glm::translate(translate_matrix, g_quick_solver.coord_);
	model_motion_matrix=translate_matrix*model_motion_matrix;

	return model_motion_matrix;
}

void GetOFOutput()
{
	g_quick_solver.GetOFOutput();
}

void computeMatricesFromInputs(){
	if(g_program_mode==1||g_program_mode==4)
	{
	glm::vec3 position = glm::vec3( -g_model_length,0 ,0 ); 
	static double lastTime = glfwGetTime();
	double currentTime = glfwGetTime();
	float deltaTime = float(currentTime - lastTime);
	g_render_deltaT=deltaTime;
	cout<<"deltaTime: "<<deltaTime<<endl;

	glm::vec3 direction(1,0,0);
	glm::vec3 up = glm::vec3(0,1,0);	
	ViewMatrix       = glm::lookAt(position,position+direction,up);
	//lookat matrix get new coord represent
	ProjectionMatrix = glm::ortho(-g_model_length,g_model_length,
		-g_model_length,g_model_length,0.0,2*g_model_length);

	lastTime = currentTime;
		return;
	}


	glm::vec3 position = glm::vec3( 0,0 ,2*g_model_length ); 
	static double lastTime = glfwGetTime();
	double currentTime = glfwGetTime();
	float deltaTime = float(currentTime - lastTime);
	g_render_deltaT=deltaTime;
	cout<<"deltaTime: "<<deltaTime<<endl;

	glm::vec3 direction(0,0,-1);
	glm::vec3 right = glm::vec3(1,0,0);	
	glm::vec3 up = glm::cross(right,direction);
	ViewMatrix       = glm::lookAt(position,position+direction,up);

	ProjectionMatrix = glm::ortho(-g_model_length,g_model_length,
		-g_model_length,g_model_length,0.0,2*g_model_length);

	if(g_program_mode!=1)
	{
		ViewMatrix       = glm::lookAt(position,position+direction,up);
		double rate=8;
		g_model_length*=rate;
		ProjectionMatrix = glm::ortho(-g_model_length,g_model_length,
			-1*g_model_length,g_model_length,-g_model_length,g_model_length);
		g_model_length/=rate;		
	}
	lastTime = currentTime;
}


void computeMatricesFromInputs2(){
	glm::vec3 position=g_quick_solver.coord_;
	glm::vec3 direction=g_quick_solver.velocity_;
	position+=(float)g_model_length*2.0f*glm::normalize(direction);
	direction=-glm::normalize(direction);

	glm::vec3 right = glm::vec3(1,0,0);
	glm::vec3 up = glm::cross( right, direction );
	ViewMatrix       = glm::lookAt(position,position+direction,up);

	ProjectionMatrix = glm::ortho(-g_model_length,g_model_length,
		-g_model_length,g_model_length,g_model_length,3*g_model_length);

}