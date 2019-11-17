#include "Viewer/Viewer.h"

GLFWwindow* window;
GLFWwindow* window2;
using namespace std;

Viewer::Viewer()
{
    RenderCase();
}


void PPMWriter(unsigned char *in,char *name,int dimx, int dimy)
{
    int i, j;
    FILE *fp = fopen(name, "wb"); /* b - binary mode */
    FILE * saveTxt=fopen((std::string(name)+".tlh").c_str(),"w");
    (void) fprintf(fp, "P6\t%d %d\t255\t", dimx, dimy);
    //CCCCCCCCCCCCCCCompute the global descriptor
    Eigen::Vector3d avg_normal(0,0,0);
    int pixel_count=0;

    for (i = dimx-1; i >= 0; --i)
    {
        for (j = 0; j < dimy; ++j)
        {
        static unsigned char color[3];
        color[0] = in[3*j+3*i*dimx];  /* red */
        color[1] = in[3*j+3*i*dimx+1];  /* green */
        color[2] = in[3*j+3*i*dimx+2];  /* blue */
        fprintf(saveTxt,"%d  %d  %d\n",color[0],color[1],color[2]);
        (void) fwrite(color, 1, 3, fp);

        //CCCCCCCCCCCCCCCompute the global descriptor
        if(color[0]==0&&color[1]==0&color[2]==0)
            continue;
		//std::cout<<(int)color[0]<<" "<<(int)color[1]<<" "<<(int)color[2]<<std::endl;
        Eigen::Vector3d this_normal=Eigen::Vector3d(((double)((int)color[0]+1))/128.0-1,((double)((int)color[1]+1))/(double)128.0-1,
        ((double)((int)color[2]+1))/(double)128.0-1);
		avg_normal+=this_normal;
		//std::cout<<this_normal<<std::endl<<std::endl;
        pixel_count++;
        }
    }
	g_avg_normal=avg_normal/(double)pixel_count;
	g_projection_area=pixel_count/(double)(g_render_resolution*g_render_resolution)
		*(2*g_model_length)*(2*g_model_length);
    std::cout<<"average normal: "<<g_avg_normal<<std::endl;
    std::cout<<"colored pixel count: "<<pixel_count<<std::endl;
	std::cout<<"projection area: "<<g_projection_area<<std::endl;

    fclose(saveTxt);
    (void) fclose(fp);

}

void saveImage(int index)
{
    unsigned char* image = (unsigned char*)malloc(sizeof(unsigned char) * 3
		* g_render_resolution * g_render_resolution);
    glReadPixels(0, 0, g_render_resolution, g_render_resolution, GL_RGB, GL_UNSIGNED_BYTE, image);
    // Warning : enregistre de bas en haut
	if(index==1)
	{
		char buffer [200];
		sprintf(buffer,"../Render/%s_%.2lf_%.2lf_%.2lf_%.2lf.ppm",g_model_name.c_str(),g_rotation_angle,
			g_rotation_axis(0),g_rotation_axis(1),g_rotation_axis(2));

		PPMWriter(image,buffer,g_render_resolution, g_render_resolution);
	}
	if(index==2)
	{
		char buffer [200];
		sprintf(buffer,"../Render/temp.ppm");
		PPMWriter(image,buffer,g_render_resolution, g_render_resolution);
	}
}



int Viewer::RenderCase()
{
    int window_width=g_render_resolution;
    int window_height=g_render_resolution;
	// Initialise GLFW
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow( window_width, window_height, "Fixed View", NULL, NULL);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwSetWindowPos(window, 0, 0);
	// Open a window and create its OpenGL context
	glfwMakeContextCurrent(window);
		
	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}
	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark background
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS); 

	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	std::string shader_path=g_project_path+"/Sources/Viewer";
	GLuint programID = LoadShaders((shader_path+"/NormalMap.vertexshader").c_str(),(shader_path+"/NormalMap.fragmentshader").c_str());

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");
	GLuint ViewMatrixID = glGetUniformLocation(programID, "V");
	GLuint ModelMatrixID = glGetUniformLocation(programID, "M");
	
	// Read our .obj file
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;

    for(int i=0;i<g_rotated_model.n_faces();i++)
    {
        FaceHandle fh=g_rotated_model.face_handle(i);
        FaceVertexIter fv_it=g_rotated_model.fv_begin(fh);
        FaceVertexIter fv_it_end=g_rotated_model.fv_end(fh);
        for(; fv_it!=fv_it_end;fv_it++)
        {
            VertexHandle vh=fv_it.handle();
            OpenMesh::Vec3d v_coord=g_rotated_model.point(vh);
            OpenMesh::Vec3d v_normal=g_rotated_model.normal(vh);
            glm::vec3 render_vertex(v_coord[0],v_coord[1],v_coord[2]);
            glm::vec3 render_normal(v_normal[0],v_normal[1],v_normal[2]);
            vertices.push_back(render_vertex);
            normals.push_back(render_normal);
        }
    }
	// Load it into a VBO

	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);

	GLuint normalbuffer;
	glGenBuffers(1, &normalbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
	glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);

//////////////start init window2
	// Open a window and create its OpenGL context
	window2 = glfwCreateWindow( window_width, window_height, "Projection Area", NULL, NULL);
	if( window2 == NULL ){
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		getchar();
		glfwTerminate();
		return -1;
	}
	// Open a window and create its OpenGL context
	glfwMakeContextCurrent(window2);
	glfwSetWindowPos(window2, 1000, 500);
		
	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window2, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark background
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS); 

	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	GLuint VertexArrayID2;
	glGenVertexArrays(1, &VertexArrayID2);
	glBindVertexArray(VertexArrayID2);

	// Create and compile our GLSL program from the shaders
	GLuint programID2 = LoadShaders((shader_path+"/NormalMap.vertexshader").c_str(),(shader_path+"/NormalMap.fragmentshader").c_str());

	// Get a handle for our "MVP" uniform
	GLuint MatrixID2 = glGetUniformLocation(programID2, "MVP");
	GLuint ViewMatrixID2 = glGetUniformLocation(programID2, "V");
	GLuint ModelMatrixID2 = glGetUniformLocation(programID2, "M");
	
	// Load it into a VBO
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);

	glGenBuffers(1, &normalbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
	glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);





	do{
		glfwMakeContextCurrent(window);
		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Use our shader
		glUseProgram(programID);

		//compute model motion
		glm::mat4 ModelMatrix(1.0);
		if(g_program_mode!=1)
			ModelMatrix=computeModelMotion();
		// Compute the MVP matrix from keyboard and mouse input
		computeMatricesFromInputs();
		glm::mat4 ProjectionMatrix = getProjectionMatrix();
		glm::mat4 ViewMatrix = getViewMatrix();
		glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(
			0,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 3rd attribute buffer : normals
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
		glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,(void*)0);
		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, vertices.size() );
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);

		// Swap buffers
		glfwSwapBuffers(window);
		
		if(g_program_mode==1)
		{
    		saveImage(g_program_mode); break;
		}

/////////////////////loop of window2
		glfwMakeContextCurrent(window2);

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// Use our shader
		glUseProgram(programID2);

		// Compute the MVP matrix from keyboard and mouse input
		computeMatricesFromInputs2();
		glm::mat4 ProjectionMatrix2 = getProjectionMatrix();
		glm::mat4 ViewMatrix2 = getViewMatrix();
		glm::mat4 MVP2 = ProjectionMatrix2 * ViewMatrix2 * ModelMatrix;


		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID2, 1, GL_FALSE, &MVP2[0][0]);
		glUniformMatrix4fv(ModelMatrixID2, 1, GL_FALSE, &ModelMatrix[0][0]);
		glUniformMatrix4fv(ViewMatrixID2, 1, GL_FALSE, &ViewMatrix2[0][0]);

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,(void*)0);

		// 3rd attribute buffer : normals
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, normalbuffer);
		glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,0,(void*)0);
		// Draw the triangles !
		glDrawArrays(GL_TRIANGLES, 0, vertices.size() );
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);

		// Swap buffers
		glfwSwapBuffers(window2);

		if(g_program_mode==2)
    		saveImage(g_program_mode);

		glfwPollEvents();
		
		GetOFOutput();
	} // Check if the ESC key was pressed or the window was closed
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 &&
		   glfwGetKey(window2, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window2) == 0 );

	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1, &normalbuffer);
	glDeleteProgram(programID);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteProgram(programID2);
	glDeleteVertexArrays(1, &VertexArrayID2);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}