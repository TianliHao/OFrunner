#include "OfInterpolator/of_interpolator.h"
using namespace std;

void OFInterpolator::voxel_to_obj()
{
	float x,y,z;
    float gap=gap_;
	x=y=z=-0.5*g_model_length-gap/2.0;

    vector<vector<vector<int>>> cube_index=voxel_model_;
    vector<vector<vector<int>>> node_index(voxel_res_+1,vector<vector<int>>(voxel_res_+1,vector<int>(voxel_res_+1,0)));

    string filename=g_project_path+"/"+g_model_folder+"/Voxels/"+g_model_name;
	//char filename[100];	//sprintf(filename, (filename1+"%d.obj").c_str(), nameindex);
	FILE *meshfp = fopen(filename.c_str(),"w");
	if (meshfp == (FILE *)NULL) {
		printf("File I/O Error:  Cannot create file \n");
		return;
	}
	int CUBE_SIZE=voxel_res_;
	int new_index=1;
	for(int a=0;a<CUBE_SIZE;a++)
	{
		for(int b=0;b<CUBE_SIZE;b++)
		{
			for(int c=0;c<CUBE_SIZE;c++)
			{
				if(cube_index[a][b][c]==1)
				{
					float coor_x;
					float coor_y;
					float coor_z;
					if(a==0||cube_index[a-1][b][c]==0)
					{	

						node_index[a][b][c]=new_index++;
						save_node_info(meshfp,a,b,c,new_index-1,x,y,z);
						node_index[a][b][c+1]=new_index++;
						save_node_info(meshfp,a,b,c+1,new_index-1,x,y,z);
						node_index[a][b+1][c+1]=new_index++;
						save_node_info(meshfp,a,b+1,c+1,new_index-1,x,y,z);
						node_index[a][b+1][c]=new_index++;
						save_node_info(meshfp,a,b+1,c,new_index-1,x,y,z);
					}
					if(a==CUBE_SIZE-1||cube_index[a+1][b][c]==0)
					{	
						node_index[a+1][b][c]=new_index++;
						save_node_info(meshfp,a+1,b,c,new_index-1,x,y,z);
						node_index[a+1][b+1][c]=new_index++;
						save_node_info(meshfp,a+1,b+1,c,new_index-1,x,y,z);
						node_index[a+1][b+1][c+1]=new_index++;
						save_node_info(meshfp,a+1,b+1,c+1,new_index-1,x,y,z);
						node_index[a+1][b][c+1]=new_index++;
						save_node_info(meshfp,a+1,b,c+1,new_index-1,x,y,z);
					}
					if(b==0||cube_index[a][b-1][c]==0)
					{
						node_index[a][b][c]=new_index++;
						save_node_info(meshfp,a,b,c,new_index-1,x,y,z);
						node_index[a+1][b][c]=new_index++;
						save_node_info(meshfp,a+1,b,c,new_index-1,x,y,z);
						node_index[a+1][b][c+1]=new_index++;
						save_node_info(meshfp,a+1,b,c+1,new_index-1,x,y,z);
						node_index[a][b][c+1]=new_index++;
						save_node_info(meshfp,a,b,c+1,new_index-1,x,y,z);
					}
					if(b==CUBE_SIZE-1||cube_index[a][b+1][c]==0)
					{
						node_index[a][b+1][c]=new_index++;
						save_node_info(meshfp,a,b+1,c,new_index-1,x,y,z);
						node_index[a+1][b+1][c]=new_index++;
						save_node_info(meshfp,a+1,b+1,c,new_index-1,x,y,z);
						node_index[a+1][b+1][c+1]=new_index++;
						save_node_info(meshfp,a+1,b+1,c+1,new_index-1,x,y,z);
						node_index[a][b+1][c+1]=new_index++;
						save_node_info(meshfp,a,b+1,c+1,new_index-1,x,y,z);
					}
					if(c==0||cube_index[a][b][c-1]==0)
					{
						node_index[a][b][c]=new_index++;
						save_node_info(meshfp,a,b,c,new_index-1,x,y,z);
						node_index[a+1][b][c]=new_index++;
						save_node_info(meshfp,a+1,b,c,new_index-1,x,y,z);
						node_index[a+1][b+1][c]=new_index++;
						save_node_info(meshfp,a+1,b+1,c,new_index-1,x,y,z);
						node_index[a][b+1][c]=new_index++;
						save_node_info(meshfp,a,b+1,c,new_index-1,x,y,z);
					}
					if(c==CUBE_SIZE-1||cube_index[a][b][c+1]==0)
					{
						node_index[a][b][c+1]=new_index++;
						save_node_info(meshfp,a,b,c+1,new_index-1,x,y,z);
						node_index[a+1][b][c+1]=new_index++;
						save_node_info(meshfp,a+1,b,c+1,new_index-1,x,y,z);
						node_index[a+1][b+1][c+1]=new_index++;
						save_node_info(meshfp,a+1,b+1,c+1,new_index-1,x,y,z);
						node_index[a][b+1][c+1]=new_index++;
						save_node_info(meshfp,a,b+1,c+1,new_index-1,x,y,z);
					}
				}
			}
		}
	}

	for(int a=0;a<CUBE_SIZE;a++)
	{
		for(int b=0;b<CUBE_SIZE;b++)
		{
			for(int c=0;c<CUBE_SIZE;c++)
			{
				if(cube_index[a][b][c]==1)
				{
					int index1,index2,index3;
					if(a==0||cube_index[a-1][b][c]==0)
					{	
						index1=node_index[a][b+1][c+1];
						index2=node_index[a][b][c];
						index3=node_index[a][b][c+1];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						index1=node_index[a][b+1][c+1];
						index2=node_index[a][b+1][c];
						index3=node_index[a][b][c];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						//[a][b+1][c+1] [a][b][c+1] [a][b][c]
						//[a][b+1][c+1] [a][b][c] [a][b+1][c]
					}
					if(a==CUBE_SIZE-1||cube_index[a+1][b][c]==0)
					{	
						index1=node_index[a+1][b+1][c+1];
						index2=node_index[a+1][b][c+1];
						index3=node_index[a+1][b][c];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						index1=node_index[a+1][b+1][c+1];
						index2=node_index[a+1][b][c];
						index3=node_index[a+1][b+1][c];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						//[a+1][b+1][c+1] [a+1][b][c] [a+1][b][c+1]
						//[a+1][b+1][c+1] [a+1][b+1][c] [a+1][b][c]
					}
					if(b==0||cube_index[a][b-1][c]==0)
					{
						index1=node_index[a+1][b][c];
						index2=node_index[a][b][c+1];
						index3=node_index[a][b][c];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						index1=node_index[a+1][b][c];
						index2=node_index[a+1][b][c+1];
						index3=node_index[a][b][c+1];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						//[a+1][b][c] [a][b][c] [a][b][c+1]
						//[a+1][b][c] [a][b][c+1] [a+1][b][c+1]
					}
					if(b==CUBE_SIZE-1||cube_index[a][b+1][c]==0)
					{
						index1=node_index[a+1][b+1][c];
						index2=node_index[a][b+1][c+1];
						index3=node_index[a+1][b+1][c+1];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						index1=node_index[a+1][b+1][c];
						index2=node_index[a][b+1][c];
						index3=node_index[a][b+1][c+1];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						//[a+1][b+1][c] [a+1][b+1][c+1] [a][b+1][c+1]
						//[a+1][b+1][c] [a][b+1][c+1] [a][b+1][c]
					}
					if(c==0||cube_index[a][b][c-1]==0)
					{
						index1=node_index[a+1][b][c];
						index2=node_index[a][b+1][c];
						index3=node_index[a+1][b+1][c];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						index1=node_index[a+1][b][c];
						index2=node_index[a][b][c];
						index3=node_index[a][b+1][c];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						//[a+1][b][c] [a+1][b+1][c] [a][b+1][c]
						//[a+1][b][c] [a][b+1][c] [a][b][c]
					}
					if(c==CUBE_SIZE-1||cube_index[a][b][c+1]==0)
					{
						index1=node_index[a+1][b][c+1];
						index2=node_index[a][b+1][c+1];
						index3=node_index[a][b][c+1];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						index1=node_index[a+1][b][c+1];
						index2=node_index[a+1][b+1][c+1];
						index3=node_index[a][b+1][c+1];
						fprintf(meshfp,"f %d %d %d\n",index1,index2,index3);
						//[a+1][b][c+1] [a][b][c+1] [a][b+1][c+1]
						//[a+1][b][c+1] [a][b+1][c+1] [a+1][b+1][c+1]
					}
				}
			}
		}
	}
	fclose(meshfp);
}

void OFInterpolator::save_node_info(FILE* meshfp,int a,int b,int c,int new_index,float x,float y,float z)
{
    
	float coor_x=gap_*a+(float)x;
	float coor_y=gap_*b+(float)y;
	float coor_z=gap_*c+(float)z;
	fprintf(meshfp,"v %.5f %.5f %.5f\n",coor_x,coor_y,coor_z);
}