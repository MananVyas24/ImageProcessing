/**
EE569 HW #3 Problem 2:Perspective Transformation & Imaging Geometry
    Name:	Manan Vyas
	USC ID:	7483-8632-00
	email:	mvyas@usc.edu
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "IO.h"
#include "Filters.h"
#include "ImageProcessingFunctions.h"
#include "opencv2/core/core.hpp" //use opencv to calculate inverse matrix in B(b)
#include <math.h>

#define PI 3.1415926

using namespace cv;
using namespace std;

// Part (A):Pre-processing
bool getRGBfromLocation(RGB* color, unsigned char* src1,unsigned char* src2,unsigned char* src3,unsigned char* src4,unsigned char* src5,int Size,float* loc)
{
    color->R = 0;
    color->G = 0;
    color->B = 0;
    float x = loc[0];
    float y = loc[1];
    float z = loc[2];
    unsigned char* src = NULL;
    if(z==1 && x>=-1 && x<=1 && y>=-1 && y<=1)
    {
        src = src1;
        int X = (x+1)/2*(Size-1);
        int Y = (y+1)/2*(Size-1);
        color->R = src[3*(X*Size+Y)];
        color->G = src[3*(X*Size+Y)+1];
        color->B = src[3*(X*Size+Y)+2];
        return true;
    }
    if(x==1 && z>=-1 && z<=1 && y>=-1 && y<=1)
    {
        src = src2;
        int X = (1-z)/2*(Size-1);
        int Y = (y+1)/2*(Size-1);
        color->R = src[3*(X*Size+Y)];
        color->G = src[3*(X*Size+Y)+1];
        color->B = src[3*(X*Size+Y)+2];
        return true;
    }
    if(y==1 && x>=-1 && x<=1 && z>=-1 && z<=1)
    {
        src = src3;
        int Y = (1-x)/2*(Size-1);
        int X = (1-z)/2*(Size-1);
        color->R = src[3*(X*Size+Y)];
        color->G = src[3*(X*Size+Y)+1];
        color->B = src[3*(X*Size+Y)+2];
        return true;
    }
    if(x==-1 && z>=-1 && z<=1 && y>=-1 && y<=1)
    {
        src = src4;
        int X = (1-z)/2*(Size-1);
        int Y = (1-y)/2*(Size-1);
        color->R = src[3*(X*Size+Y)];
        color->G = src[3*(X*Size+Y)+1];
        color->B = src[3*(X*Size+Y)+2];
        return true;
    }
    if(y==-1 && x>=-1 && x<=1 && z>=-1 && z<=1)
    {
        src = src5;
        int X = (1-z)/2*(Size-1);
        int Y = (x+1)/2*(Size-1);
        color->R = src[3*(X*Size+Y)];
        color->G = src[3*(X*Size+Y)+1];
        color->B = src[3*(X*Size+Y)+2];
        return true;
    }
    return false;
}
void partA(unsigned char* src1,unsigned char* src2,unsigned char* src3,unsigned char* src4,unsigned char* src5,int Size)
{
    printf(" A begins...\n");

    RGB rgb;
    float xyz[3];
    cout<<"(X Y Z):";
    cin>>xyz[0]>>xyz[1]>>xyz[2];
    if(getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,xyz))
    {
        cout<<"(R G B):";
        cout<<(int)rgb.R<<" "<<(int)rgb.G<<" "<<(int)rgb.B<<endl;
    }
    else
    {
        cout<<"The location is not on any surface!\n";
    }
    printf(" A ends...\n");
}

//(B) Capturing 3D scene
float getPosOfCapture(int pDensity,int SizeCapture,float* T,float*inPos,int* outPos)
{
    //
    float out[3];
    MatrixMul(T,3,4,inPos,4,1,out);
    outPos[0] = out[0]/out[2]*pDensity;
    outPos[1] = out[1]/out[2]*pDensity;
    /*
                    Y=-pDensity
                     _______
        X=-pDensity |___|___|X=pDensity
                    |___|___|
                    Y=pDensity
    */
    int offset = SizeCapture/2.0;
    outPos[0] += offset;
    outPos[1] += offset;
    int tmp = outPos[0];
    outPos[0] = outPos[1];
    outPos[1] = tmp;
    /*
            X=0
             _______
        Y=0 |___|___|Y=2*pDensity
            |___|___|
            X=2*pDensity
    */

    return out[2]; // return Z
}
unsigned char* Capture(unsigned char* src1,unsigned char* src2,unsigned char* src3,unsigned char* src4,unsigned char* src5,int Size,int pDensity,int SizeCapture,float* T,float* cameraLoc,float* rM,float* depth)
{
    // depth is a returned depth-image
    // rM is rotation matrix
    //
    unsigned char* capture = (unsigned char*)malloc(sizeof(unsigned char)*SizeCapture*SizeCapture*3);
    int *checker = (int*)malloc(sizeof(int)*SizeCapture*SizeCapture);
    memset(checker,0,sizeof(int)*SizeCapture*SizeCapture);
    float x,y,z;
    RGB rgb;

    for(int i=0;i<SizeCapture*SizeCapture;i++)
    {
        depth[i] = -1; //initialize depth as -1
    }

    // 1
    z = 1;
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            x = 2.0*i/Size-1;
            y = 2.0*j/Size-1;
            float inPos_t[] = {x,y,z,1};
            float inPos[4];
            MatrixMul(rM,4,4,inPos_t,4,1,inPos);
            int outPos[2];
            float dist = getPosOfCapture(pDensity,SizeCapture,T,inPos,outPos);
            if(outPos[0]>=0 && outPos[0]<SizeCapture && outPos[1]>=0 && outPos[1]<SizeCapture)
            {
                if(checker[outPos[0]*SizeCapture+outPos[1]]==0)
                {
                    getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                    capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                    depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    checker[outPos[0]*SizeCapture+outPos[1]] = 1;
                }
                else
                {
                    if(dist<depth[outPos[0]*SizeCapture+outPos[1]])
                    {
                        getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                        capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                        depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    }
                }
            }
        }
    }

    // 2
    x = 1;
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            y = 2.0*i/Size-1;
            z = 1-2.0*j/Size;
            float inPos_t[] = {x,y,z,1};
            float inPos[4];
            MatrixMul(rM,4,4,inPos_t,4,1,inPos);
            int outPos[2];
            float dist = getPosOfCapture(pDensity,SizeCapture,T,inPos,outPos);
            if(outPos[0]>=0 && outPos[0]<SizeCapture && outPos[1]>=0 && outPos[1]<SizeCapture)
            {
                if(checker[outPos[0]*SizeCapture+outPos[1]]==0)
                {
                    getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                    capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                    depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    checker[outPos[0]*SizeCapture+outPos[1]] = 1;
                }
                else
                {
                    if(dist<depth[outPos[0]*SizeCapture+outPos[1]])
                    {
                        getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                        capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                        depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    }
                }
            }
        }
    }
    // 3
    y = 1;
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            x = 1-2.0*i/Size;
            z = 1-2.0*j/Size;
            float inPos_t[] = {x,y,z,1};
            float inPos[4];
            MatrixMul(rM,4,4,inPos_t,4,1,inPos);
            int outPos[2];
            float dist = getPosOfCapture(pDensity,SizeCapture,T,inPos,outPos);
            if(outPos[0]>=0 && outPos[0]<SizeCapture && outPos[1]>=0 && outPos[1]<SizeCapture)
            {
                if(checker[outPos[0]*SizeCapture+outPos[1]]==0)
                {
                    getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                    capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                    depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    checker[outPos[0]*SizeCapture+outPos[1]] = 1;
                }
                else
                {
                    if(dist<depth[outPos[0]*SizeCapture+outPos[1]])
                    {
                        getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                        capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                        depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    }
                }
            }
        }
    }
    // 4
    x = -1;
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            y = 1-2.0*i/Size;
            z = 1-2.0*j/Size;
            float inPos_t[] = {x,y,z,1};
            float inPos[4];
            MatrixMul(rM,4,4,inPos_t,4,1,inPos);
            int outPos[2];
            float dist = getPosOfCapture(pDensity,SizeCapture,T,inPos,outPos);
            if(outPos[0]>=0 && outPos[0]<SizeCapture && outPos[1]>=0 && outPos[1]<SizeCapture)
            {
                if(checker[outPos[0]*SizeCapture+outPos[1]]==0)
                {
                    getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                    capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                    depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    checker[outPos[0]*SizeCapture+outPos[1]] = 1;
                }
                else
                {
                    if(dist<depth[outPos[0]*SizeCapture+outPos[1]])
                    {
                        getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                        capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                        depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    }
                }
            }
        }
    }
    // 5
    y = -1;
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            z = 1-2.0*i/Size;
            x = 2.0*j/Size-1;
            float inPos_t[] = {x,y,z,1};
            float inPos[4];
            MatrixMul(rM,4,4,inPos_t,4,1,inPos);
            int outPos[2];
            float dist = getPosOfCapture(pDensity,SizeCapture,T,inPos,outPos);
            if(outPos[0]>=0 && outPos[0]<SizeCapture && outPos[1]>=0 && outPos[1]<SizeCapture)
            {
                if(checker[outPos[0]*SizeCapture+outPos[1]]==0)
                {
                    getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                    capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                    capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                    depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    checker[outPos[0]*SizeCapture+outPos[1]] = 1;
                }
                else
                {
                    if(dist<depth[outPos[0]*SizeCapture+outPos[1]])
                    {
                        getRGBfromLocation(&rgb,src1,src2,src3,src4,src5,Size,inPos_t);
                        capture[3*(outPos[0]*SizeCapture+outPos[1])] = rgb.R;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+1] = rgb.G;
                        capture[3*(outPos[0]*SizeCapture+outPos[1])+2] = rgb.B;
                        depth[outPos[0]*SizeCapture+outPos[1]] = dist;
                    }
                }
            }
        }
    }

    delete checker;
    return capture;
}
void partB(unsigned char* src1,unsigned char* src2,unsigned char* src3,unsigned char* src4,unsigned char* src5,int Size,int pDensity,int SizeCapture)
{
    printf(" B begins...\n");

    float cameraLoc[3];
    float cameraDirX[3];
    float cameraDirY[3];
    float cameraDirZ[3];
    /*
    cout<<"Camera (X Y Z):";
    cin>>cameraLoc.X>>cameraLoc.Y>>cameraLoc.Z;
    */
    cameraLoc[0] = 5;
    cameraLoc[1] = 5;
    cameraLoc[2] = 5;

    cameraDirX[0] = -1.0/sqrt(2);
    cameraDirX[1] = 1.0/sqrt(2);
    cameraDirX[2] = 0;

    cameraDirY[0] = 1.0/sqrt(6);
    cameraDirY[1] = 1.0/sqrt(6);
    cameraDirY[2] = -2.0/sqrt(6);

    cameraDirZ[0] = -1.0/sqrt(3);
    cameraDirZ[1] = -1.0/sqrt(3);
    cameraDirZ[2] = -1.0/sqrt(3);
    // Get transform matrices:
    // Extrinsic
    float* R =  new float[12];
    R[0] = cameraDirX[0];
    R[1] = cameraDirX[1];
    R[2] = cameraDirX[2];
    MatrixMul(cameraLoc,1,3,cameraDirX,3,1,&R[3]);
    R[3] = -R[3];
    R[4] = cameraDirY[0];
    R[5] = cameraDirY[1];
    R[6] = cameraDirY[2];
    MatrixMul(cameraLoc,1,3,cameraDirY,3,1,&R[7]);
    R[7] = -R[7];
    R[8] = cameraDirZ[0];
    R[9] = cameraDirZ[1];
    R[10] = cameraDirZ[2];
    MatrixMul(cameraLoc,1,3,cameraDirZ,3,1,&R[11]);
    R[11] = -R[11];

    // Intrinsic
    float f = sqrt(3);
    float cx = 0;
    float cy = 0;
    float* K = new float[9];
    K[0] = f;
    K[1] = 0;
    K[2] = cx;
    K[3] = 0;
    K[4] = f;
    K[5] = cy;
    K[6] = 0;
    K[7] = 0;
    K[8] = 1;

    // transformation matrix:
    float *T = new float[12];
    MatrixMul(K,3,3,R,3,4,T);
    float *rM = new float[16]; // additional rotation matrix
    float theta = 0;
    rM[0]=cos(theta);rM[1]=-sin(theta);rM[2]=0;rM[3]=0;
    rM[4]=sin(theta);rM[5]=cos(theta);rM[6]=0;rM[7]=0;
    rM[8]=0;rM[9]=0;rM[10]=1;rM[11]=0;
    rM[12]=0;rM[13]=0;rM[14]=0;rM[15]=1;

    float* depth = new float[sizeof(float)*SizeCapture*SizeCapture];
    unsigned char* capture = Capture(src1,src2,src3,src4,src5,Size,pDensity,SizeCapture,T,cameraLoc,rM,depth);
    WriteImageToFile("B_capture_out.raw",capture,SizeCapture,SizeCapture,3);

    // Reverse mapping:
    float R_n[16];
    for(int i=0;i<12;i++)
        R_n[i] = R[i];
    R_n[12] = 0;R_n[13] = 0;R_n[14] = 0;R_n[15] = 1;
    Mat R_nM(4,4,CV_32F,&R_n);
    Mat RinvM = R_nM.inv();
    float* Rinv = (float*)RinvM.ptr();

    printf("Reverse mapping (i,j)->(X,Y,Z):\n");
    while(1)
    {
        printf("(i,j)(%dx%d): ",SizeCapture,SizeCapture);
        int i,j;
        cin>>i>>j;
        int offset = SizeCapture/2.0;
        float Pos_in[4];
        float Pos_out[4];
        Pos_in[2] = depth[i*SizeCapture+j];
        Pos_in[0] = 1.0*(j-offset)/pDensity/f*Pos_in[2];
        Pos_in[1] = 1.0*(i-offset)/pDensity/f*Pos_in[2];
        Pos_in[3] = 1;
        MatrixMul(Rinv,4,4,Pos_in,4,1,Pos_out);

        printf("(X,Y,Z): %.1f, %.1f, %.1f\n",Pos_out[0],Pos_out[1],Pos_out[2]);
    }

    R_nM.release();
    RinvM.release();
    delete capture;
    delete R;
    delete K;
    delete T;
    delete rM;
    delete depth;
    printf(" B ends...\n");
}

// (C) Rotating Cube
void partC(unsigned char* src1,unsigned char* src2,unsigned char* src3,unsigned char* src4,unsigned char* src5,int Size,int pDensity,int SizeCapture)
{
    printf(" C begins...\n");

    float cameraLoc[3];
    float cameraDirX[3];
    float cameraDirY[3];
    float cameraDirZ[3];

    cameraLoc[0] = 5;
    cameraLoc[1] = 5;
    cameraLoc[2] = 5;

    cameraDirX[0] = -1.0/sqrt(2);
    cameraDirX[1] = 1.0/sqrt(2);
    cameraDirX[2] = 0;

    cameraDirY[0] = 1.0/sqrt(6);
    cameraDirY[1] = 1.0/sqrt(6);
    cameraDirY[2] = -2.0/sqrt(6);

    cameraDirZ[0] = -1.0/sqrt(3);
    cameraDirZ[1] = -1.0/sqrt(3);
    cameraDirZ[2] = -1.0/sqrt(3);
    // Get transform matrices:
    // Extrinsic
    float* R =  new float[12];
    R[0] = cameraDirX[0];
    R[1] = cameraDirX[1];
    R[2] = cameraDirX[2];
    MatrixMul(cameraLoc,1,3,cameraDirX,3,1,&R[3]);
    R[3] = -R[3];
    R[4] = cameraDirY[0];
    R[5] = cameraDirY[1];
    R[6] = cameraDirY[2];
    MatrixMul(cameraLoc,1,3,cameraDirY,3,1,&R[7]);
    R[7] = -R[7];
    R[8] = cameraDirZ[0];
    R[9] = cameraDirZ[1];
    R[10] = cameraDirZ[2];
    MatrixMul(cameraLoc,1,3,cameraDirZ,3,1,&R[11]);
    R[11] = -R[11];

    // Intrinsic
    float f = sqrt(3);
    float cx = 0;
    float cy = 0;
    float* K = new float[9];
    K[0] = f;
    K[1] = 0;
    K[2] = cx;
    K[3] = 0;
    K[4] = f;
    K[5] = cy;
    K[6] = 0;
    K[7] = 0;
    K[8] = 1;

    // Final transformation matrix:
    float *T = new float[12];
    MatrixMul(K,3,3,R,3,4,T);
    float *rM = new float[16]; // additional rotation matrix
    rM[2]=0;rM[3]=0;
    rM[6]=0;rM[7]=0;
    rM[8]=0;rM[9]=0;rM[10]=1;rM[11]=0;
    rM[12]=0;rM[13]=0;rM[14]=0;rM[15]=1;
    unsigned char* capture;
    char filename[30];
    float* depth = new float[sizeof(float)*SizeCapture*SizeCapture];
    for(int theta=1;theta<=90;theta++)
    {
        sprintf(filename,"C_capture_%d_out.raw",theta);
        printf("%s\n",filename);
        fflush(stdout);

        float theta_r = PI*theta/180;
        rM[0] = cos(theta_r);
        rM[1] = -sin(theta_r);
        rM[4] = sin(theta_r);
        rM[5] = cos(theta_r);

        capture = Capture(src1,src2,src3,src4,src5,Size,pDensity,SizeCapture,T,cameraLoc,rM,depth);
        WriteImageToFile(filename,capture,SizeCapture,SizeCapture,3);
    }

    delete depth;
    delete capture;
    delete R;
    delete K;
    delete T;
    delete rM;
    printf(" C ends...\n");
}


int main()
{
    int Size = 200;

    unsigned char *src1 = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    ReadImageFromFile("image1.raw",src1,Size,Size,3);
    unsigned char *src2 = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    ReadImageFromFile("image2.raw",src2,Size,Size,3);
    unsigned char *src3 = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    ReadImageFromFile("image3.raw",src3,Size,Size,3);
    unsigned char *src4 = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    ReadImageFromFile("image4.raw",src4,Size,Size,3);
    unsigned char *src5 = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    ReadImageFromFile("image5.raw",src5,Size,Size,3);

    // A:
    partA(src1,src2,src3,src4,src5,Size);

    int pDensity = 200;
    int SizeCapture = 200;
    // B:
    //partB(src1,src2,src3,src4,src5,Size,pDensity,SizeCapture); //Change pDensity, output image size = 2*pDensity

    // C:
    //partC(src1,src2,src3,src4,src5,Size,pDensity,SizeCapture);

    free (src1);
    free (src2);
    free (src3);
    free (src4);
    free (src5);
    return 0;
}

