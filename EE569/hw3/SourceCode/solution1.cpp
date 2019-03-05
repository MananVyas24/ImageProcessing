/**
EE569 HW #3 Problem 1:Spatial Warping
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
#include <math.h>
#include "opencv2/core/core.hpp"

#define PI 3.1415926
#define CP_NUM 3

using namespace cv;
using namespace std;


// Part (a):Warping to Diamond Shape
// (1) using 8 triangles
void partA1()
{
    printf(" A1 begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    unsigned char *tar = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    ReadImageFromFile("puppy.raw",src,Size,Size,3);

    float center = 1.0*(Size-1)/2;
    float ed = Size-1;

    Mat Tinv[8];

    float arr_i[8][CP_NUM*2] = {{center,center,0,0,0,center},
                                {center,center,0,ed,0,center},
                                {center,center,0,ed,center,ed},
                                {center,center,center,ed,ed,ed},
                                {center,center,ed,ed,ed,center},
                                {center,center,ed,center,ed,0},
                                {center,center,ed,0,center,0},
                                {center,center,0,0,center,0}};

    float arr_o[8][CP_NUM*2] = {{center,center,center/2,center/2,0,center},
                                {center,center,center/2,(center+ed)/2,0,center},
                                {center,center,center/2,(center+ed)/2,center,ed},
                                {center,center,center,ed,(center+ed)/2,(center+ed)/2},
                                {center,center,(center+ed)/2,(center+ed)/2,ed,center},
                                {center,center,ed,center,(center+ed)/2,center/2},
                                {center,center,(center+ed)/2,center/2,center,0},
                                {center,center,center/2,center/2,center,0}};

    for(int k=0;k<8;k++)
    {
        Mat input(CP_NUM,2,CV_32F,&arr_i[k][0]);

        Mat A(CP_NUM*2,6,CV_32F,Scalar(0));
        for(int i=0; i<CP_NUM; i++)
        {
            A.at<float>(i*2,0) = input.at<float>(i*2);
            A.at<float>(i*2,1) = input.at<float>(i*2+1);
            A.at<float>(i*2,2) = 1;
            A.at<float>(i*2+1,3) = input.at<float>(i*2);
            A.at<float>(i*2+1,4) = input.at<float>(i*2+1);
            A.at<float>(i*2+1,5) = 1;
        }

        Mat B(CP_NUM*2,1,CV_32F,&arr_o[k][0]);

        Mat T = A.inv(DECOMP_SVD) * B;
        T = T.reshape(0,2);
        float tmp[] = {0,0,1};
        Mat tM(1,3,CV_32F,&tmp);
        T.push_back(tM);

        Tinv[k] = T.inv(DECOMP_SVD);

        input.release();
        A.release();
        T.release();
        B.release();
        tM.release();
    }
    // TMatrixInv is 3x3 and is used to find input point pos
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            float X = 0;
            float Y = 0;
            Mat TMatrixInv;

            if(i<=Size/2 && j<=Size/2 && j>=i)
            {
                TMatrixInv = Tinv[0];
            }
            if(i<=Size/2 && j>=Size/2 && i+j<Size)
            {
                TMatrixInv = Tinv[1];
            }
            if(i<=Size/2 && j>=Size/2 && i+j>=Size)
            {
                TMatrixInv = Tinv[2];
            }
            if(i>=Size/2 && j>=Size/2 && j>=i)
            {
                TMatrixInv = Tinv[3];
            }
            if(i>=Size/2 && j>=Size/2 && j<i)
            {
                TMatrixInv = Tinv[4];
            }
            if(i>=Size/2 && j<=Size/2 && j+i>=Size)
            {
                TMatrixInv = Tinv[5];
            }
            if(i>=Size/2 && j<=Size/2 && j+i<Size)
            {
                TMatrixInv = Tinv[6];
            }
            if(i<=Size/2 && j<=Size/2 && j<i)
            {
                TMatrixInv = Tinv[7];
            }

            X = i*TMatrixInv.at<float>(0)+j*TMatrixInv.at<float>(1)+TMatrixInv.at<float>(2);
            Y = i*TMatrixInv.at<float>(3)+j*TMatrixInv.at<float>(4)+TMatrixInv.at<float>(5);
            if(X<0 || X>Size-1 || Y<0 || Y>Size-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*Size+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                int p,q;
                p = (int)X;
                q = (int)Y;
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    float a = X-p;
                    float b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        float tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("puppy_1_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" A1 ends...\n");
}

void partA1_Reverse()
{
    printf(" A1 reverse begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    unsigned char *tar = (unsigned char*)malloc(sizeof(unsigned char)*Size*Size*3);
    ReadImageFromFile("puppy_1_out.raw",src,Size,Size,3);

    float center = 1.0*(Size-1)/2;
    float ed = Size-1;

    Mat Tinv[8];

    float arr_i[8][CP_NUM*2] = {{center,center,0,0,0,center},
                                {center,center,0,ed,0,center},
                                {center,center,0,ed,center,ed},
                                {center,center,center,ed,ed,ed},
                                {center,center,ed,ed,ed,center},
                                {center,center,ed,center,ed,0},
                                {center,center,ed,0,center,0},
                                {center,center,0,0,center,0}};

    float arr_o[8][CP_NUM*2] = {{center,center,center/2,center/2,0,center},
                                {center,center,center/2,(center+ed)/2,0,center},
                                {center,center,center/2,(center+ed)/2,center,ed},
                                {center,center,center,ed,(center+ed)/2,(center+ed)/2},
                                {center,center,(center+ed)/2,(center+ed)/2,ed,center},
                                {center,center,ed,center,(center+ed)/2,center/2},
                                {center,center,(center+ed)/2,center/2,center,0},
                                {center,center,center/2,center/2,center,0}};

    for(int k=0;k<8;k++)
    {
        Mat input(CP_NUM,2,CV_32F,&arr_i[k][0]);

        Mat A(CP_NUM*2,6,CV_32F,Scalar(0));
        for(int i=0; i<CP_NUM; i++)
        {
            A.at<float>(i*2,0) = input.at<float>(i*2);
            A.at<float>(i*2,1) = input.at<float>(i*2+1);
            A.at<float>(i*2,2) = 1;
            A.at<float>(i*2+1,3) = input.at<float>(i*2);
            A.at<float>(i*2+1,4) = input.at<float>(i*2+1);
            A.at<float>(i*2+1,5) = 1;
        }

        Mat B(CP_NUM*2,1,CV_32F,&arr_o[k][0]);

        Mat T = A.inv(DECOMP_SVD) * B;
        T = T.reshape(0,2);
        float tmp[] = {0,0,1};
        Mat tM(1,3,CV_32F,&tmp);
        T.push_back(tM);

        Tinv[k] = T;

        input.release();
        A.release();
        T.release();
        B.release();
        tM.release();
    }
    // TMatrixInv is 3x3 and is used to find input point pos
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            float X = 0;
            float Y = 0;
            Mat TMatrixInv;

            if(i<=Size/2 && j<=Size/2 && j>=i)
            {
                TMatrixInv = Tinv[0];
            }
            if(i<=Size/2 && j>=Size/2 && i+j<Size)
            {
                TMatrixInv = Tinv[1];
            }
            if(i<=Size/2 && j>=Size/2 && i+j>=Size)
            {
                TMatrixInv = Tinv[2];
            }
            if(i>=Size/2 && j>=Size/2 && j>=i)
            {
                TMatrixInv = Tinv[3];
            }
            if(i>=Size/2 && j>=Size/2 && j<i)
            {
                TMatrixInv = Tinv[4];
            }
            if(i>=Size/2 && j<=Size/2 && j+i>=Size)
            {
                TMatrixInv = Tinv[5];
            }
            if(i>=Size/2 && j<=Size/2 && j+i<Size)
            {
                TMatrixInv = Tinv[6];
            }
            if(i<=Size/2 && j<=Size/2 && j<i)
            {
                TMatrixInv = Tinv[7];
            }

            X = i*TMatrixInv.at<float>(0)+j*TMatrixInv.at<float>(1)+TMatrixInv.at<float>(2);
            Y = i*TMatrixInv.at<float>(3)+j*TMatrixInv.at<float>(4)+TMatrixInv.at<float>(5);
            if(X<0 || X>Size-1 || Y<0 || Y>Size-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*Size+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                int p,q;
                p = (int)X;
                q = (int)Y;
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("puppy_1_reverse_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" A1 reverse ends...\n");
}
void partA2()
{
    printf(" A2 begins...\n");
    int Size = 500;
    float center = 1.0*(Size-1)/2;
    float centerl = center-1;
    float centerr = center+1;
    float ed = Size-1;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];

    ReadImageFromFile("puppy.raw",src,Size,Size,3);

    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            for(int color=0; color<BytesPerPixel; color++)
                tar[BytesPerPixel*(i*Size+j)+color] = 0;
            if(i<center && centerl*i+center*j>=center*centerl && (ed-centerr)*i>=(j-centerr)*center)
            {
                double lb = centerl-i*centerl/center; // left boundary
                double len = 2*(center - lb);
                double t = Size*(j-lb)/len;
                int q = t;
                double a = t-q;
                if(q==ed)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(i*Size+q)+color];
                }
                else
                {
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*src[(i*Size+q)*BytesPerPixel+color] + a*src[(i*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
            if(i>=center && centerl*(i-center)<=(ed-center)*j && (i-center)*(ed-centerr)<=(ed-j)*(ed-center))
            {
                double lb = centerl*(i-center)/(ed-center); // left boundary
                double len = 2*(center - lb);
                double t = Size*(j-lb)/len;
                int q = t;
                double a = t-q;
                if(q==ed)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(i*Size+q)+color];
                }
                else
                {
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*src[(i*Size+q)*BytesPerPixel+color] + a*src[(i*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }

    WriteImageToFile("puppy_2_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" A2 ends...\n");
}
void partA2_Reverse()
{
    printf(" A2 reverse begins...\n");
    int Size = 500;
    float center = 1.0*(Size-1)/2;
    float centerl = center-1;
    float ed = Size-1;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];

    ReadImageFromFile("puppy_2_out.raw",src,Size,Size,3);

    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            for(int color=0; color<BytesPerPixel; color++)
                tar[BytesPerPixel*(i*Size+j)+color] = 0;
            if(i<center)
            {
                double lb = centerl-i*centerl/center; // left boundary
                double len = 2*(center - lb);
                double t = 1.0*j/Size*len+lb;
                int q = t;
                double a = t-q;
                if(q==ed)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(i*Size+q)+color];
                }
                else
                {
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*src[(i*Size+q)*BytesPerPixel+color] + a*src[(i*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
            if(i>=center)
            {
                double lb = centerl*(i-center)/(ed-center); // left boundary
                double len = 2*(center - lb);
                double t = 1.0*j/Size*len+lb;
                int q = t;
                double a = t-q;
                if(q==ed)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(i*Size+q)+color];
                }
                else
                {
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*src[(i*Size+q)*BytesPerPixel+color] + a*src[(i*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }

    WriteImageToFile("puppy_2_reverse_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" A2 reverse ends...\n");
}

// Part (b):Warping to Pentagon Shape
void partB1()
{
    printf(" B1 begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];

    ReadImageFromFile("cowboy.raw",src,Size,Size,3);

    double TMatrixInv_1[9] = {1,1,-250,0,2,-250,0,0,1};
    double TMatrixInv_2[9] = {1,-1,250,0,1.992,-248,0,0,1};
    double TMatrixInv_3[9] = {2,0,-250,-0.992,1,248,0,0,1};
    double TMatrixInv_4[9] = {1.24,0,-61.25,0.62,1,-155,0,0,1};
    double TMatrixInv_5[9] = {1.245,0,-61.25,-0.0025,1.996,-248.375,0,0,1};
    double TMatrixInv_6[9] = {1.245,0,-61.25,-0.625,1,156.25,0,0,1};
    double TMatrixInv_7[9] = {2,0,-250,1,1,-250,0,0,1};

    // TMatrixInv is 3x3 and is used to find input point pos
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            double X = 0;
            double Y = 0;
            double *TMatrixInv;

            if(i<=Size/2 && j<=Size/2 && j>=i)
            {
                TMatrixInv = TMatrixInv_1;
            }
            if(i<=Size/2 && j>=Size/2 && i+j<Size)
            {
                TMatrixInv = TMatrixInv_2;
            }
            if(i<=Size/2 && j>=Size/2 && i+j>=Size)
            {
                TMatrixInv = TMatrixInv_3;
            }
            if(i>Size/2 && 5*i-8*j+750<0)
            {
                TMatrixInv = TMatrixInv_4;
            }
            if(i>Size/2 && 5*i-8*j+750>=0 && 5*i+8*j-3250>=0)
            {
                TMatrixInv = TMatrixInv_5;
            }
            if(i>Size/2 && 5*i+8*j-3250<0)
            {
                TMatrixInv = TMatrixInv_6;
            }
            if(i<=Size/2 && j<=Size/2 && j<i)
            {
                TMatrixInv = TMatrixInv_7;
            }

            X = i*TMatrixInv[0]+j*TMatrixInv[1]+TMatrixInv[2];
            Y = i*TMatrixInv[3]+j*TMatrixInv[4]+TMatrixInv[5];
            if(X<0 || X>Size-1 || Y<0 || Y>Size-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*Size+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                int p,q;
                p = (int)X;
                q = (int)Y;
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("cowboy_1_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" B1 ends...\n");
}
void partB1_Reverse()
{
    printf(" B1 reverse begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];

    ReadImageFromFile("cowboy_1_out.raw",src,Size,Size,3);

    double TMatrixInv_1[9] = {1,-0.5,125,0,0.5,125,0,0,1};
    double TMatrixInv_2[9] = {1,0.502,-125.502,0,0.502,124.498,0,0,1};
    double TMatrixInv_3[9] = {0.5,0,125,0.496,1,-124,0,0,1};
    double TMatrixInv_4[9] = {0.8032,0,49.1968,-0.498,1,124.498,0,0,1};
    double TMatrixInv_5[9] = {0.8032,0,49.1968,0.001,0.501,124.498,0,0,1};
    double TMatrixInv_6[9] = {0.8032,0,49.1968,0.502,1,-125.502,0,0,1};
    double TMatrixInv_7[9] = {0.5,0,125,-0.5,1,125,0,0,1};

    // TMatrixInv is 3x3 and is used to find input point pos
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            double X = 0;
            double Y = 0;
            double *TMatrixInv;

            if(i<=Size/2 && j<=Size/2 && j>=i)
            {
                TMatrixInv = TMatrixInv_1;
            }
            if(i<=Size/2 && j>=Size/2 && i+j<Size)
            {
                TMatrixInv = TMatrixInv_2;
            }
            if(i<=Size/2 && j>=Size/2 && i+j>=Size)
            {
                TMatrixInv = TMatrixInv_3;
            }
            if(i>Size/2 && j>i)
            {
                TMatrixInv = TMatrixInv_4;
            }
            if(i>Size/2 && j<=i && j+i>=Size)
            {
                TMatrixInv = TMatrixInv_5;
            }
            if(i>Size/2 && j+i<Size)
            {
                TMatrixInv = TMatrixInv_6;
            }
            if(i<=Size/2 && j<=Size/2 && j<i)
            {
                TMatrixInv = TMatrixInv_7;
            }

            X = i*TMatrixInv[0]+j*TMatrixInv[1]+TMatrixInv[2];
            Y = i*TMatrixInv[3]+j*TMatrixInv[4]+TMatrixInv[5];
            if(X<0 || X>Size-1 || Y<0 || Y>Size-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*Size+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                int p,q;
                p = (int)X;
                q = (int)Y;
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("cowboy_1_inverse_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" B1 reverse ends...\n");
}

void partB2()
{
    printf(" B2 begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];
    ReadImageFromFile("cowboy.raw",src,Size,Size,3);

    float center = 1.0*(Size-1)/2;
    float cl = center-1;
    float cr = center+1;
    float ed = Size-1;
    float margin = Size/25.0;
    float tb = margin; //top boundary
    float bb = ed-margin;
    float lb = margin;
    float rb = ed-margin;
    float bl = (center+lb)/2; // bottom left
    float br = (center+ed)/2;

    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            for(int color=0; color<BytesPerPixel; color++)
                tar[BytesPerPixel*(i*Size+j)+color] = 0;
            if(i>=tb && i<=center && (i-tb)*(lb-cl)<=(j-cl)*(center-tb) && (i-tb)*(rb-cr)>=(j-cr)*(center-tb))
            {
                double base = (i-tb)*(lb-cl)/(center-tb)+cl; // left boundary
                double len = (i-tb)*(rb-cr)/(center-tb)+cr - base;
                float X = (i-tb)/(bb-tb)*Size;
                float Y = Size*(j-base)/len;
                int p = X;
                int q = Y;
                // bilinear interpolation
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }

            }
            if(i>center && i<=bb && (i-center)*(bl-lb)<=(bb-center)*(j-lb) && (i-center)*(br-rb)>=(bb-center)*(j-rb))
            {
                double base = (i-center)*(bl-lb)/(bb-center)+lb; // left boundary
                double len = (i-center)*(br-rb)/(bb-center)+rb - base;
                float X = (i-tb)/(bb-tb)*Size;
                float Y = Size*(j-base)/len;
                int p = X;
                int q = Y;
                // bilinear interpolation
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("cowboy_2_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" B2 ends...\n");
}
void partB2_Reverse()
{
    printf(" B2 reverse begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];
    ReadImageFromFile("cowboy_2_out.raw",src,Size,Size,3);

    float center = 1.0*(Size-1)/2;
    float cl = center-1;
    float cr = center+1;
    float ed = Size-1;
    float margin = Size/25.0;
    float tb = margin; //top boundary
    float bb = ed-margin;
    float lb = margin;
    float rb = ed-margin;
    float bl = (center+lb)/2; // bottom left
    float br = (center+ed)/2;

    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            for(int color=0; color<BytesPerPixel; color++)
                tar[BytesPerPixel*(i*Size+j)+color] = 0;
            if(i<=center)
            {
                double t = 1.0*i/Size*(bb-tb)+tb;
                double base = (t-tb)*(lb-cl)/(center-tb)+cl; // left boundary
                double len = (t-tb)*(rb-cr)/(center-tb)+cr - base;
                float X = t;
                float Y = 1.0*j/Size*len + base;
                int p = X;
                int q = Y;
                // bilinear interpolation
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }

            }
            if(i>center)
            {
                double t = 1.0*i/Size*(bb-tb)+tb;
                double base = (t-center)*(bl-lb)/(bb-center)+lb; // left boundary
                double len = (t-center)*(br-rb)/(bb-center)+rb - base;
                float X = t;
                float Y = 1.0*j/Size*len + base;
                int p = X;
                int q = Y;
                // bilinear interpolation
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("cowboy_2_reverse_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" B2 reverse ends...\n");
}

// Part (c):Warping to Circle Shape
// 2nd transform
void partC1()
{
    printf(" C1 begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];

    ReadImageFromFile("transformmer.raw",src,Size,Size,3);

    double TMatrixInv_1[12] = {-110.5695,1.4169,0.0415,-0.0000,-0.0001,0.0000,
                               -110.5695,-0.8185,2.2769,0.0016,-0.0017,-0.0000
                              };
    double TMatrixInv_2[12] = {-110.5695,2.2769,-0.8185,-0.0017,0.0016,0.0000,
                               -110.5695,0.0415,1.4169,-0.0001,-0.0000,0.0000
                              };
    double TMatrixInv_3[12] = {-104.7530,1.4235,0.0181,-0.0000,-0.0001,-0.0000,
                               -104.7530,0.8397,0.6018,-0.0017,0.0016,0
                              };
    double TMatrixInv_4[12] = {-104.7530,0.6018,0.8397,0.0016,-0.0017,0.0000,
                               -104.7530,0.0181,1.4235,-0.0001,-0.0000,-0.0000
                              };


    // TMatrixInv is 3x3 and is used to find input point pos
    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            double X = 0;
            double Y = 0;
            double *TMatrixInv;

            if(j>=Size/2 && i>= 250-(j-Size/2) && i<= 250+(j-Size/2))
            {
                TMatrixInv = TMatrixInv_1;
            }

            if(i>=Size/2 && j>= 250-(i-Size/2) && j<= 250+(i-Size/2))
            {
                TMatrixInv = TMatrixInv_2;
            }
            if(j<Size/2 && i>= 250-(Size/2-j) && i<= 250+(Size/2-j))
            {
                TMatrixInv = TMatrixInv_3;
            }
            if(i<Size/2 && j>= 250-(Size/2-i) && j<= 250+(Size/2-i))
            {
                TMatrixInv = TMatrixInv_4;
            }
            X = 1*TMatrixInv[0]+i*TMatrixInv[1]+j*TMatrixInv[2]+i*i*TMatrixInv[3]+j*j*TMatrixInv[4]+i*j*TMatrixInv[5];
            Y = 1*TMatrixInv[6]+i*TMatrixInv[7]+j*TMatrixInv[8]+i*i*TMatrixInv[9]+j*j*TMatrixInv[10]+i*j*TMatrixInv[11];
            if(X<0 || X>Size-1 || Y<0 || Y>Size-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*Size+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                int p,q;
                p = (int)X;
                q = (int)Y;
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("transformmer_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" C1 ends...\n");
}
// polar coordinates
void partC2()
{
    printf(" C2 begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];
    ReadImageFromFile("transformmer.raw",src,Size,Size,3);

    float center = 1.0*(Size-1)/2;
    float ed = Size-1;

    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            double l_uv,theta,L,l_xy,x,y,X,Y;
            double u = j-center;
            double v = center-i;

            l_uv = sqrt(u*u+v*v);
            if(l_uv==0)
            {
                X = 250;
                Y = 250;
            }
            else
            {
                theta = acos(1.0*u/l_uv);
                if(theta!=0)
                {
                    if(theta>PI/4 && theta<3*PI/4)
                        L = (ed/2)/fabs(sin(theta));
                    else
                        L = (ed/2)/fabs(cos(theta));
                    l_xy = L*l_uv/(Size/2);
                    x = fabs(l_xy*cos(theta));
                    y = fabs(l_xy*sin(theta));
                    if(i<=center&&j<=center)
                    {
                        X = center-y;
                        Y = center-x;
                    }
                    if(i>center&&j<=center)
                    {
                        X = center+y;
                        Y = center-x;
                    }
                    if(i<=center&&j>center)
                    {
                        X = center-y;
                        Y = center+x;
                    }
                    if(i>center&&j>center)
                    {
                        X = center+y;
                        Y = center+x;
                    }

                }
                else
                {
                    X = i;
                    Y = j;
                }
            }

            if(X<0 || X>Size-1 || Y<0 || Y>Size-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*Size+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                int p,q;
                p = (int)X;
                q = (int)Y;
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("transformmer_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" C2 ends...\n");
}

void partC2_Reverse()
{
    printf(" C2 reverse begins...\n");

    int Size = 500;
    int BytesPerPixel = 3;
    unsigned char *src = new unsigned char[Size*Size*3];
    unsigned char *tar = new unsigned char[Size*Size*3];

    ReadImageFromFile("transformmer_out.raw",src,Size,Size,3);

    float center = 1.0*(Size-1)/2;
    float ed = Size-1;

    for(int i=0; i<Size; i++)
    {
        for(int j=0; j<Size; j++)
        {
            double l_uv,theta,L,l_xy,x,y,X,Y;
            double u = j-center;
            double v = center-i;

            l_uv = sqrt(u*u+v*v);
            if(l_uv==0)
            {
                X = 250;
                Y = 250;
            }
            else
            {
                theta = acos(1.0*u/l_uv);
                if(theta!=0)
                {
                    if(theta>PI/4 && theta<3*PI/4)
                        L = (ed/2)/fabs(sin(theta));
                    else
                        L = (ed/2)/fabs(cos(theta));
                    l_xy = (ed/2)*l_uv/L;
                    x = fabs(l_xy*cos(theta));
                    y = fabs(l_xy*sin(theta));
                    if(i<=center&&j<=center)
                    {
                        X = center-y;
                        Y = center-x;
                    }
                    if(i>center&&j<=center)
                    {
                        X = center+y;
                        Y = center-x;
                    }
                    if(i<=center&&j>center)
                    {
                        X = center-y;
                        Y = center+x;
                    }
                    if(i>center&&j>center)
                    {
                        X = center+y;
                        Y = center+x;
                    }

                }
                else
                {
                    X = i;
                    Y = j;
                }
            }

            if(X<0 || X>Size-1 || Y<0 || Y>Size-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*Size+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                int p,q;
                p = (int)X;
                q = (int)Y;
                if(p==Size-1 || q==Size-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*Size+j)+color] = src[BytesPerPixel*(p*Size+q)+color];
                }
                else
                {
                    double a = X-p;
                    double b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        double tmp = (1-a)*(1-b)*src[(p*Size+q)*BytesPerPixel+color]+(1-a)*b*src[(p*Size+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*Size+q)*BytesPerPixel+color]+a*b*src[((p+1)*Size+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*Size+j)+color] = tmp;
                    }
                }
            }
        }
    }


    WriteImageToFile("transformmer_reverse_out.raw",tar,Size,Size,3);

    delete src;
    delete tar;
    printf(" C2 reverse ends...\n");
}


int main()
{
//    partA1();
//    partA1_Reverse();
//    partA2();
//    partA2_Reverse();
//    partB1();
//    partB1_Reverse();
//    partB2();
//    partB2_Reverse();
//    partC1();
    partC2();
    partC2_Reverse();
    return 0;
}
