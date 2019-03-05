/**
EE569 HW #4 Problem 1:Face Warping
    Name:	Manan Vyas
	USC ID:	7483-8632-00
	email:	mvyas@usc.edu
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "IO.h"
#include "Filters.h"
#include "ImageProcessingFunctions.h"
#include <math.h>
#include <math.h>
#include "opencv2/core/core.hpp"


#define PI 3.1415926

using namespace cv;
using namespace std;

typedef struct pointPair{
    float fromI;
    float fromJ;
    float toI;
    float toJ;
} pointPair;
typedef struct fPoint{
    float x;
    float y;
} fPoint;

// Part (a) Back to Baby
float sign(fPoint p1, fPoint p2, fPoint p3)
{
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}
bool PointInTriangle(fPoint pt, fPoint v1, fPoint v2, fPoint v3)
{
    bool b1, b2, b3;

    b1 = sign(pt, v1, v2) < 0.0f;
    b2 = sign(pt, v2, v3) < 0.0f;
    b3 = sign(pt, v3, v1) < 0.0f;

    return ((b1 == b2) && (b2 == b3));
}
void getTargetCPs(char *srcfile1, char *srcfile2, float alpha, char* filename)
{
    ofstream outfile(filename);
    ifstream infile1(srcfile1);
    ifstream infile2(srcfile2);
    int PointsNUM;
    infile1>>PointsNUM;
    infile2>>PointsNUM;
    outfile<<PointsNUM<<endl;
    pointPair *ControlPoints = new pointPair[PointsNUM];
    for(int k=0;k<PointsNUM;k++)
    {
        infile1>>ControlPoints[k].fromI;
        infile1>>ControlPoints[k].fromJ;
        infile2>>ControlPoints[k].toI;
        infile2>>ControlPoints[k].toJ;
        outfile<<1.0*(1-alpha)*ControlPoints[k].fromI + 1.0*alpha*ControlPoints[k].toI<<" ";
        outfile<<1.0*(1-alpha)*ControlPoints[k].fromJ + 1.0*alpha*ControlPoints[k].toJ<<endl;
    }
    infile1.close();
    infile2.close();
    outfile.close();
    delete ControlPoints;
}
unsigned char* imageWarping(unsigned char* src, int SizeI, int SizeJ,char* srccpfile,char* tarcpfile,char* trianglefile)
{
    int CP_NUM = 3;
    int BytesPerPixel = 3;
    unsigned char *tar = new unsigned char[SizeI*SizeJ*3];
    ifstream infileP(srccpfile);
    ifstream infilePt(tarcpfile);
    int PointsNUM;
    infileP>>PointsNUM;
    infilePt>>PointsNUM;
    pointPair *ControlPoints = new pointPair[PointsNUM];
    for(int k=0;k<PointsNUM;k++)
    {
        infileP>>ControlPoints[k].fromI;
        infileP>>ControlPoints[k].fromJ;
        infilePt>>ControlPoints[k].toI;
        infilePt>>ControlPoints[k].toJ;
    }
    infileP.close();
    infilePt.close();

    ifstream infileT(trianglefile);
    int TRI_NUM;
    infileT>>TRI_NUM;
    Mat Tinv[TRI_NUM];
    // Input output point pairs
    float arr_i[TRI_NUM][6];
    float arr_o[TRI_NUM][6];
    for(int k=0;k<TRI_NUM;k++)
    {
        int a,b,c;
        infileT>>a>>b>>c;
        arr_i[k][0] = ControlPoints[a].fromI; arr_i[k][1] = ControlPoints[a].fromJ;
        arr_o[k][0] = ControlPoints[a].toI; arr_o[k][1] = ControlPoints[a].toJ;
        arr_i[k][2] = ControlPoints[b].fromI; arr_i[k][3] = ControlPoints[b].fromJ;
        arr_o[k][2] = ControlPoints[b].toI; arr_o[k][3] = ControlPoints[b].toJ;
        arr_i[k][4] = ControlPoints[c].fromI; arr_i[k][5] = ControlPoints[c].fromJ;
        arr_o[k][4] = ControlPoints[c].toI; arr_o[k][5] = ControlPoints[c].toJ;
    }
    infileT.close();
    for(int k=0;k<TRI_NUM;k++)
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

    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {

            float X = 0;
            float Y = 0;
            Mat TMatrixInv;

            fPoint pt;
            pt.x = i; pt.y = j;
            for(int k=0;k<TRI_NUM;k++)
            {
                fPoint p1,p2,p3;
                p1.x = arr_o[k][0]; p1.y = arr_o[k][1];
                p2.x = arr_o[k][2]; p2.y = arr_o[k][3];
                p3.x = arr_o[k][4]; p3.y = arr_o[k][5];
                if(PointInTriangle(pt,p1,p2,p3))
                {
                    TMatrixInv = Tinv[k];
                }
            }
            if(!TMatrixInv.empty())
            {   X = i*TMatrixInv.at<float>(0)+j*TMatrixInv.at<float>(1)+TMatrixInv.at<float>(2);
                Y = i*TMatrixInv.at<float>(3)+j*TMatrixInv.at<float>(4)+TMatrixInv.at<float>(5);
            }
            else
            {
                X = i;
                Y = j;
            }

            int p,q;
            p = (int)X;
            q = (int)Y;
            if(p<0 || p>SizeI-1 || q<0 || q>SizeJ-1)
            {
                for(int color=0; color<BytesPerPixel; color++)
                    tar[BytesPerPixel*(i*SizeJ+j)+color] = 0;
            }
            else
            {
                // Bilinear interpolation
                if(p==SizeI-1 || q==SizeJ-1)
                {
                    for(int color=0; color<BytesPerPixel; color++)
                        tar[BytesPerPixel*(i*SizeJ+j)+color] = src[BytesPerPixel*(p*SizeJ+q)+color];
                }
                else
                {
                    float a = X-p;
                    float b = Y-q;
                    for(int color=0; color<BytesPerPixel; color++)
                    {
                        float tmp = (1-a)*(1-b)*src[(p*SizeJ+q)*BytesPerPixel+color]+(1-a)*b*src[(p*SizeJ+q+1)*BytesPerPixel+color]+a*(1-b)*src[((p+1)*SizeJ+q)*BytesPerPixel+color]+a*b*src[((p+1)*SizeJ+q+1)*BytesPerPixel+color];
                        tar[BytesPerPixel*(i*SizeJ+j)+color] = tmp;
                    }
                }
            }
            TMatrixInv.release();
        }
    }

    delete ControlPoints;
    return tar;
}

void partA()
{
    printf(" A begins...\n");
    int SizeI = 350;
    int SizeJ = 300;
    char filename0[30] = "";
    sprintf(filename0,"controlPoints_old.txt");
    char filename1[30] = "";
    sprintf(filename1,"controlPoints_target.txt");
    char filename2[30] = "";
    sprintf(filename2,"controlTriangles.txt");
    char filename3[30] = "";
    sprintf(filename3,"controlPoints_young.txt");

    getTargetCPs(filename0,filename3,0.5,filename1);

    unsigned char* src = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("old_drew.raw",src,SizeI,SizeJ,3);
    unsigned char* tar1 = imageWarping(src,SizeI,SizeJ,filename0,filename1,filename2);
    WriteImageToFile("old_drew_warping.raw",tar1,SizeI,SizeJ,3);


    ReadImageFromFile("young_drew.raw",src,SizeI,SizeJ,3);
    unsigned char* tar2 = imageWarping(src,SizeI,SizeJ,filename3,filename1,filename2);
    WriteImageToFile("young_drew_warping.raw",tar2,SizeI,SizeJ,3);

    delete src;
    delete tar1;
    delete tar2;
    printf(" A ends...\n");
}

void partB()
{
    printf(" B begins...\n");
    int SizeI = 350;
    int SizeJ = 300;
    unsigned char* src1 = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("old_drew_warping.raw",src1,SizeI,SizeJ,3);
    unsigned char* src2 = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("young_drew_warping.raw",src2,SizeI,SizeJ,3);

    unsigned char* tar = new unsigned char[SizeI*SizeJ*3];

    for(int k=0;k<=100;k++)
    {
        float alpha = 1.0*k/100;
        for(int i=0;i<SizeI;i++)
        {
            for(int j=0;j<SizeJ;j++)
            {
                for(int color=0;color<3;color++)
                    tar[3*(i*SizeJ+j)+color] = (1-alpha)*src1[3*(i*SizeJ+j)+color]+alpha*src2[3*(i*SizeJ+j)+color];
            }
        }
        char filename[10];
        sprintf(filename,"%d.raw",k);
        WriteImageToFile(filename,tar,SizeI,SizeJ,3);
    }

    delete src1;
    delete src2;
    delete tar;
    printf(" B ends...\n");
}

void partC1()
{
    printf(" C1 begins...\n");
    int SizeI = 400;
    int SizeJ = 400;
    char filename0[30] = "";
    sprintf(filename0,"controlPoints_bruce.txt");
    char filename1[30] = "";
    sprintf(filename1,"controlPoints_target.txt");
    char filename2[30] = "";
    sprintf(filename2,"controlTriangles.txt");
    char filename3[30] = "";
    sprintf(filename3,"controlPoints_hulk.txt");

    getTargetCPs(filename0,filename3,0.5,filename1);

    unsigned char* src = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("bruce_banner.raw",src,SizeI,SizeJ,3);
    unsigned char* tar1 = imageWarping(src,SizeI,SizeJ,filename0,filename1,filename2);
    WriteImageToFile("bruce_warping.raw",tar1,SizeI,SizeJ,3);

    ReadImageFromFile("hulk.raw",src,SizeI,SizeJ,3);
    unsigned char* tar2 = imageWarping(src,SizeI,SizeJ,filename3,filename1,filename2);
    WriteImageToFile("hulk_warping.raw",tar2,SizeI,SizeJ,3);

    delete src;
    delete tar1;
    delete tar2;
    printf(" C1 ends...\n");
}
void partC2()
{
    printf(" C2 begins...\n");
    int SizeI = 400;
    int SizeJ = 400;
    unsigned char* src1 = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("bruce_warping.raw",src1,SizeI,SizeJ,3);
    unsigned char* src2 = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("hulk_warping.raw",src2,SizeI,SizeJ,3);

    unsigned char* tar = new unsigned char[SizeI*SizeJ*3];

    for(int k=0;k<=100;k++)
    {
        float alpha = 1.0*k/100;
        for(int i=0;i<SizeI;i++)
        {
            for(int j=0;j<SizeJ;j++)
            {
                for(int color=0;color<3;color++)
                    tar[3*(i*SizeJ+j)+color] = (1-alpha)*src1[3*(i*SizeJ+j)+color]+alpha*src2[3*(i*SizeJ+j)+color];
            }
        }
        char filename[10];
        sprintf(filename,"%d.raw",k);
        WriteImageToFile(filename,tar,SizeI,SizeJ,3);
    }

    delete src1;
    delete src2;
    delete tar;
    printf(" C2 ends...\n");
}

void partC3()
{
    printf(" C3 begins...\n");
    /*
    //old-young
    int SizeI = 350;
    int SizeJ = 300;
    unsigned char* srcface = new unsigned char[SizeI*SizeJ*3];
    unsigned char* tarface = new unsigned char[SizeI*SizeJ*3];
    unsigned char* tar = new unsigned char[SizeI*SizeJ*3];
    char filename0[30] = "";
    sprintf(filename0,"controlPoints_old.txt");
    char filename1[30] = "";
    sprintf(filename1,"controlPoints_target_C.txt");
    char filename2[30] = "";
    sprintf(filename2,"controlTriangles.txt");
    char filename3[30] = "";
    sprintf(filename3,"controlPoints_young.txt");

    ReadImageFromFile("old_drew.raw",srcface,SizeI,SizeJ,3);
    ReadImageFromFile("young_drew.raw",tarface,SizeI,SizeJ,3);
    */

    // bruce-hulk
    int SizeI = 400;
    int SizeJ = 400;
    unsigned char* srcface = new unsigned char[SizeI*SizeJ*3];
    unsigned char* tarface = new unsigned char[SizeI*SizeJ*3];
    unsigned char* tar = new unsigned char[SizeI*SizeJ*3];
    char filename0[30] = "";
    sprintf(filename0,"controlPoints_bruce.txt");
    char filename1[30] = "";
    sprintf(filename1,"controlPoints_target_C.txt");
    char filename2[30] = "";
    sprintf(filename2,"controlTriangles.txt");
    char filename3[30] = "";
    sprintf(filename3,"controlPoints_hulk.txt");

    ReadImageFromFile("bruce_banner.raw",srcface,SizeI,SizeJ,3);
    ReadImageFromFile("hulk.raw",tarface,SizeI,SizeJ,3);


    for(int k=0;k<=100;k++)
    {
        cout<<k<<endl;
        fflush(stdout);
        float alpha = 1.0*k/100;
        getTargetCPs(filename0,filename3,alpha,filename1);
        unsigned char* tar1 = imageWarping(srcface,SizeI,SizeJ,filename0,filename1,filename2);
        unsigned char* tar2 = imageWarping(tarface,SizeI,SizeJ,filename3,filename1,filename2);
        for(int i=0;i<SizeI;i++)
        {
            for(int j=0;j<SizeJ;j++)
            {
                for(int color=0;color<3;color++)
                    tar[3*(i*SizeJ+j)+color] = (1-alpha)*tar1[3*(i*SizeJ+j)+color]+alpha*tar2[3*(i*SizeJ+j)+color];
            }
        }
        char filename[10];
        sprintf(filename,"%d.raw",k);
        WriteImageToFile(filename,tar,SizeI,SizeJ,3);
        delete tar1;
    }

    delete srcface;
    delete tarface;
    delete tar;
    printf(" C3 ends...\n");
}
int main()
{
//    partA();
//    partB();
    partC1();
    partC2();
//    partC3();
    return 0;
}
