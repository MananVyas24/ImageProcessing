/**
EE569 HW #3 Problem 3: Texture Analysis and Segmentation using Laws Filters
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
#include "opencv2/highgui/highgui.hpp"

#define PI 3.1415926

using namespace cv;
using namespace std;

extern "C" int *k_means(double**, int, int, int, double, double**);

// Part (A):Texture Image Clustering
void partA()
{
    printf(" A begins...\n");
    float F[3][5] = {{-1.0/6,-2.0/6,0,2.0/6,1.0/6},
        {-1.0/4,0,2.0/4,0,-1.0/4},
        {-1.0/6,2.0/6,0,-2.0/6,1.0/6}
    };
    float L[9][25]; // 9 laws filter
    for(int m=0; m<3; m++)
    {
        for(int n=0; n<3; n++)
        {
            for(int i=0; i<5; i++)
            {
                for(int j=0; j<5; j++)
                {
                    L[m*3+n][i*5+j] = F[m][i]*F[n][j];
                }
            }
        }
    }

    int Size = 128;
    unsigned char* src = new unsigned char[Size*Size];
    float *B = new float[Size*Size];

    int W = 15;
    int DataNum = 12;
    int DataDim = 9;
    int ClusterNum = 4;

    Mat featureVector(DataNum,DataDim,CV_32F);
    float feature[12][9];

    for(int f=0; f<12; f++)
    {
        char filename[20];
        sprintf(filename,"texture%d.raw",f+1);
        ReadImageFromFile(filename,src,Size,Size,1);

        for(int k=0; k<9; k++)
        {
            featureVector.at<float>(f,k) = 0;

            Filter2(&L[k][0],5,src,B,Size,Size,1);
            for(int i=0; i<Size; i++)
            {
                for(int j=0; j<Size; j++)
                {
                    double tmp = 0;
                    for(int m=-W/2; m<=W/2; m++)
                    {
                        for(int n=-W/2; n<=W/2; n++)
                        {
                            if(i+m>=0&&i+m<Size&&j+n>=0&&j+n<Size)
                            {
                                tmp += B[(i+m)*Size+j+n]*B[(i+m)*Size+j+n];
                            }
                        }
                    }
                    tmp /= W*W;
                    featureVector.at<float>(f,k) += tmp;
                }
            }
            featureVector.at<float>(f,k) /= (Size*Size);
            feature[f][k] = featureVector.at<float>(f,k);
        }
    }

/*
    while(1)
    {
        int i,j;
        cout<<"Input:";
        cin>>i>>j;
        cout<<"D= "<<GetDistance(&feature[i-1][0],&feature[j-1][0],9)<<endl;
    }
*/
//    int l[] = {0,1,2,0,0,0,3,0,0,0,0,0};
//    Mat labels(12,1,CV_32S,l);
    Mat labels;
    int attempts = 5;
    Mat centers;

    //kmeans(featureVector,ClusterNum,labels,TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001),attempts, KMEANS_USE_INITIAL_LABELS, centers);
    kmeans(featureVector,ClusterNum,labels,TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001),attempts, KMEANS_RANDOM_CENTERS, centers);
    for (int i = 0; i < DataNum; i++)
    {
        printf("Data point %d is in cluster %d\n", i+1, labels.at<int>(i));
    }

    delete src;
    delete B;
    featureVector.release();
    labels.release();
    centers.release();

    printf(" A ends...\n");
}

// Part (B):Texture Segmentation
void partB()
{
    printf(" B begins...\n");

    float F[3][3] = {{1.0/4,2.0/4,1.0/4},
        {-1.0/2,0,1.0/2},
        {1.0/4,-2.0/4,1.0/4}
    };
    float L[9][9]; // 9 laws filter
    for(int m=0; m<3; m++)
    {
        for(int n=0; n<3; n++)
        {
            for(int i=0; i<3; i++)
            {
                for(int j=0; j<3; j++)
                {
                    L[m*3+n][i*3+j] = F[m][i]*F[n][j];
                }
            }
        }
    }
    int SizeI = 450;
    int SizeJ = 600;
    unsigned char* src = new unsigned char[SizeI*SizeJ];
    unsigned char* output = new unsigned char[SizeI*SizeJ];
    float **B = new float*[9];
    for(int i=0; i<9; i++)
    {
        B[i] = new float[SizeI*SizeJ];
    }
    int DataNum = SizeI*SizeJ;
    int DataDim = 9;
    int ClusterNum = 6;

    Mat featureVector(DataNum,DataDim,CV_32F);

    ReadImageFromFile("comb.raw",src,SizeI,SizeJ,1);

    for(int i=0; i<9; i++)
    {
        Filter2(&L[i][0],3,src,B[i],SizeI,SizeJ,1);
    }

    unsigned char graylevel[6] = {0, 51, 102, 153, 204, 255};

    int W = 31;
    cout<<"Segmenting..."<<W<<endl;
    fflush(stdout);
    for(int i=0; i<SizeI; i++)
    {
        cout<<"line: "<<i<<endl;
        fflush(stdout);
        for(int j=0; j<SizeJ; j++)
        {
            for(int k=0; k<9; k++)
            {
                featureVector.at<float>(i*SizeJ+j,k) = 0;
                float tmp = 0;
                for(int m=-W/2; m<=W/2; m++)
                {
                    for(int n=-W/2; n<=W/2; n++)
                    {
                        if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                        {
                            tmp += B[k][(i+m)*SizeJ+j+n]*B[k][(i+m)*SizeJ+j+n];
                        }
                    }
                }
                tmp /= W*W;
                featureVector.at<float>(i*SizeJ+j,k) += tmp;
            }
        }
    }

    Mat labels;
    int attempts = 5;
    Mat centers;
    // Use kmeans function in Opencv
    kmeans(featureVector,ClusterNum,labels,TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001),attempts, KMEANS_RANDOM_CENTERS, centers);

    for(int i=0; i<DataNum; i++)
    {
        output[i] = graylevel[labels.at<int>(i)];
    }

    char filename[10];
    sprintf(filename,"B_%d.raw",W);
    WriteImageToFile(filename,output,SizeI,SizeJ,1);

    delete src;
    for(int i=0; i<9; i++)
        delete B[i];
    delete B;
    delete output;
    featureVector.release();
    labels.release();
    centers.release();

    printf(" B ends...\n");
}

void partC1()
{
    printf(" C1 begins...\n");

    float F[3][3] = {{1.0/4,2.0/4,1.0/4},
        {-1.0/2,0,1.0/2},
        {1.0/4,-2.0/4,1.0/4}
    };
    float L[9][9]; // 9 laws filter
    for(int m=0; m<3; m++)
    {
        for(int n=0; n<3; n++)
        {
            for(int i=0; i<3; i++)
            {
                for(int j=0; j<3; j++)
                {
                    L[m*3+n][i*3+j] = F[m][i]*F[n][j];
                }
            }
        }
    }
    int SizeI = 450;
    int SizeJ = 600;
    unsigned char* src = new unsigned char[SizeI*SizeJ];
    unsigned char* output = new unsigned char[SizeI*SizeJ];
    float **B = new float*[9];
    for(int i=0; i<9; i++)
    {
        B[i] = new float[SizeI*SizeJ];
    }
    int DataNum = SizeI*SizeJ;
    int FeatureDim = 11; // 9 for Laws filter
    int ClusterNum = 6;

    Mat featureVector(DataNum,FeatureDim,CV_32F);

    ReadImageFromFile("comb.raw",src,SizeI,SizeJ,1);

    for(int i=0; i<9; i++)
    {
        Filter2(&L[i][0],3,src,B[i],SizeI,SizeJ,1);
    }

    unsigned char graylevel[6] = {0, 51, 102, 153, 204, 255};

    int W = 51;
    cout<<"Segmenting..."<<W<<endl;
    fflush(stdout);
    float sums[11];
    for(int i=0;i<FeatureDim;i++)
    {
        sums[i] = 0;
    }
    for(int i=0; i<SizeI; i++)
    {
        cout<<"line: "<<i<<endl;
        fflush(stdout);
        for(int j=0; j<SizeJ; j++)
        {
            for(int k=0; k<9; k++)
            {
                featureVector.at<float>(i*SizeJ+j,k) = 0;
                float tmp = 0;
                for(int m=-W/2; m<=W/2; m++)
                {
                    for(int n=-W/2; n<=W/2; n++)
                    {
                        if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                        {
                            tmp += B[k][(i+m)*SizeJ+j+n]*B[k][(i+m)*SizeJ+j+n];
                        }
                    }
                }
                tmp /= W*W;
                featureVector.at<float>(i*SizeJ+j,k) += tmp;
                sums[k] += tmp;
            }
            featureVector.at<float>(i*SizeJ+j,9) = i;
            sums[9] += featureVector.at<float>(i*SizeJ+j,9);
            featureVector.at<float>(i*SizeJ+j,10) = j;
            sums[10] += featureVector.at<float>(i*SizeJ+j,10);

        }
    }
    for(int i=0;i<FeatureDim;i++)
    {
        sums[i] /= SizeI*SizeJ;
    }
    int weight = 10;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            for(int k=0; k<FeatureDim; k++)
            {
                if(sums[k]!=0)
                    featureVector.at<float>(i*SizeJ+j,k) /= sums[k];
                if(k==0)
                    featureVector.at<float>(i*SizeJ+j,k) *= weight;
            }
        }
    }

    //cout<<featureVector.row(1)<<endl;
    Mat labels;
    int attempts = 5;
    Mat centers;
    // Use kmeans function in Opencv
    kmeans(featureVector,ClusterNum,labels,TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001),attempts, KMEANS_RANDOM_CENTERS, centers);
    for(int i=0; i<DataNum; i++)
    {
        output[i] = graylevel[labels.at<int>(i)];
    }

    char filename[10];
    sprintf(filename,"C1_%d.raw",W);
    WriteImageToFile(filename,output,SizeI,SizeJ,1);

    delete src;
    for(int i=0; i<9; i++)
        delete B[i];
    delete B;
    delete output;
    featureVector.release();
    labels.release();
    centers.release();

    printf(" C1 ends...\n");
}

void partC2()
{
    printf(" C2 begins...\n");

    float F[5][5] = {{-1.0/6,-2.0/6,0,2.0/6,1.0/6},
        {-1.0/4,0,2.0/4,0,-1.0/4},
        {-1.0/6,2.0/6,0,-2.0/6,1.0/6},
        {1.0/16,4.0/16,6.0/16,4.0/16,1.0/16},
        {1.0/16,-4.0/16,6.0/16,-4.0/16,1.0/16}
    };
    float L[25][25]; // 25 laws filter
    for(int m=0; m<5; m++)
    {
        for(int n=0; n<5; n++)
        {
            for(int i=0; i<5; i++)
            {
                for(int j=0; j<5; j++)
                {
                    L[m*5+n][i*5+j] = F[m][i]*F[n][j];
                }
            }
        }
    }

    int SizeI = 450;
    int SizeJ = 600;
    unsigned char* src = new unsigned char[SizeI*SizeJ];
    unsigned char* output = new unsigned char[SizeI*SizeJ];
    float **B = new float*[25];
    for(int i=0; i<25; i++)
    {
        B[i] = new float[SizeI*SizeJ];
    }
    int DataNum = SizeI*SizeJ;
    int FeatureDim = 25; //25 laws filters
    int ClusterNum = 6;

    Mat featureVector(DataNum,FeatureDim,CV_32F);

    ReadImageFromFile("comb.raw",src,SizeI,SizeJ,1);

    for(int i=0; i<25; i++)
    {
        Filter2(&L[i][0],5,src,B[i],SizeI,SizeJ,1);
    }

    unsigned char graylevel[6] = {0, 51, 102, 153, 204, 255};

    int W = 31;
    cout<<"Segmenting..."<<W<<endl;
    fflush(stdout);
    for(int i=0; i<SizeI; i++)
    {
        cout<<"line: "<<i<<endl;
        fflush(stdout);
        for(int j=0; j<SizeJ; j++)
        {
            for(int k=0; k<25; k++)
            {
                featureVector.at<float>(i*SizeJ+j,k) = 0;
                float tmp = 0;
                for(int m=-W/2; m<=W/2; m++)
                {
                    for(int n=-W/2; n<=W/2; n++)
                    {
                        if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                        {
                            tmp += B[k][(i+m)*SizeJ+j+n]*B[k][(i+m)*SizeJ+j+n];
                        }
                    }
                }
                tmp /= W*W;
                featureVector.at<float>(i*SizeJ+j,k) += tmp;
            }
        }
    }

    Mat reducedFeatureVector;
    int maxComponents = 1;
    PCA pca(featureVector,noArray(),CV_PCA_DATA_AS_ROW,maxComponents);
    pca.project(featureVector,reducedFeatureVector);
    cout<<reducedFeatureVector.row(0)<<endl;

    Mat labels;
    int attempts = 2;
    Mat centers;
    // Use kmeans function in Opencv
    kmeans(reducedFeatureVector,ClusterNum,labels,TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001),attempts, KMEANS_RANDOM_CENTERS, centers);

    for(int i=0; i<DataNum; i++)
    {
        output[i] = graylevel[labels.at<int>(i)];
    }

    char filename[10];
    sprintf(filename,"C2_%d.raw",W);
    WriteImageToFile(filename,output,SizeI,SizeJ,1);

    delete src;
    for(int i=0; i<25; i++)
        delete B[i];
    delete B;
    delete output;
    featureVector.release();
    reducedFeatureVector.release();
    labels.release();
    centers.release();

    printf(" C2 ends...\n");
}

int main()
{
    //partA();
    //partB();
    //partC1();
    partC2();
    return 0;
}
