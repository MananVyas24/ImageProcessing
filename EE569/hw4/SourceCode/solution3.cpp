/**
EE569 HW #4 Problem 3: Image Segmentation
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

// Part (A):Image segmentation using K-means
// using K-means function provided by opencv
float max3(float a,float b,float c)
{
    if(a>b)
    {
        if(a>c)
            return a;
        else
            return c;
    }
    else
    {
        if(b>c)
            return b;
        else
            return c;
    }
}
void partA()
{
    printf(" A begins...\n");
    int SizeI = 321;
    int SizeJ = 481;
    unsigned char* src = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("seg3.raw",src,SizeI,SizeJ,3);
    Mat featureVector(SizeI*SizeJ,3,CV_32F);
    for(int i=0;i<SizeI;i++)
    {
        for(int j=0;j<SizeJ;j++)
        {
            featureVector.at<float>(i*SizeJ+j,0) = src[3*(i*SizeJ+j)];
            featureVector.at<float>(i*SizeJ+j,1) = src[3*(i*SizeJ+j)+1];
            featureVector.at<float>(i*SizeJ+j,2) = src[3*(i*SizeJ+j)+2];
            //featureVector.at<float>(i*SizeJ+j,3) = 0.299*src[3*(i*SizeJ+j)]+0.587*src[3*(i*SizeJ+j)+1]+0.114*src[3*(i*SizeJ+j)+2];
            //featureVector.at<float>(i*SizeJ+j,3) = max3(src[3*(i*SizeJ+j)],src[3*(i*SizeJ+j)+1],src[3*(i*SizeJ+j)+2]);

        }
    }

    int K = 10;
    Mat labels;
    Mat centers;
    kmeans(featureVector,K,labels,TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001),5, KMEANS_PP_CENTERS, centers);

    unsigned char* tar = new unsigned char[SizeI*SizeJ*3];
    for(int i=0;i<SizeI;i++)
    {
        for(int j=0;j<SizeJ;j++)
        {
            tar[3*(i*SizeJ+j)] = centers.at<float>(labels.at<int>(i*SizeJ+j),0);
            tar[3*(i*SizeJ+j)+1] = centers.at<float>(labels.at<int>(i*SizeJ+j),1);
            tar[3*(i*SizeJ+j)+2] = centers.at<float>(labels.at<int>(i*SizeJ+j),2);
        }
    }
    WriteImageToFile("seg3_kmeans.raw",tar,SizeI,SizeJ,3);


    featureVector.release();
    labels.release();
    centers.release();
    delete src;
    delete tar;
    printf(" A ends...\n");
}

//Part (BC):Image segmentation using mean shift filtering && clustering
void partB()
{
    printf(" B begins...\n");
    int SizeI = 321;
    int SizeJ = 481;
    unsigned char* src = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("seg1.raw",src,SizeI,SizeJ,3);
    unsigned char* tar = new unsigned char[SizeI*SizeJ*3];

//    MeanShiftFiltering(src,tar,SizeI,SizeJ,30,20,10);
//    WriteImageToFile("seg3_msfiltered.raw",tar,SizeI,SizeJ,3);

    MeanShiftSegmentation(src,tar,SizeI,SizeJ,30,30,20,20);
    WriteImageToFile("seg1_mssegmented.raw",tar,SizeI,SizeJ,3);

    delete src;
    delete tar;
    printf(" B ends...\n");
}


int main()
{
    //partA();
    partB();
    return 0;
}

