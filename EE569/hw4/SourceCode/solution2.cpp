/**
EE569 HW #4 Problem 2:OCR
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

#define PI 3.1415926

using namespace std;

typedef struct Point
{
    int i;
    int j;
} Point;

typedef struct Rect
{
    int l;
    int r;
    int t;
    int b;
} Rect;

typedef struct Features
{
    char ch;
    Point center;
    int height;
    int width;
    Rect bbox;
    float NormalArea;
    float NormalPerimeter;
    int EulerNum;
    float Circularity;
    float SpatialMoment_1st_R;
    float SpatialMoment_1st_C;
    float AspectRatio;
    float LeftRatio;
    float RightRatio;
    float TopRatio;
    float BottomRatio;
} Features;

// Part (a) OCR Segmentation and Training
// when abstracting features, object in image should be 255, background is 0
void GetBoundaryBox(unsigned char* img, int SizeI,int SizeJ,Rect* bbox)
{
    bbox->l = SizeJ;
    bbox->r = 0;
    bbox->t = SizeI;
    bbox->b = 0;

    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(img[i*SizeJ+j]==255)
            {
                if(i<bbox->t)
                    bbox->t = i;
                if(i>bbox->b)
                    bbox->b = i;
                if(j<bbox->l)
                    bbox->l = j;
                if(j>bbox->r)
                    bbox->r = j;
            }
        }
    }
}
float NormalArea(unsigned char* img, int SizeI,int SizeJ,int bboxArea)
{
    float area = 0;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(img[i*SizeJ+j] == 255)
                area ++;
        }
    }
    return area/bboxArea;
}
float NormalPerimeter(unsigned char* img, int SizeI, int SizeJ,int bboxPerimeter)
{
    float pmt = 0;
    int neighbori[] = {-1,0,0,1};
    int neighborj[] = {0,-1,1,0};
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(img[i*SizeJ+j] == 255)
            {
                for(int k=0; k<4; k++)
                {
                    if(i+neighbori[k]>=0 && i+neighbori[k]<SizeI && j+neighborj[k]>=0 && j+neighborj[k]<SizeJ)
                    {
                        if(img[(i+neighbori[k])*SizeJ+j+neighborj[k]]==0)
                            pmt ++;
                    }
                }
            }
        }
    }
    return pmt/bboxPerimeter;
}
// 4-connected Euler Number
int EulerNumber(unsigned char* img, int SizeI, int SizeJ)
{
    int EN = 0;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(img[i*SizeJ+j] == 255)
            {
                if(j+1<SizeJ && i+1<SizeI)
                {
                    // Q1:
                    if(img[i*SizeJ+j+1]==0 && img[(i+1)*SizeJ+j]==0 && img[(i+1)*SizeJ+j+1]==0)
                        EN++;
                    //Q3
                    if(img[i*SizeJ+j+1]==255 && img[(i+1)*SizeJ+j]==0 && img[(i+1)*SizeJ+j+1]==255)
                        EN--;

                    if(img[i*SizeJ+j+1]==0 && img[(i+1)*SizeJ+j]==255 && img[(i+1)*SizeJ+j+1]==255)
                        EN--;

                    if(img[i*SizeJ+j+1]==255 && img[(i+1)*SizeJ+j]==255 && img[(i+1)*SizeJ+j+1]==0)
                        EN--;

                    //Qd
                    if(img[i*SizeJ+j+1]==0 && img[(i+1)*SizeJ+j]==0 && img[(i+1)*SizeJ+j+1]==255)
                        EN-=2;
                }


            }
            else
            {
                if(j+1<SizeJ && i+1<SizeI)
                {
                    // Q1
                    if(img[i*SizeJ+j+1]==255 && img[(i+1)*SizeJ+j]==0 && img[(i+1)*SizeJ+j+1]==0)
                        EN++;

                    if(img[i*SizeJ+j+1]==0 && img[(i+1)*SizeJ+j]==0 && img[(i+1)*SizeJ+j+1]==255)
                        EN++;

                    if(img[i*SizeJ+j+1]==0 && img[(i+1)*SizeJ+j]==255 && img[(i+1)*SizeJ+j+1]==0)
                        EN++;

                    // Q3
                    if(img[i*SizeJ+j+1]==255 && img[(i+1)*SizeJ+j]==255 && img[(i+1)*SizeJ+j+1]==255)
                        EN--;

                    // Qd
                    if(img[i*SizeJ+j+1]==255 && img[(i+1)*SizeJ+j]==255 && img[(i+1)*SizeJ+j+1]==0)
                        EN-=2;
                }
            }
        }
    }
    return EN/4;
}
// Circularity can be calculated from area and perimeter.

// Spatial Central Moment
float SpatialMoment(unsigned char* img, int SizeI, int SizeJ, int m, int n)
{
    float M = 0;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            float x = i+0.5;
            float y = j+0.5;
            M += pow(x,m)*pow(y,n)*img[i*SizeJ+j]/255;
        }
    }
    return M;
}
float SpatialCentralMoment(unsigned char* img, int SizeI, int SizeJ, int m, int n, float xc, float yc)
{
    float M = 0;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            float x = i+0.5;
            float y = j+0.5;
            M += pow(x-xc,m)*pow(y-yc,n)*img[i*SizeJ+j]/255;
        }
    }
    return M;
}
void GetRatios(unsigned char* img, int SizeI, int SizeJ,int bboxwidth, int bboxheight, int l, float* lr, int r, float* rr, int t, float* tr, int b, float* br)
{
    int i,j;
    i = t;
    for(j=l;j<=r;j++)
    {
        if(img[i*SizeJ+j]==255)
            *tr += 1;
    }
    *tr = *tr/bboxwidth;
    i = b;
    for(j=l;j<=r;j++)
    {
        if(img[i*SizeJ+j]==255)
            *br += 1;
    }
    *br = *br/bboxwidth;
    j=l;
    for(i=t;i<=b;i++)
    {
        if(img[i*SizeJ+j]==255)
            *lr += 1;
    }
    *lr = *lr/bboxheight;
    j=r;
    for(i=t;i<=b;i++)
    {
        if(img[i*SizeJ+j]==255)
            *rr += 1;
    }
    *rr = *rr/bboxheight;
}

Features* processTrain()
{
    printf("Training begins...\n");
    int SizeI = 300;
    int SizeJ = 600;
    unsigned char* trainImage = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("training-bill.raw",trainImage,SizeI,SizeJ,3);
    unsigned char* imgbuff1 = new unsigned char[SizeI*SizeJ];
    unsigned char* imgbuff2 = new unsigned char[SizeI*SizeJ];
    Color2Gray(trainImage,imgbuff1,SizeI,SizeJ);
    Gray2Binary(imgbuff1,imgbuff2,SizeI,SizeJ,127,'L');
    //WriteImageToFile("test.raw",imgbuff2,SizeI,SizeJ,1);
    int *mapping = new int[SizeI*SizeJ];
    int *labels = new int[100];
    int labelNum = ConnectedComponentLabeling(imgbuff2,mapping,SizeI,SizeJ,SizeI*SizeJ,labels);

    // see segmentation
    /*
    unsigned char* segoutImage = new unsigned char[SizeI*SizeJ*3];
    RGB *rgbs = new RGB[100];
    rgbs[0].R = 255;
    rgbs[0].G = 255;
    rgbs[0].B = 255;
    for(int k=1;k<100;k++)
    {
        rgbs[k].R = 255.0*(rand()%RAND_MAX)/RAND_MAX;
        rgbs[k].G = 255.0*(rand()%RAND_MAX)/RAND_MAX;
        rgbs[k].B = 255.0*(rand()%RAND_MAX)/RAND_MAX;
    }
    for(int i=0;i<SizeI;i++)
    {
        for(int j=0;j<SizeJ;j++)
        {
            segoutImage[3*(i*SizeJ+j)] = rgbs[mapping[i*SizeJ+j]].R;
            segoutImage[3*(i*SizeJ+j)+1] = rgbs[mapping[i*SizeJ+j]].G;
            segoutImage[3*(i*SizeJ+j)+2] = rgbs[mapping[i*SizeJ+j]].B;
        }

    }
    WriteImageToFile("SegoutImage.raw",segoutImage,SizeI,SizeJ,3);
    delete segoutImage;
    ///
    */

    Features *features = new Features[labelNum];
    unsigned char* charImage = new unsigned char[SizeI*SizeJ];
    // Get this alphabet from image that obtained blow
    char characters[] = {'1','2','3','6','8','9','0','4','5','7',
                         '.','T','I','D','B','t','o','a','e','n',
                         'c','u','$','S','O','U','B','A','L'
                        };
    for(int k=0; k<labelNum; k++)
    {
        printf("Processing character %c...\n", characters[k]);
        fflush(stdout);
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                if(mapping[i*SizeJ+j]==labels[k])
                {
                    charImage[i*SizeJ+j] = 255;
                }
                else
                {
                    charImage[i*SizeJ+j] = 0;
                }
            }
        }

        //May use thinning
        // Object is white and background is black before thinning
        // Thinning
        //Thinning(charImage,imgbuff1,SizeI,SizeJ);
        //GrayCopy(imgbuff1,charImage,SizeI,SizeJ);
        //char filename[20];
        //sprintf(filename,"char_thin_%d.raw",k);
        //WriteImageToFile(filename,imgbuff1,SizeI,SizeJ,1);

        // abstract features
        Rect bbox;
        GetBoundaryBox(charImage,SizeI,SizeJ,&bbox);
        /*
        for(int i=bbox.t;i<=bbox.b;i++)
        {
            int j = bbox.l;
            charImage[i*SizeJ+j] = 255;
            j = bbox.r;
            charImage[i*SizeJ+j] = 255;
        }
        for(int j=bbox.l;j<=bbox.r;j++)
        {
            int i = bbox.t;
            charImage[i*SizeJ+j] = 255;
            i = bbox.b;
            charImage[i*SizeJ+j] = 255;
        }
        char filename[10];
        sprintf(filename,"char_%d.raw",k);
        WriteImageToFile(filename,charImage,SizeI,SizeJ,1);
        */
        int bboxArea = (bbox.r-bbox.l+1)*(bbox.b-bbox.t+1);
        int bboxPerimeter = 2*((bbox.r-bbox.l+1)+(bbox.b-bbox.t+1));

        features[k].ch = characters[k];
        features[k].height = bbox.b - bbox.t + 1;
        features[k].width = bbox.r - bbox.l + 1;
        features[k].bbox.b = bbox.b;
        features[k].bbox.t = bbox.t;
        features[k].bbox.l = bbox.l;
        features[k].bbox.r = bbox.r;
        features[k].NormalArea = NormalArea(charImage,SizeI,SizeJ,bboxArea);
        features[k].NormalPerimeter = NormalPerimeter(charImage,SizeI,SizeJ,bboxPerimeter);
        features[k].Circularity = (4*PI*features[k].NormalArea*bboxArea)/(features[k].NormalPerimeter*bboxPerimeter*features[k].NormalPerimeter*bboxPerimeter);
        features[k].EulerNum = EulerNumber(charImage,SizeI,SizeJ);
        features[k].center.i = SpatialMoment(charImage,SizeI,SizeJ,1,0)/SpatialMoment(charImage,SizeI,SizeJ,0,0);
        features[k].center.j = SpatialMoment(charImage,SizeI,SizeJ,0,1)/SpatialMoment(charImage,SizeI,SizeJ,0,0);
        features[k].SpatialMoment_1st_R = SpatialCentralMoment(charImage,SizeI,SizeJ,1,0,features[k].center.i,features[k].center.j);
        features[k].SpatialMoment_1st_C = SpatialCentralMoment(charImage,SizeI,SizeJ,0,1,features[k].center.i,features[k].center.j);
        features[k].AspectRatio = 1.0*features[k].height/features[k].width;
        //GetRatios(charImage,SizeI,SizeJ,features[k].width,features[k].height,features[k].bbox.l,&(features[k].LeftRatio),features[k].bbox.r,&(features[k].RightRatio),
        //          features[k].bbox.t,&(features[k].TopRatio),features[k].bbox.b,&(features[k].BottomRatio));

    }

    delete trainImage;
    delete imgbuff1;
    delete imgbuff2;
    delete mapping;
    delete labels;
    delete charImage;
    printf("Training ends...\n");
    return features;
}

float featureDist(Features *feature1,Features *feature2)
{
    float dist = 0;
    dist += (feature1->AspectRatio-feature2->AspectRatio)*(feature1->AspectRatio-feature2->AspectRatio);
    dist += (feature1->Circularity-feature2->Circularity)*(feature1->Circularity-feature2->Circularity);
    dist += (feature1->NormalArea-feature2->NormalArea)*(feature1->NormalArea-feature2->NormalArea);
    dist += (feature1->NormalPerimeter-feature2->NormalPerimeter)*(feature1->NormalPerimeter-feature2->NormalPerimeter);
    //dist += (feature1->LeftRatio-feature2->LeftRatio)*(feature1->LeftRatio-feature2->LeftRatio);
    //dist += (feature1->RightRatio-feature2->RightRatio)*(feature1->RightRatio-feature2->RightRatio);
    //dist += (feature1->TopRatio-feature2->TopRatio)*(feature1->TopRatio-feature2->TopRatio);
    //dist += (feature1->BottomRatio-feature2->BottomRatio)*(feature1->BottomRatio-feature2->BottomRatio);
    return sqrt(dist);
}

char recognization(unsigned char* img, int SizeI, int SizeJ, Features *featurebase, int basesize)
{
    Features feature;
    Rect bbox;
    GetBoundaryBox(img,SizeI,SizeJ,&bbox);
    int bboxArea = (bbox.r-bbox.l+1)*(bbox.b-bbox.t+1);
    int bboxPerimeter = 2*((bbox.r-bbox.l+1)+(bbox.b-bbox.t+1));

    feature.height = bbox.b - bbox.t + 1;
    feature.width = bbox.r - bbox.l + 1;
    feature.bbox.b = bbox.b;
    feature.bbox.t = bbox.t;
    feature.bbox.l = bbox.l;
    feature.bbox.r = bbox.r;
    feature.NormalArea = NormalArea(img,SizeI,SizeJ,bboxArea);
    feature.NormalPerimeter = NormalPerimeter(img,SizeI,SizeJ,bboxPerimeter);
    feature.Circularity = (4*PI*feature.NormalArea*bboxArea)/(feature.NormalPerimeter*bboxPerimeter*feature.NormalPerimeter*bboxPerimeter);
    feature.EulerNum = EulerNumber(img,SizeI,SizeJ);
    feature.center.i = SpatialMoment(img,SizeI,SizeJ,1,0)/SpatialMoment(img,SizeI,SizeJ,0,0);
    feature.center.j = SpatialMoment(img,SizeI,SizeJ,0,1)/SpatialMoment(img,SizeI,SizeJ,0,0);
    feature.SpatialMoment_1st_R = SpatialCentralMoment(img,SizeI,SizeJ,1,0,feature.center.i,feature.center.j);
    feature.SpatialMoment_1st_C = SpatialCentralMoment(img,SizeI,SizeJ,0,1,feature.center.i,feature.center.j);
    feature.AspectRatio = 1.0*feature.height/feature.width;
    //GetRatios(img,SizeI,SizeJ,feature.width,feature.height,feature.bbox.l,&(feature.LeftRatio),feature.bbox.r,&(feature.RightRatio),
    //              feature.bbox.t,&(feature.TopRatio),feature.bbox.b,&(feature.BottomRatio));


    char r = ' ';

    if(feature.height*feature.width<10)
        r = '.';

    // Decision Tree
    if(feature.EulerNum == -1)
    {
        // 8, B, $
        float minD = 9999;
        for(int b=0; b<basesize; b++)
        {
            if(featurebase[b].EulerNum==-1)
            {
                float tmp = featureDist(&featurebase[b],&feature);
                if(tmp<minD)
                {
                    r = featurebase[b].ch;
                    minD = tmp;
                }
            }
        }
    }
    else
    {
        if(feature.EulerNum == 0)
        {
            // 6,9,0,4,D,o,a,e,O,A
            float minD = 9999;
            for(int b=0; b<basesize; b++)
            {
                if(featurebase[b].EulerNum==0)
                {
                    float tmp = featureDist(&featurebase[b],&feature);
                    if(tmp<minD)
                    {
                        r = featurebase[b].ch;
                        minD = tmp;
                    }
                }
            }
        }
        else
        {
            // 1,2,3,5,7,T,I,t,n,c,u,S,U,L,.
            float minD = 9999;
            for(int b=0; b<basesize; b++)
            {
                if(featurebase[b].EulerNum==1)
                {
                    float tmp = featureDist(&featurebase[b],&feature);
                    if(tmp<minD)
                    {
                        r = featurebase[b].ch;
                        minD = tmp;
                    }
                }
            }
        }
    }
    if(r!=' '&&r!='.')
        printf("Location: (x=%d, y=%d) :",feature.center.j,feature.center.i);
    return r;
}

void RemoveSingleDots(unsigned char* src, unsigned char* tar,int SizeI,int SizeJ,char type)
{
    int obj,bkg;
    if(type=='D')
    {
        obj = 0;
        bkg = 255;
    }
    else
    {
        if(type=='L')
        {
            obj = 255;
            bkg = 0;
        }
        else
            return;
    }
    int neighbori[]= {-1,-1,-1,0,0,1,1,1};
    int neighborj[]= {-1,0,1,-1,1,-1,0,1};
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[i*SizeJ+j] = src[i*SizeJ+j];
            if(src[i*SizeJ+j]==obj)
            {
                int flag = 0;
                for(int k=0; k<8; k++)
                {
                    if(i+neighbori[k]>=0 && i+neighbori[k]<SizeI && j+neighborj[k]>=0 && j+neighborj[k]<SizeJ)
                    if(src[(i+neighbori[k])*SizeJ+j+neighborj[k]] != bkg)
                        flag = 1;
                }
                if(flag == 0)
                    tar[i*SizeJ+j] = bkg;
            }
        }
    }
}
void Reconnect(unsigned char* src, unsigned char* tar,int SizeI,int SizeJ,char type)
{
    int obj,bkg;
    if(type=='D')
    {
        obj = 0;
        bkg = 255;
    }
    else
    {
        if(type=='L')
        {
            obj = 255;
            bkg = 0;
        }
        else
            return;
    }
    for(int i=0;i<SizeI;i++)
    {
        for(int j=0;j<SizeJ;j++)
        {
            tar[i*SizeJ+j] = src[i*SizeJ+j];
            if(src[i*SizeJ+j]==bkg)
            {
                if(i-1>=0 && i+1<SizeI && j-1>=0 && j+1<SizeJ)
                {
                    if(src[(i-1)*SizeJ+j]==obj && src[(i+1)*SizeJ+j]==obj)
                        tar[i*SizeJ+j] = obj;
                    if(src[i*SizeJ+j-1]==obj && src[i*SizeJ+j+1]==obj)
                        tar[i*SizeJ+j] = obj;
                }
            }
        }
    }
}

void processTest(Features *featurebase,int basesize)
{
    printf("Testing begins...\n");
    //int SizeI = 82;
    //int SizeJ = 272;
    //int SizeI = 316;
    //int SizeJ = 901;
    int SizeI = 564;
    int SizeJ = 397;
    unsigned char* testImage = new unsigned char[SizeI*SizeJ*3];
    ReadImageFromFile("restaurant-bill.raw",testImage,SizeI,SizeJ,3);
    unsigned char* imgbuff1 = new unsigned char[SizeI*SizeJ];
    unsigned char* imgbuff2 = new unsigned char[SizeI*SizeJ];
    Color2Gray(testImage,imgbuff1,SizeI,SizeJ);

    //Gray2BinaryAdaptive(imgbuff1,imgbuff2,SizeI,SizeJ,5,5);
    Gray2Binary(imgbuff1,imgbuff2,SizeI,SizeJ,180,'L');
    WriteImageToFile("tmp_bin.raw",imgbuff2,SizeI,SizeJ,1);

    int *mapping = new int[SizeI*SizeJ];
    int *labels = new int[1000];
    int labelNum = ConnectedComponentLabeling(imgbuff2,mapping,SizeI,SizeJ,SizeI*SizeJ,labels);
    unsigned char* charImage = new unsigned char[SizeI*SizeJ];
    /*
    // see segmentation
    unsigned char* segoutImage = new unsigned char[SizeI*SizeJ*3];
    RGB *rgbs = new RGB[100];
    rgbs[0].R = 255;
    rgbs[0].G = 255;
    rgbs[0].B = 255;
    for(int k=1;k<100;k++)
    {
        rgbs[k].R = 255.0*(rand()%RAND_MAX)/RAND_MAX;
        rgbs[k].G = 255.0*(rand()%RAND_MAX)/RAND_MAX;
        rgbs[k].B = 255.0*(rand()%RAND_MAX)/RAND_MAX;
    }
    for(int i=0;i<SizeI;i++)
    {
        for(int j=0;j<SizeJ;j++)
        {
            segoutImage[3*(i*SizeJ+j)] = rgbs[mapping[i*SizeJ+j]].R;
            segoutImage[3*(i*SizeJ+j)+1] = rgbs[mapping[i*SizeJ+j]].G;
            segoutImage[3*(i*SizeJ+j)+2] = rgbs[mapping[i*SizeJ+j]].B;
        }

    }
    WriteImageToFile("SegoutImage.raw",segoutImage,SizeI,SizeJ,3);
    delete segoutImage;
    ///
    */
    for(int k=0; k<labelNum; k++)
    {
        fflush(stdout);
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                if(mapping[i*SizeJ+j]==labels[k])
                {
                    charImage[i*SizeJ+j] = 255;
                }
                else
                {
                    charImage[i*SizeJ+j] = 0;
                }
            }
        }
        /*
        // Object is white and background is black before thinning
        // Thinning
        Thinning(charImage,imgbuff1,SizeI,SizeJ);
        GrayCopy(imgbuff1,charImage,SizeI,SizeJ);
        */
        //Reconnect(charImage,imgbuff1,SizeI,SizeJ,'L');
        GrayCopy(charImage,imgbuff1,SizeI,SizeJ);

        char filename[20];
        sprintf(filename,"char_test_%d.raw",k);
        WriteImageToFile(filename,imgbuff1,SizeI,SizeJ,1);

        char c = recognization(imgbuff1,SizeI,SizeJ,featurebase,basesize);
        if(c!=' ' && c!='.')
        {
            printf("Processing character %d...\n", k);
            printf("%c\n",c);
        }
    }

    delete testImage;
    delete imgbuff1;
    delete imgbuff2;
    delete mapping;
    delete labels;
    delete charImage;
    printf("Testing ends...\n");
}

int main()
{
    Features *featurebase = processTrain();
    int basesize = 29;
    processTest(featurebase,basesize);

    delete featurebase;
    return 0;
}
