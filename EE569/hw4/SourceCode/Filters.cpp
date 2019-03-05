/**
    Name:	Manan Vyas
	USC ID:	7483-8632-00
	email:	mvyas@usc.edu
*/

#include "Filters.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "IO.h"

void ImageExpandRepeat(int Margin, const unsigned char* ImageData, unsigned char* ImageDataOut, int SizeI, int SizeJ, int BytesPerPixel)
{
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;

    for(int color=0; color<BytesPerPixel; color++)
    {
        // Construct Expanded Image
        // Center
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                ImageDataOut[((i+Margin)*ExpandedJ+j+Margin)*BytesPerPixel+color] = ImageData[BytesPerPixel*(i*SizeJ+j)+color];
            }
        }
        // 4 corners
        for(int i=0; i<Margin; i++)
        {
            for(int j=0; j<Margin; j++)
            {
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+j)+color] = ImageData[color];
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+Margin+SizeJ+j)+color] = ImageData[BytesPerPixel*(SizeJ-1)+color];
                ImageDataOut[BytesPerPixel*((Margin+SizeI+i)*ExpandedJ+j)+color] = ImageData[BytesPerPixel*((SizeI-1)*SizeJ)+color];
                ImageDataOut[BytesPerPixel*((Margin+SizeI+i)*ExpandedJ+Margin+SizeJ+j)+color] = ImageData[BytesPerPixel*(SizeI*SizeJ-1)+color];
            }
        }
        // Top and bottom
        for(int i=0; i<Margin; i++)
        {
            for(int j=Margin; j<ExpandedJ-Margin; j++)
            {
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+j)+color] = ImageData[BytesPerPixel*(j-Margin)+color];
                ImageDataOut[BytesPerPixel*((Margin+SizeI+i)*ExpandedJ+j)+color] = ImageData[BytesPerPixel*((SizeI-1)*SizeJ+j-Margin)+color];
            }
        }
        // Left and right
        for(int i=Margin; i<ExpandedI-Margin; i++)
        {
            for(int j=0; j<Margin; j++)
            {
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+j)+color] = ImageData[BytesPerPixel*(i-Margin)*SizeJ+color];
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+Margin+SizeJ+j)+color] = ImageData[BytesPerPixel*((i-Margin)*SizeJ+SizeJ-1)+color];
            }
        }
    }
}

void ImageExpandSymmetry(int Margin, const unsigned char* ImageData, unsigned char* ImageDataOut, int SizeI, int SizeJ, int BytesPerPixel)
{

}

void ImageExpandBlack(int Margin, const unsigned char* ImageData, unsigned char* ImageDataOut, int SizeI, int SizeJ, int BytesPerPixel)
{
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;

    for(int color=0; color<BytesPerPixel; color++)
    {
        // Construct Expanded Image
        // Center
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                ImageDataOut[((i+Margin)*ExpandedJ+j+Margin)*BytesPerPixel+color] = ImageData[BytesPerPixel*(i*SizeJ+j)+color];
            }
        }
        // 4 corners
        for(int i=0; i<Margin; i++)
        {
            for(int j=0; j<Margin; j++)
            {
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+j)+color] = 0;
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+Margin+SizeJ+j)+color] = 0;
                ImageDataOut[BytesPerPixel*((Margin+SizeI+i)*ExpandedJ+j)+color] = 0;
                ImageDataOut[BytesPerPixel*((Margin+SizeI+i)*ExpandedJ+Margin+SizeJ+j)+color] = 0;
            }
        }
        // Top and bottom
        for(int i=0; i<Margin; i++)
        {
            for(int j=Margin; j<ExpandedJ-Margin; j++)
            {
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+j)+color] = 0;
                ImageDataOut[BytesPerPixel*((Margin+SizeI+i)*ExpandedJ+j)+color] = 0;
            }
        }
        // Left and right
        for(int i=Margin; i<ExpandedI-Margin; i++)
        {
            for(int j=0; j<Margin; j++)
            {
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+j)+color] = 0;
                ImageDataOut[BytesPerPixel*(i*ExpandedJ+Margin+SizeJ+j)+color] = 0;
            }
        }
    }
}

void Filter(const double *FilterData, int FilterSize, const unsigned char* ImageData, unsigned char* ImageDataOut, int SizeI, int SizeJ, int BytesPerPixel)
{
    int Margin = FilterSize/2;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char ExpandedImage[ExpandedI][ExpandedJ][BytesPerPixel];
    // Construct Expanded Image
    ImageExpandRepeat(Margin, ImageData, &ExpandedImage[0][0][0], SizeI, SizeJ, BytesPerPixel);

    for(int color=0; color<BytesPerPixel; color++)
    {
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                double sum = 0;
                for(int m=-Margin; m<=Margin; m++)
                {
                    for(int n=-Margin; n<=Margin; n++)
                    {
                        sum += ExpandedImage[i+Margin+m][j+Margin+n][color] * FilterData[(m+Margin)*FilterSize+n+Margin];
                    }
                }
                ImageDataOut[BytesPerPixel*(i*SizeJ+j)+color] = sum;
                if(sum-(int)sum > 0.5)
                    ImageDataOut[BytesPerPixel*(i*SizeJ+j)+color] ++;
            }
        }
    }
}

int compare (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b );
}
void MedianFilter(char FilterShape, int FilterSize, const unsigned char* ImageData, unsigned char* ImageDataOut, int SizeI, int SizeJ, int BytesPerPixel)
{

    int Margin = FilterSize/2;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char ExpandedImage[ExpandedI][ExpandedJ][BytesPerPixel];
    // Construct Expanded Image
    ImageExpandRepeat(Margin, ImageData, &ExpandedImage[0][0][0], SizeI, SizeJ, BytesPerPixel);

    WriteImageToFile("ExpandedImage.raw",&ExpandedImage[0][0][0],ExpandedI,ExpandedJ,BytesPerPixel);

    int SortedArray[FilterSize*FilterSize];

    for(int color=0; color<BytesPerPixel; color++)
    {
        switch(FilterShape)
        {
        case 'S':
            for(int i=0; i<SizeI; i++)
            {
                for(int j=0; j<SizeJ; j++)
                {
                    for(int m=-Margin; m<=Margin; m++)
                    {
                        for(int n=-Margin; n<=Margin; n++)
                        {
                            SortedArray[(m+Margin)*FilterSize+n+Margin]= ExpandedImage[i+Margin+m][j+Margin+n][color];
                        }
                    }
                    qsort(SortedArray,FilterSize*FilterSize,sizeof(int),compare);
                    ImageDataOut[BytesPerPixel*(i*SizeJ+j)+color] = SortedArray[FilterSize*FilterSize/2];
                }
            }
            break;
        case 'C':
            for(int i=0; i<SizeI; i++)
            {
                for(int j=0; j<SizeJ; j++)
                {
                    int n = 0;
                    for(int m=-Margin; m<=Margin; m++)
                    {
                        SortedArray[n]= ExpandedImage[i+Margin+m][j+Margin][color];
                        n ++;
                    }
                    for(int m=-Margin; m<=Margin; m++)
                    {
                        if(m!=0)
                        {
                            SortedArray[Margin+m+Margin]= ExpandedImage[i+Margin][j+Margin+m][color];
                            n ++;
                        }
                    }
                    qsort(SortedArray,FilterSize+FilterSize-1,sizeof(int),compare);
                    ImageDataOut[BytesPerPixel*(i*SizeJ+j)+color] = SortedArray[(FilterSize+FilterSize-1)/2];
                }
            }
            break;
        }
    }
}

// BilateralFilter for Gray Image
double Function_C(int ita_i,int ita_j,int x_i,int x_j, double sigma_d)
{
    return exp(-((ita_i-x_i)*(ita_i-x_i)+(ita_j-x_j)*(ita_j-x_j))/(2*sigma_d*sigma_d));
}
double Function_S(int ita_gray,int x_gray,double sigma_r)
{
    return exp(-((ita_gray-x_gray)*(ita_gray-x_gray))/(2*sigma_r*sigma_r));
}
void BilateralFilter(int FilterSize,double sigma_d, double sigma_r, const unsigned char* ImageData, unsigned char* ImageDataOut, int SizeI, int SizeJ)
{
    int Margin = FilterSize/2;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char ExpandedGrayImage[ExpandedI][ExpandedJ];

    double DomainFilter[FilterSize][FilterSize];
    double RangeFilter[256];

    // Construct Expanded Image
    ImageExpandRepeat(Margin,ImageData,&ExpandedGrayImage[0][0],SizeI,SizeJ,1);
    ////

    // Construct closeness functions
    for(int m=-Margin; m<=Margin; m++)
    {
        for(int n=-Margin; n<=Margin; n++)
        {
            DomainFilter[m+Margin][n+Margin] = Function_C(m,n,0,0,sigma_d);
        }
    }
    for(int i=0; i<256; i++)
    {
        RangeFilter[i] = Function_S(0,i,sigma_r);
    }

    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            double sum = 0;
            double K = 0;
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    double k = DomainFilter[m+Margin][n+Margin] * RangeFilter[abs(ExpandedGrayImage[i+Margin+m][j+Margin+n]-ImageData[i*SizeJ+j])];
                    K += k;
                    sum += ExpandedGrayImage[i+Margin+m][j+Margin+n] * k;
                }
            }
            ImageDataOut[i*SizeJ+j] = sum / K;
        }
    }
}

// Non Local Means Denoising Filter
void NLMFilter(int NeighborSize, int SearchWindowSize, double h, const unsigned char* ImageData, unsigned char* ImageDataOut, int SizeI, int SizeJ)
{
    int SearchMargin = SearchWindowSize/2;
    int Margin = NeighborSize/2;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char ExpandedGrayImage[ExpandedI][ExpandedJ];

    // Construct Expanded Image
    ImageExpandRepeat(Margin,ImageData,&ExpandedGrayImage[0][0],SizeI,SizeJ,1);
    /////

    printf("Processing line: ");

    for(int i=0; i<SizeI; i++)
    {
        printf("%d ",i);
        fflush(stdout);
        for(int j=0; j<SizeJ; j++)
        {

            double Z = 0;
            double sum = 0;
            for(int m=i-SearchMargin; m<i+SearchMargin; m++)
            {
                if(m<0 || m>=SizeI)
                    continue;
                for(int n=j-SearchMargin; n<j+SearchMargin; n++)
                {
                    if(n<0 || n>=SizeJ)
                        continue;
                    double W = 0;
                    double d = 0;
                    for(int p=-Margin; p<=Margin; p++)
                    {
                        for(int q=-Margin; q<=Margin; q++)
                        {
                            int gray1 = ExpandedGrayImage[i+Margin+p][j+Margin+q];
                            int gray2 = ExpandedGrayImage[m+Margin+p][n+Margin+q];
                            d += (gray1-gray2)*(gray1-gray2);
                        }
                    }
                    W = exp(-d/(h*h));
                    Z += W;
                    sum += W * ExpandedGrayImage[m+Margin][n+Margin];
                }
            }
            ImageDataOut[i*SizeJ+j] = sum / Z;

        }

    }
}

void UniformFilter(int FilterSize,const unsigned char* src, unsigned char* tar, int SizeI, int SizeJ)
{
    int Margin = FilterSize/2;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char ExpandedGrayImage[ExpandedI][ExpandedJ];
    ImageExpandRepeat(Margin,src,&ExpandedGrayImage[0][0],SizeI,SizeJ,1);

    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            double sum = 0;
            double K = 0;
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    double k = 1;
                    K += k;
                    sum += 1.0*ExpandedGrayImage[i+Margin+m][j+Margin+n] * k;
                }
            }
            tar[i*SizeJ+j] = sum / K;
        }
    }
}

void GaussianFilter(int FilterSize,double sigma, const unsigned char* src, unsigned char* tar, int SizeI, int SizeJ)
{
    double PI = 3.14159;
    int Margin = FilterSize/2;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char ExpandedGrayImage[ExpandedI][ExpandedJ];
    ImageExpandRepeat(Margin,src,&ExpandedGrayImage[0][0],SizeI,SizeJ,1);

    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            double sum = 0;
            double K = 0;
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    double k = 1/(2*PI*sigma*sigma)*exp(-(m*m+n*n)/(2*sigma*sigma));
                    K += k;
                    sum += 1.0*ExpandedGrayImage[i+Margin+m][j+Margin+n] * k;
                }
            }
            tar[i*SizeJ+j] = sum / K;
        }
    }
}
void Filter2(const float *FilterData, int FilterSize, const unsigned char* ImageData, float* DataOut, int SizeI, int SizeJ, int BytesPerPixel)
{
    int Margin = FilterSize/2;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char ExpandedImage[ExpandedI][ExpandedJ][BytesPerPixel];
    // Construct Expanded Image
    ImageExpandRepeat(Margin, ImageData, &ExpandedImage[0][0][0], SizeI, SizeJ, BytesPerPixel);

    for(int color=0; color<BytesPerPixel; color++)
    {
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                float sum = 0;
                for(int m=-Margin; m<=Margin; m++)
                {
                    for(int n=-Margin; n<=Margin; n++)
                    {
                        sum += ExpandedImage[i+Margin+m][j+Margin+n][color] * FilterData[(m+Margin)*FilterSize+n+Margin];
                    }
                }
                DataOut[BytesPerPixel*(i*SizeJ+j)+color] = sum;
            }
        }
    }
}
