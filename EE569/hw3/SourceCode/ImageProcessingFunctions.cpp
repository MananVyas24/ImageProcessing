#include "IO.h"
#include "Filters.h"
#include "ImageProcessingFunctions.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

void HistogramEqualization(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ)
{
    int histogram[256];
    int histCumu[256];
    int Transfer[256];

    GetHistogram(src,SizeI,SizeJ,histogram);
    histCumu[0] = histogram[0];
    for(int i=1; i<256; i++)
    {
        histCumu[i] = histCumu[i-1]+histogram[i];
    }
    const int totalN = SizeI*SizeJ;
    // Construct Transfer function
    for(int i=0; i<256; i++)
    {
        Transfer[i] = 255*histCumu[i]/totalN;
    }
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[i*SizeJ+j] = Transfer[src[i*SizeJ+j]];
        }
    }
}

void GrayCopy(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ)
{
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[i*SizeJ+j] = src[i*SizeJ+j];
        }
    }
}

void Color2Gray(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ)
{
    //Y = 0.299R + 0.587G + 0.114B
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[i*SizeJ+j] = (unsigned char)(0.299*src[3*(i*SizeJ+j)]+0.587*src[3*(i*SizeJ+j)+1]+0.114*src[3*(i*SizeJ+j)+2]);
        }
    }
}

void Gray2Binary(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int Threshold,char type)
{
    // type='D':object is dark ; type='L': object is light

    if(type=='D')
    {
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                if(src[i*SizeJ+j]<Threshold)
                    tar[i*SizeJ+j] = 255;
                else
                    tar[i*SizeJ+j] = 0;
            }
        }
    }
    if(type=='L')
    {
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                if(src[i*SizeJ+j]>Threshold)
                    tar[i*SizeJ+j] = 255;
                else
                    tar[i*SizeJ+j] = 0;
            }
        }
    }
}

void RGB2YUV(const unsigned char *src, unsigned char *Y, int *U, int *V,int SizeI, int SizeJ)
{
    //Y = 0.299R + 0.587G + 0.114B
    //U = -0.147R - 0.289G + 0.436B
    //V = 0.615R - 0.515G - 0.100B
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            Y[i*SizeJ+j] = (unsigned char)(0.299*src[3*(i*SizeJ+j)]+0.587*src[3*(i*SizeJ+j)+1]+0.114*src[3*(i*SizeJ+j)+2]);
            U[i*SizeJ+j] = -0.147*src[3*(i*SizeJ+j)]-0.289*src[3*(i*SizeJ+j)+1]+0.436*src[3*(i*SizeJ+j)+2];
            V[i*SizeJ+j] = 0.615*src[3*(i*SizeJ+j)]-0.515*src[3*(i*SizeJ+j)+1]-0.1*src[3*(i*SizeJ+j)+2];
        }
    }
}

void YUV2RGB(const unsigned char *Y, const int *U, const int *V,unsigned char *tar, int SizeI, int SizeJ)
{
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[3*(i*SizeJ+j)] = Y[i*SizeJ+j]+1.14*V[i*SizeJ+j];
            tar[3*(i*SizeJ+j)+1] = Y[i*SizeJ+j]-0.39*U[i*SizeJ+j]-0.58*V[i*SizeJ+j];
            tar[3*(i*SizeJ+j)+2] = Y[i*SizeJ+j]+2.03*U[i*SizeJ+j];
        }
    }
}

void SplitRGB(const unsigned char *src, unsigned char *R,unsigned char *G,unsigned char *B,int SizeI, int SizeJ)
{
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            R[i*SizeJ+j] = src[3*(i*SizeJ+j)];
            G[i*SizeJ+j] = src[3*(i*SizeJ+j)+1];
            B[i*SizeJ+j] = src[3*(i*SizeJ+j)+2];
        }
    }
}

void GrayScaleInverse(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ)
{
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[i*SizeJ+j] = 255 - src[i*SizeJ+j];
        }
    }
}

void EdgeDetectionGradient(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int OperatorType,int thres)
{
    // OperatorType:
    // 0: Pixel difference; 1: Roberts; 2: Sobel;
    int Margin = 1;
    int ExpandedI = SizeI + 2;
    int ExpandedJ = SizeJ + 2;
    unsigned char ExpandedImage[ExpandedI][ExpandedJ];
    ImageExpandRepeat(Margin,src,&ExpandedImage[0][0],SizeI,SizeJ,1);
    int oprow[9];
    int opcol[9];
    switch(OperatorType)
    {
    case 0:
    {
        oprow[0] = 0;
        oprow[1] = 0;
        oprow[2] = 0;
        oprow[3] = -1;
        oprow[4] = 1;
        oprow[5] = 0;
        oprow[6] = 0;
        oprow[7] = 0;
        oprow[8] = 0;

        opcol[0] = 0;
        opcol[1] = 0;
        opcol[2] = 0;
        opcol[3] = 0;
        opcol[4] = 1;
        opcol[5] = 0;
        opcol[6] = 0;
        opcol[7] = -1;
        opcol[8] = 0;
        break;
    }
    case 1:
    {
        oprow[0] = 0;
        oprow[1] = 0;
        oprow[2] = 0;
        oprow[3] = 0;
        oprow[4] = 1;
        oprow[5] = 0;
        oprow[6] = 0;
        oprow[7] = 0;
        oprow[8] = -1;

        opcol[0] = 0;
        opcol[1] = 0;
        opcol[2] = 0;
        opcol[3] = 0;
        opcol[4] = 1;
        opcol[5] = 0;
        opcol[6] = -1;
        opcol[7] = 0;
        opcol[8] = 0;
        break;
    }
    case 2:
    {
        oprow[0] = -1;
        oprow[1] = 0;
        oprow[2] = 1;
        oprow[3] = -2;
        oprow[4] = 0;
        oprow[5] = 2;
        oprow[6] = -1;
        oprow[7] = 0;
        oprow[8] = 1;

        opcol[0] = 1;
        opcol[1] = 2;
        opcol[2] = 1;
        opcol[3] = 0;
        opcol[4] = 0;
        opcol[5] = 0;
        opcol[6] = -1;
        opcol[7] = -2;
        opcol[8] = -1;
        break;
    }
    default:
        ;
    }

    int GradientMap[SizeI][SizeJ];
    unsigned char Edge[SizeI][SizeJ];
    unsigned char GradientImage[SizeI][SizeJ];
    unsigned char GradientImageInv[SizeI][SizeJ];
    int max=-999999,min=999999;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            float r_row = 0;
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    r_row += ExpandedImage[i+Margin+m][j+Margin+n] * oprow[(m+1)*3+n+1];
                }
            }
            float r_col = 0;
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    r_col += ExpandedImage[i+Margin+m][j+Margin+n] * opcol[(m+1)*3+n+1];
                }
            }
            GradientMap[i][j] = sqrt(r_row*r_row + r_col*r_col);
            if(GradientMap[i][j]>max)
                max = GradientMap[i][j];
            if(GradientMap[i][j]<min)
                min = GradientMap[i][j];
            if(GradientMap[i][j] > thres)
            {
                Edge[i][j] = 255;
            }
            else
            {
                Edge[i][j] = 0;
            }
        }
    }

    WriteImageToFile("EdgeImage_Grad.raw",&Edge[0][0],SizeI,SizeJ,1);
    WriteDataToFile("GradientData_Grad.dat",&GradientMap[0][0],SizeI*SizeJ);

    memset(tar,255,SizeI*SizeJ);
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            GradientImage[i][j] = 255.0*(GradientMap[i][j]-min)/max;
        }
    }
    GrayScaleInverse(&GradientImage[0][0],&GradientImageInv[0][0],SizeI,SizeJ);
    WriteImageToFile("GradientImage_1st_out.raw",&GradientImageInv[0][0],SizeI,SizeJ,1);
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(Edge[i][j]==255)
                tar[i*SizeJ+j] = GradientImageInv[i][j];
        }
    }
}

void EdgeDetectionLaplacian(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int thres)
{
    int Margin = 1;
    int ExpandedI = SizeI + 2;
    int ExpandedJ = SizeJ + 2;
    unsigned char *ExpandedImage = new unsigned char[ExpandedI*ExpandedJ];
    ImageExpandRepeat(Margin,src,ExpandedImage,SizeI,SizeJ,1);
    int op[9];
    op[0] = 0;
    op[1] = -1;
    op[2] = 0;
    op[3] = -1;
    op[4] = 4;
    op[5] = -1;
    op[6] = 0;
    op[7] = -1;
    op[8] = 0;

    int *GradientMap = new int[SizeI*SizeJ];
    unsigned char *Edge = new unsigned char[SizeI*SizeJ];
    unsigned char *GradientImage = new unsigned char[SizeI*SizeJ];
    unsigned char *GradientImageInv = new unsigned char[SizeI*SizeJ];
    int *threeLevelImage = new int[ExpandedI*ExpandedJ];

    int max=-999999,min=999999;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            float r = 0;
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    r += ExpandedImage[(i+Margin+m)*ExpandedJ+j+Margin+n] * op[(m+1)*3+n+1];
                }
            }
            GradientMap[i*SizeJ+j] = r;
            if(GradientMap[i*SizeJ+j]>max)
                max = GradientMap[i*SizeJ+j];
            if(GradientMap[i*SizeJ+j]<min)
                min = GradientMap[i*SizeJ+j];
            if(r < -thres)
            {
                threeLevelImage[(i+Margin)*ExpandedJ+j+Margin] = -1;
            }
            else
            {
                if(r > thres)
                {
                    threeLevelImage[(i+Margin)*ExpandedJ+j+Margin] = 1;
                }
                else
                    threeLevelImage[(i+Margin)*ExpandedJ+j+Margin] = 0;
            }
        }
    }
    // Pattern matching
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            int c[9];
            c[0] = threeLevelImage[(Margin+i-1)*ExpandedJ+Margin+j-1];
            c[1] = threeLevelImage[(Margin+i-1)*ExpandedJ+Margin+j];
            c[2] = threeLevelImage[(Margin+i-1)*ExpandedJ+Margin+j+1];
            c[3] = threeLevelImage[(Margin+i)*ExpandedJ+Margin+j-1];
            c[4] = threeLevelImage[(Margin+i)*ExpandedJ+Margin+j];
            c[5] = threeLevelImage[(Margin+i)*ExpandedJ+Margin+j+1];
            c[6] = threeLevelImage[(Margin+i+1)*ExpandedJ+Margin+j-1];
            c[7] = threeLevelImage[(Margin+i+1)*ExpandedJ+Margin+j];
            c[8] = threeLevelImage[(Margin+i+1)*ExpandedJ+Margin+j+1];
            // Actually if we don't only check if zero , we will get a better result...
            if(threeLevelImage[(Margin+i)*ExpandedJ+Margin+j] == 0)
            {
                /*
                                // \ | /
                                // - x -
                                // / | \
                            */
                if((c[0]+c[8]==0 && c[0] != 0) || (c[1]+c[7]==0 && c[1] != 0) || (c[2]+c[6]==0 && c[2] != 0) || (c[3]+c[5]==0 && c[3] != 0))
                {
                    Edge[i*SizeJ+j] = 255;
                }
                else
                    Edge[i*SizeJ+j] = 0;
            }
            else
            {
                if(threeLevelImage[(Margin+i)*ExpandedJ+Margin+j] == 1)
                {
                    if((c[0] == -1) || (c[1] == -1) || (c[2] == -1) || (c[3] == -1) || (c[4] == -1) ||
                            (c[5] == -1) || (c[6] == -1) || (c[7] == -1) || (c[8] == -1))
                    {
                        Edge[i*SizeJ+j] = 255;
                    }
                    else
                        Edge[i*SizeJ+j] = 0;
                }
                else
                {
                    Edge[i*SizeJ+j] = 0;
                }
            }
        }
    }

    WriteImageToFile("EdgeImage_Lapl.raw",Edge,SizeI,SizeJ,1);
    WriteDataToFile("GradientData_Lapl.dat",GradientMap,SizeI*SizeJ);

    memset(tar,255,SizeI*SizeJ);
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            //double tmp = 255.0*(GradientMap[i*SizeJ+j]-min)/(max-min);
            double tmp = abs(GradientMap[i*SizeJ+j]);
            GradientImage[i*SizeJ+j] = tmp;
            if(Edge[i*SizeJ+j]==255)
                tar[i*SizeJ+j] = GradientImage[i*SizeJ+j];
        }
    }
    GrayScaleInverse(GradientImage,GradientImageInv,SizeI,SizeJ);
    WriteImageToFile("GradientImage_2nd_out.raw",GradientImageInv,SizeI,SizeJ,1);
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(Edge[i*SizeJ+j]==255)
                tar[i*SizeJ+j] = GradientImageInv[i*SizeJ+j];
        }
    }

    delete ExpandedImage;
    delete GradientMap;
    delete GradientImage;
    delete threeLevelImage;
    delete Edge;
    delete GradientImageInv;
}

void EdgeDetectionDOG(double sigma, double k, int FilterSize, const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int thres)
{
    unsigned char *TmpImage1 = new unsigned char[SizeI*SizeJ];
    unsigned char *TmpImage2 = new unsigned char[SizeI*SizeJ];
    int *DifferenceMap = new int[SizeI*SizeJ];
    unsigned char *DifferenceImage = new unsigned char[SizeI*SizeJ];
    unsigned char *DifferenceImageInv = new unsigned char[SizeI*SizeJ];
    unsigned char *Edge = new unsigned char[SizeI*SizeJ];
    GaussianFilter(FilterSize,sigma,src,TmpImage1,SizeI,SizeJ);
    GaussianFilter(FilterSize,k*sigma,src,TmpImage2,SizeI,SizeJ);
    int min = 999999;
    int max = -999999;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            DifferenceMap[i*SizeJ+j] = (int)TmpImage1[i*SizeJ+j] - TmpImage2[i*SizeJ+j];
            if(DifferenceMap[i*SizeJ+j]>max)
                max = DifferenceMap[i*SizeJ+j];
            if(DifferenceMap[i*SizeJ+j]<min)
                min = DifferenceMap[i*SizeJ+j];
            if(DifferenceMap[i*SizeJ+j]>thres)
                Edge[i*SizeJ+j] = 255;
            else
                Edge[i*SizeJ+j] = 0;

        }
    }
    WriteDataToFile("GradientData_Dog.dat",DifferenceMap,SizeI*SizeJ);

    memset(tar,255,SizeI*SizeJ);
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            double tmp = 255.0*(DifferenceMap[i*SizeJ+j]-min)/(max-min);
            DifferenceImage[i*SizeJ+j] = tmp;
            if(Edge[i*SizeJ+j]==255)
                tar[i*SizeJ+j] = DifferenceImage[i*SizeJ+j];
        }
    }
    GrayScaleInverse(DifferenceImage,DifferenceImageInv,SizeI,SizeJ);
    WriteImageToFile("GradientImage_dog_out.raw",DifferenceImageInv,SizeI,SizeJ,1);
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(Edge[i*SizeJ+j]==255)
                tar[i*SizeJ+j] = DifferenceImageInv[i*SizeJ+j];
        }
    }
    delete TmpImage1;
    delete TmpImage2;
    delete DifferenceMap;
    delete DifferenceImage;
    delete DifferenceImageInv;
    delete Edge;
}

int checkConditionalPattern(int pattern, char type)
{
    int Sb = 0, Se = 57;
    int Tb = 16, Te = 61;
    int Kb = 26, Ke = 65;
    int CPatterns[]= {0b001010000,0b100010000,0b000010100,0b000010001, //S
                      0b000011000,0b010010000,0b000110000,0b000010010, //S
                      0b001011000,0b011010000,0b110010000,0b100110000, //S
                      0b000110100,0b000010110,0b000010011,0b000011001, //S
                      0b110011000,0b010011001,0b011110000,0b001011010, //ST
                      0b011011000,0b110110000,0b000110110,0b000011011, //ST
                      0b110011001,0b011110100, //ST

                      0b011011011,0b111111000,0b110110110,0b000111111, //STK
                      0b111011011,0b011011111,0b111111100,0b111111001, //STK
                      0b111110110,0b110110111,0b100111111,0b001111111, //STK
                      0b111011111,0b111111101,0b111110111,0b101111111, //STK
                      0b001011001,0b111010000,0b100110100,0b000010111, //STK
                      0b111011000,0b011011001,0b111110000,0b110110100, //STK
                      0b100110110,0b000110111,0b000011111,0b001011011, //STK
                      0b111011001,0b111110100,0b100110111,0b001011111, //STK
                      0b010011000,0b010110000,0b000110010,0b000011010, //TK
                      0b111111011,0b111111110,0b110111111,0b011111111
                     };// K
    switch (type)
    {
    case 'S':
        for(int i=Sb; i<=Se; i++)
        {
            if(pattern == CPatterns[i])
                return 1;
        }
        break;
    case 'T':
        for(int i=Tb; i<=Te; i++)
        {
            if(pattern == CPatterns[i])
                return 1;
        }
        break;
    case 'K':
        for(int i=Kb; i<=Ke; i++)
        {
            if(pattern == CPatterns[i])
                return 1;
        }
        break;
    }

    return 0;
}

int checkUnconditionalPattern(int pattern, char type)
{
    int UCPatternsSTWithoutD[]= {0b001010000,0b100010000, // Spur
                                 0b000010010,0b000011000, // Single 4-connection
                                 0b001011000,0b011010000,0b110010000,0b100110000, //LCluster
                                 0b000110100,0b000010110,0b000010011,0b000011001,
                                 0b011110000,0b110011000,0b010011001,0b001011010, // 4-connected offset
                                 0b011011100,0b001011100,0b011010100, // Spur corner Cluster
                                 0b110110001,0b100110001,0b110010001,
                                 0b001110110,0b001110100,0b001010110,
                                 0b100011011,0b100011001,0b100010011
                                };
    int UCPatternsSTWithD[]= {0b110110000, // Corner Cluster
                              0b010111000,// TeeBranch
                              0b010111000,
                              0b000111010,
                              0b000111010,
                              0b010110010,
                              0b010110010,
                              0b010011010,
                              0b010011010,
                              0b101010001,0b101010010,0b101010011,0b101010100,0b101010101,0b101010110,0b101010111,// VeeBranch
                              0b100010101,0b100011100,0b100011101,0b101010100,0b101010101,0b101011100,0b101011101,
                              0b001010101,0b010010101,0b011010101,0b100010101,0b101010101,0b110010101,0b111010101,
                              0b001010101,0b001110001,0b001110101,0b101010001,0b101010101,0b101110001,0b101110101,
                              0b010011100,// DiagonalBranch
                              0b010110001,
                              0b001110010,
                              0b100011010
                             };
    int UCPatternsSTDMask[]= {0b110110000, // Corner Cluster
                              0b011111011,// TeeBranch
                              0b110111110,
                              0b110111110,
                              0b011111011,
                              0b010111111,
                              0b111111010,
                              0b111111010,
                              0b010111111,
                              0b101010111,0b101010111,0b101010111,0b101010111,0b101010111,0b101010111,0b101010111,// VeeBranch
                              0b101011111,0b101011111,0b101011111,0b101011111,0b101011111,0b101011111,0b101011111,
                              0b111010101,0b111010101,0b111010101,0b111010101,0b111010101,0b111010101,0b111010101,
                              0b101110101,0b101110101,0b101110101,0b101110101,0b101110101,0b101110101,0b101110101,
                              0b011111110,// DiagonalBranch
                              0b110111011,
                              0b011111110,
                              0b110111011
                             };
    int UCPatternsKWithoutD[]= {0b000010001,0b000010100,0b001010000,0b100010000, // Spur
                                0b000010010,0b000011000,0b000110000,0b010010000, // Single 4-connection
                                0b010011000,0b010110000,0b000011010,0b000110010// LCorner
                               };
    int UCPatternsKWithD[]=
    {
        0b110110000,// Corner Cluster
        0b000011011,
        0b010111000,// TeeBranch
        0b010110010,
        0b000111010,
        0b010011010,
        0b101010001,0b101010010,0b101010011,0b101010100,0b101010101,0b101010110,0b101010111,// VeeBranch
        0b100010101,0b100011100,0b100011101,0b101010100,0b101010101,0b101011100,0b101011101,
        0b001010101,0b010010101,0b011010101,0b100010101,0b101010101,0b110010101,0b111010101,
        0b001010101,0b001110001,0b001110101,0b101010001,0b101010101,0b101110001,0b101110101,
        0b010011100,// DiagonalBranch
        0b010110001,
        0b001110010,
        0b100011010,
    };
    int UCPatternsKDMask[]= {0b110110000, // Corner Cluster
                             0b000011011,
                             0b010111000,// TeeBranch
                             0b010110010,
                             0b000111010,
                             0b010011010,
                             0b101010111,0b101010111,0b101010111,0b101010111,0b101010111,0b101010111,0b101010111,// VeeBranch
                             0b101011111,0b101011111,0b101011111,0b101011111,0b101011111,0b101011111,0b101011111,
                             0b111010101,0b111010101,0b111010101,0b111010101,0b111010101,0b111010101,0b111010101,
                             0b101110101,0b101110101,0b101110101,0b101110101,0b101110101,0b101110101,0b101110101,
                             0b011111110,// DiagonalBranch
                             0b110111011,
                             0b011111110,
                             0b110111011
                            };
    int LenSTWithoutD = sizeof(UCPatternsSTWithoutD)/sizeof(int);
    int LenSTWithD = sizeof(UCPatternsSTWithD)/sizeof(int);
    int LenKWithoutD = sizeof(UCPatternsKWithoutD)/sizeof(int);
    int LenKWithD = sizeof(UCPatternsKWithD)/sizeof(int);
    switch (type)
    {
    case 'S':
    case 'T':
        for(int i=0; i<LenSTWithoutD; i++)
        {
            if(pattern == UCPatternsSTWithoutD[i])
                return 1;
        }
        for(int i=0; i<LenSTWithD; i++)
        {
            int tmp = pattern & UCPatternsSTDMask[i];
            if(tmp == UCPatternsSTWithD[i])
                return 1;
        }
        break;
    case 'K':
        for(int i=0; i<LenKWithoutD; i++)
        {
            if(pattern == UCPatternsKWithoutD[i])
                return 1;
        }
        for(int i=0; i<LenKWithD; i++)
        {
            int tmp = pattern & UCPatternsKDMask[i];
            if(tmp == UCPatternsKWithD[i])
                return 1;
        }
        break;
    }
    return 0;
}

int MorphOperation(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ,char type)
{
    // S: Shrinking, T: Thinning, K: Skeletonizing
    int remainPixel = 0;
    int Margin = 1;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char *ExpandedImage = new unsigned char[ExpandedI*ExpandedJ];

    ImageExpandBlack(Margin,src,ExpandedImage,SizeI,SizeJ,1);

    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            ExpandedImage[(i+Margin)*ExpandedJ+j+Margin] = src[i*SizeJ+j];
        }
    }
    unsigned char*ToBeRemoved = new unsigned char[SizeI*SizeJ];
    unsigned char *ExpandedToBeRemoved = new unsigned char[ExpandedI*ExpandedJ];
    unsigned char*FinalRemove = new unsigned char[SizeI*SizeJ];
    memset(ToBeRemoved,0,SizeI*SizeJ);
    memset(FinalRemove,0,SizeI*SizeJ);
    // Stage 1
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(ExpandedImage[(i+Margin)*ExpandedJ+j+Margin] == 255)
            {
                int pattern = 0;
                for(int m=-1; m<=1; m++)
                {
                    for(int n=-1; n<=1; n++)
                    {
                        pattern = pattern<<1;
                        if(ExpandedImage[(i+Margin+m)*ExpandedJ+j+Margin+n]==255)
                            pattern += 1;
                    }
                }
                int ifhit = checkConditionalPattern(pattern,type);
                if(ifhit==1)
                    ToBeRemoved[i*SizeJ+j] = 255;
                else
                    ToBeRemoved[i*SizeJ+j] = 0;
            }
        }
    }
    WriteImageToFile("TobeRemoved.raw",ToBeRemoved,SizeI,SizeJ,1);
    ImageExpandBlack(Margin,ToBeRemoved,ExpandedToBeRemoved,SizeI,SizeJ,1);
    // Stage 2
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(ToBeRemoved[i*SizeJ+j]==255)
            {
                int pattern = 0;
                for(int m=-1; m<=1; m++)
                {
                    for(int n=-1; n<=1; n++)
                    {
                        pattern = pattern<<1;
                        if(ExpandedToBeRemoved[(i+Margin+m)*ExpandedJ+j+Margin+n]==255)
                            pattern += 1;
                    }
                }
                int ifhit = checkUnconditionalPattern(pattern,type);
                if(ifhit==1)
                    FinalRemove[i*SizeJ+j] = 0;
                else
                    FinalRemove[i*SizeJ+j] = 255;
            }
        }
    }
    WriteImageToFile("FinalRemove.raw",FinalRemove,SizeI,SizeJ,1);
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(FinalRemove[i*SizeJ+j]==255)
                tar[i*SizeJ+j] = 0;
            else
                tar[i*SizeJ+j] = src[i*SizeJ+j];
            if(tar[i*SizeJ+j] == 255)
                remainPixel ++;
        }
    }

    delete ExpandedImage;
    delete ToBeRemoved;
    delete ExpandedToBeRemoved;
    delete FinalRemove;

    return remainPixel;
}

void Thinning(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ)
{
    unsigned char* ImageBuff1 = new unsigned char[SizeI*SizeJ];
    unsigned char* ImageBuff2 = new unsigned char[SizeI*SizeJ];
    char opertion = 'T';

    int remain_pre = 0;
    int remain = MorphOperation(src,ImageBuff1,SizeI,SizeJ,opertion);
    int outbuff = 2;

    while(remain != remain_pre)
    {
        remain_pre = remain;
        if(outbuff == 1)
        {
            remain = MorphOperation(ImageBuff2,ImageBuff1,SizeI,SizeJ,opertion);
            outbuff = 2;
        }
        else
        {
            remain = MorphOperation(ImageBuff1,ImageBuff2,SizeI,SizeJ,opertion);
            outbuff = 1;
        }
    }
    if(outbuff == 1)
    {
        GrayCopy(ImageBuff2,tar,SizeI,SizeJ);
    }
    else
    {
        GrayCopy(ImageBuff1,tar,SizeI,SizeJ);
    }
    delete ImageBuff1;
    delete ImageBuff2;
}

void Dilation4(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ)
{
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(src[i*SizeJ+j] == 255)
            {
                tar[i*SizeJ+j] = 255;
                if(i-1>=0)
                    tar[(i-1)*SizeJ+j] = 255;
                if(i+1<SizeI)
                    tar[(i+1)*SizeJ+j] = 255;
                if(j-1>=0)
                    tar[i*SizeJ+j-1] = 255;
                if(j+1<SizeJ)
                    tar[i*SizeJ+j+1] = 255;
            }
        }
    }
}

void Erode4(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ)
{
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[i*SizeJ+j] = 0;
            if(src[i*SizeJ+j] == 255)
            {
                if(i-1>=0 && src[(i-1)*SizeJ+j] == 255)
                {
                    if(i+1<SizeI && src[(i+1)*SizeJ+j] == 255)
                    {
                        if(j-1>=0 && src[i*SizeJ+j-1] == 255)
                        {
                            if(j+1<SizeJ && src[i*SizeJ+j+1] == 255)
                            {
                                tar[i*SizeJ+j] = 255;
                            }
                        }
                    }
                }
            }
        }
    }
}

void FillSingleHole(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ)
{
    int Margin = 1;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char *ExpandedImage = new unsigned char[ExpandedI*ExpandedJ];
    ImageExpandBlack(Margin,src,ExpandedImage,SizeI,SizeJ,1);

    int neighbori[]= {-1,-1,-1,0,0,1,1,1};
    int neighborj[]= {-1,0,1,-1,1,-1,0,1};
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            tar[i*SizeJ+j] = src[i*SizeJ+j];
            if(src[i*SizeJ+j]==0)
            {
                int flag = 0;
                for(int k=0; k<8; k++)
                {
                    if(ExpandedImage[(i+Margin+neighbori[k])*ExpandedJ+j+Margin+neighborj[k]] != 255)
                        flag = 1;
                }
                if(flag == 0)
                    tar[i*SizeJ+j] = 255;
            }
        }
    }

    delete ExpandedImage;
}

int CountConnectedObjects(const unsigned char *src, int SizeI, int SizeJ, int Thres)
{
    int count = 0;
    unsigned char *checker = new unsigned char[SizeI*SizeJ];
    memset(checker,0,SizeI*SizeJ);
    int *m_stack_i = new int[SizeI*SizeJ];
    int *m_stack_j = new int[SizeI*SizeJ];
    int m_stack_p = 0;
    int current_i,current_j;
    int neighbori[]= {-1,-1,-1,0,0,1,1,1};
    int neighborj[]= {-1,0,1,-1,1,-1,0,1};
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(src[i*SizeJ+j]==255 && checker[i*SizeJ+j]==0)
            {
                int objSize = 0;
                m_stack_i[m_stack_p] = i;
                m_stack_j[m_stack_p] = j;
                checker[i*SizeJ+j] = 255;
                m_stack_p ++;
                while(m_stack_p>0)
                {
                    current_i = m_stack_i[m_stack_p-1];
                    current_j = m_stack_j[m_stack_p-1];
                    m_stack_p--;
                    objSize ++;
                    for(int k=0; k<8; k++)
                    {
                        if(current_i+neighbori[k]>=0 && current_i+neighbori[k]<SizeI && current_j+neighborj[k]>=0 && current_j+neighborj[k]<SizeJ)
                        {
                            if(src[(current_i+neighbori[k])*SizeJ+current_j+neighborj[k]]==255 && checker[(current_i+neighbori[k])*SizeJ+current_j+neighborj[k]]==0)
                            {
                                m_stack_i[m_stack_p] = current_i+neighbori[k];
                                m_stack_j[m_stack_p] = current_j+neighborj[k];
                                checker[(current_i+neighbori[k])*SizeJ+current_j+neighborj[k]] = 255;
                                m_stack_p++;
                            }
                        }
                    }
                }
                if(objSize >= Thres)
                    count ++;
            }
        }
    }

    delete checker;
    delete m_stack_i;
    delete m_stack_j;

    return count;
}

int CountSingleDots(const unsigned char *src, int SizeI, int SizeJ)
{
    int count = 0;
    int Margin = 1;
    int ExpandedI = SizeI + 2*Margin;
    int ExpandedJ = SizeJ + 2*Margin;
    unsigned char *ExpandedImage = new unsigned char[ExpandedI*ExpandedJ];
    ImageExpandBlack(Margin,src,ExpandedImage,SizeI,SizeJ,1);

    int neighbori[]= {-1,-1,-1,0,0,1,1,1};
    int neighborj[]= {-1,0,1,-1,1,-1,0,1};
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            if(src[i*SizeJ+j]==255)
            {
                int flag = 0;
                for(int k=0; k<8; k++)
                {
                    if(ExpandedImage[(i+Margin+neighbori[k])*ExpandedJ+j+Margin+neighborj[k]] != 0)
                        flag = 1;
                }
                if(flag == 0)
                    count ++;
            }
        }
    }

    delete ExpandedImage;

    return count;
}

int CountPattern(const unsigned char *src, int SizeI, int SizeJ, int pattern,int patternSize)
{
    int count = 0;
    int Margin = patternSize/2;
    for(int i=Margin; i<SizeI-Margin; i++)
    {
        for(int j=Margin; j<SizeJ-Margin; j++)
        {
            int thispattern = 0;
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    thispattern = thispattern<<1;
                    if(src[(i+m)*SizeJ+j+n]==255)
                        thispattern += 1;
                }
            }
            if(thispattern == pattern)
                count ++;
        }
    }
    return count;
}

double randGenerator(char type)
{
    double PI = 3.14159265;
    double r1 = 1.0*(rand()%RAND_MAX)/RAND_MAX;
    double r2 = 1.0*(rand()%RAND_MAX)/RAND_MAX;
    switch (type)
    {
    case 'U':
    {
        return r1;
        break;
    }
    case 'G': // Use  Box–Muller method to generate Gaussian r.v
    {
        return sqrt(-2*log(r1))*cos(2*PI*r2);
        break;
    }
    default:
        return 0;
    }
}

void Dithering(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ,char type,int *indexM,int N)
{
    // indexM is of length N*N
    switch(type)
    {
    case 0: // Fixed Threshold
    {
        int T = 127;
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                if(src[i*SizeJ+j]<T)
                    tar[i*SizeJ+j] = 0;
                else
                    tar[i*SizeJ+j] = 255;
            }
        }
        break;
    }
    case 1: // Random Threshold
    {
        int T;
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                T = randGenerator('U')*256;
                //T = 50*randGenerator('G')+127;
                if(src[i*SizeJ+j]<T)
                    tar[i*SizeJ+j] = 0;
                else
                    tar[i*SizeJ+j] = 255;
            }
        }
        break;
    }
    case 2: //Dithering Matrix
    {
        double T;
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                T = 255.0*(indexM[N*(i%N)+(j%N)]+0.5)/(N*N);
                if(src[i*SizeJ+j]<T)
                    tar[i*SizeJ+j] = 0;
                else
                    tar[i*SizeJ+j] = 255;
            }
        }
        break;
    }

    }
}

void ErrorDiffusion(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ,int Method,int Thres)
{
    double *temp = new double[SizeI*SizeJ];
    memset(temp,0,sizeof(temp));
    switch(Method)
    {
    case 1: // Floyd-Steinberg
    {
        double DiffusionM0[3][3] = {{0,0,0},{0,0,7.0/16},{3.0/16,5.0/16,1.0/16}};
        double DiffusionM1[3][3] = {{0,0,0},{7.0/16,0,0},{1.0/16,5.0/16,3.0/16}};
        int dir = 0;
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                temp[i*SizeJ+j] = src[i*SizeJ+j];
            }
        }
        for(int i=0; i<SizeI; i++)
        {
            dir = i%2;
            if(dir==0) // left to right
            {
                for(int j=0; j<SizeJ; j++)
                {
                    if(temp[i*SizeJ+j]<Thres)
                        tar[i*SizeJ+j] = 0;
                    else
                        tar[i*SizeJ+j] = 255;
                    double err = temp[i*SizeJ+j] - tar[i*SizeJ+j];
                    for(int m=-1; m<=1; m++)
                    {
                        for(int n=-1; n<=1; n++)
                        {
                            if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                            {
                                temp[(i+m)*SizeJ+j+n] += err*DiffusionM0[m+1][n+1];
                            }
                        }
                    }

                }

            }
            else
            {
                for(int j=SizeJ-1; j>=0; j--)
                {
                    if(temp[i*SizeJ+j]<Thres)
                        tar[i*SizeJ+j] = 0;
                    else
                        tar[i*SizeJ+j] = 255;
                    double err = temp[i*SizeJ+j] - tar[i*SizeJ+j];
                    for(int m=-1; m<=1; m++)
                    {
                        for(int n=-1; n<=1; n++)
                        {
                            if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                            {
                                temp[(i+m)*SizeJ+j+n] += err*DiffusionM1[m+1][n+1];
                            }
                        }
                    }

                }
            }
        }
        break;
    }
    case 2: //JJN
    {
        double DiffusionM0[5][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,7.0/48,5.0/48},{3.0/48,5.0/48,7.0/48,5.0/48,3.0/48},{1.0/48,3.0/48,5.0/48,3.0/48,1.0/48}};
        double DiffusionM1[5][5] = {{0,0,0,0,0},{0,0,0,0,0},{5.0/48,7.0/48,0,0,0},{3.0/48,5.0/48,7.0/48,5.0/48,3.0/48},{1.0/48,3.0/48,5.0/48,3.0/48,1.0/48}};
        int dir = 0;
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                temp[i*SizeJ+j] = src[i*SizeJ+j];
            }
        }
        for(int i=0; i<SizeI; i++)
        {
            dir = i%2;
            if(dir==0) // left to right
            {
                for(int j=0; j<SizeJ; j++)
                {
                    if(temp[i*SizeJ+j]<Thres)
                        tar[i*SizeJ+j] = 0;
                    else
                        tar[i*SizeJ+j] = 255;
                    double err = temp[i*SizeJ+j] - tar[i*SizeJ+j];
                    for(int m=-2; m<=2; m++)
                    {
                        for(int n=-2; n<=2; n++)
                        {
                            if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                            {
                                temp[(i+m)*SizeJ+j+n] += err*DiffusionM0[m+2][n+2];
                            }
                        }
                    }
                }

            }
            else
            {
                for(int j=SizeJ-1; j>=0; j--)
                {
                    if(temp[i*SizeJ+j]<Thres)
                        tar[i*SizeJ+j] = 0;
                    else
                        tar[i*SizeJ+j] = 255;
                    double err = temp[i*SizeJ+j] - tar[i*SizeJ+j];
                    for(int m=-2; m<=2; m++)
                    {
                        for(int n=-2; n<=2; n++)
                        {
                            if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                            {
                                temp[(i+m)*SizeJ+j+n] += err*DiffusionM1[m+2][n+2];
                            }
                        }
                    }
                }
            }
        }
        break;
    }
    case 3: //Stucki
    {
        double DiffusionM0[5][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,8.0/42,4.0/42},{2.0/42,4.0/42,8.0/42,4.0/42,2.0/42},{1.0/42,2.0/42,4.0/42,2.0/42,1.0/42}};
        double DiffusionM1[5][5] = {{0,0,0,0,0},{0,0,0,0,0},{4.0/42,8.0/42,0,0,0},{2.0/42,4.0/42,8.0/42,4.0/42,2.0/42},{1.0/42,2.0/42,4.0/42,2.0/42,1.0/42}};
        int dir = 0;
        for(int i=0; i<SizeI; i++)
        {
            for(int j=0; j<SizeJ; j++)
            {
                temp[i*SizeJ+j] = src[i*SizeJ+j];
            }
        }
        for(int i=0; i<SizeI; i++)
        {
            dir = i%2;
            if(dir==0) // left to right
            {
                for(int j=0; j<SizeJ; j++)
                {
                    if(temp[i*SizeJ+j]<Thres)
                        tar[i*SizeJ+j] = 0;
                    else
                        tar[i*SizeJ+j] = 255;
                    double err = temp[i*SizeJ+j] - tar[i*SizeJ+j];
                    for(int m=-2; m<=2; m++)
                    {
                        for(int n=-2; n<=2; n++)
                        {
                            if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                            {
                                temp[(i+m)*SizeJ+j+n] += err*DiffusionM0[m+2][n+2];
                            }
                        }
                    }
                }

            }
            else
            {
                for(int j=SizeJ-1; j>=0; j--)
                {
                    if(temp[i*SizeJ+j]<Thres)
                        tar[i*SizeJ+j] = 0;
                    else
                        tar[i*SizeJ+j] = 255;
                    double err = temp[i*SizeJ+j] - tar[i*SizeJ+j];
                    for(int m=-2; m<=2; m++)
                    {
                        for(int n=-2; n<=2; n++)
                        {
                            if(i+m>=0&&i+m<SizeI&&j+n>=0&&j+n<SizeJ)
                            {
                                temp[(i+m)*SizeJ+j+n] += err*DiffusionM1[m+2][n+2];
                            }
                        }
                    }
                }
            }
        }
        break;
    }
    }

    delete temp;
}

void InverseHalftoning(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ)
{
    // Inverse Halftoning with Cascade Algorithm
    int Margin = 2;
    int ExpandedI = SizeI+2*Margin;
    int ExpandedJ = SizeJ+2*Margin;
    unsigned char* expandedImage = new unsigned char[ExpandedI*ExpandedJ];
    unsigned char* buff1 = new unsigned char[SizeI*SizeJ];
    // First stage
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            double ABRL = 0;
            int rl = 0, cl = 0;

            // for row
            int r = i;
            while(r>=0 && src[r*SizeJ+j]==src[i*SizeJ+j])
            {
                rl ++;
                r --;
            }
            r = i+1;
            while(r<SizeI && src[r*SizeJ+j]==src[i*SizeJ+j])
            {
                rl ++;
                r ++;
            }
            // for column
            int c = j;
            while(c>=0 && src[i*SizeJ+c]==src[i*SizeJ+j])
            {
                cl ++;
                c --;
            }
            c = j+1;
            while(c<SizeJ && src[i*SizeJ+c]==src[i*SizeJ+j])
            {
                cl ++;
                c ++;
            }

            if(src[i*SizeJ+j]==0)
            {
                ABRL = 0.5*(1.0/(cl+1) + 1.0/(rl+1));
            }
            else
            {
                ABRL = 0.5*(1.0*cl/(cl+1) + 1.0*rl/(rl+1));
            }

            /*
                        int T = 5;
                        for(int m=-T;m<=T;m++)
                        {
                            if(i+m>=0 && i+m<SizeI)
                            {
                                if(src[(i+m)*SizeJ+j]==255)
                                {
                                    rl ++;
                                }
                            }
                            if(j+m>=0 && j+m<SizeJ)
                            {
                                if(src[i*SizeJ+j+m]==255)
                                {
                                    cl ++;
                                }
                            }
                        }
                        ABRL = 0.5*(1.0*cl/(2*T+1) + 1.0*rl/(2*T+1));
            */
            buff1[i*SizeJ+j] = 255*ABRL;
        }
    }
    ImageExpandRepeat(Margin,buff1,expandedImage,SizeI,SizeJ,1);

    double MAXVAR = 50;
    double beta = 50;
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            double mu = 0;
            double sigma_sq = 0;
            double K = 800;
            // Calculate mu
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    mu += expandedImage[ExpandedJ*(i+Margin+m)+j+Margin+n];
                }
            }
            mu = mu / (Margin*2+1) / (Margin*2+1);
            // Calculate sigma
            for(int m=-Margin; m<=Margin; m++)
            {
                for(int n=-Margin; n<=Margin; n++)
                {
                    sigma_sq += (expandedImage[ExpandedJ*(i+Margin+m)+j+Margin+n]-mu)*(expandedImage[ExpandedJ*(i+Margin+m)+j+Margin+n]-mu);
                }
            }
            sigma_sq = sigma_sq/24;
            double tmp = mu+(sigma_sq/(sigma_sq+K))*(buff1[i*SizeJ+j]-mu);

            tar[i*SizeJ+j] = tmp;
            if(sigma_sq <= MAXVAR)
            {
                if(tmp+beta<mu)
                {
                    tar[i*SizeJ+j] = tmp+beta;
                }
                if(tmp-beta>mu)
                {
                    tar[i*SizeJ+j] = tmp-beta;
                }
            }
        }
    }

    delete buff1;
    delete expandedImage;
}

bool MatrixMul(float* M1,int m1,int n1, float* M2, int m2, int n2, float* R)
{
    if(n1!=m2)
    {
        return false;
    }
    for(int i=0; i<m1; i++)
    {
        for(int j=0; j<n2; j++)
        {
            R[i*n2+j] = 0;
            for(int k=0; k<n1; k++)
            {
                R[i*n2+j] += M1[i*n1+k]*M2[k*n2+j];
            }
        }
    }
    return true;
}

float GetDistance(float* P1,float* P2,int len)
{
    float r = 0;
    for(int i=0; i<len; i++)
    {
        r += (P1[i]-P2[i])*(P1[i]-P2[i]);
    }
    return sqrt(r);
}

double GetDistanceD(double* P1,double* P2,int len)
{
    double r = 0;
    for(int i=0; i<len; i++)
    {
        r += (P1[i]-P2[i])*(P1[i]-P2[i]);
    }
    return sqrt(r);
}
