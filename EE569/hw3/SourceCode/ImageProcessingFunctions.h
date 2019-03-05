/**
    Name:	Manan Vyas
	USC ID:	7483-8632-00
	email:	mvyas@usc.edu
*/

#ifndef __IMAGEPROCESSINGFUNCTIONS_H__
#define __IMAGEPROCESSINGFUNCTIONS_H__

typedef struct RGB{
    unsigned char R;
    unsigned char G;
    unsigned char B;
}RGB;

void HistogramEqualization(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ);

void GrayCopy(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ);

void Color2Gray(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ);

void Gray2Binary(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int Threshold,char type);

void RGB2YUV(const unsigned char *src, unsigned char *Y, int *U, int *V,int SizeI, int SizeJ);

void YUV2RGB(const unsigned char *Y, const int *U, const int *V,unsigned char *tar, int SizeI, int SizeJ);

void SplitRGB(const unsigned char *src, unsigned char *R,unsigned char *G,unsigned char *B,int SizeI, int SizeJ);

void GrayScaleInverse(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ);

void EdgeDetectionGradient(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int OperatorType,int thres);

void EdgeDetectionLaplacian(const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int thres);

void EdgeDetectionDOG(double sigma, double k, int FilterSize, const unsigned char *src, unsigned char *tar,int SizeI, int SizeJ,int thres);

int MorphOperation(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ,char type);

void Thinning(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ);

void Dilation4(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ);

void Erode4(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ);

void FillSingleHole(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ);

int CountConnectedObjects(const unsigned char *src, int SizeI, int SizeJ, int Thres);

int CountSingleDots(const unsigned char *src, int SizeI, int SizeJ);

int CountPattern(const unsigned char *src, int SizeI, int SizeJ, int pattern,int patternSize);

void Dithering(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ,char type,int *indexM,int N);

void ErrorDiffusion(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ,int Method,int Thres);

void InverseHalftoning(const unsigned char *src, unsigned char *tar, int SizeI, int SizeJ);

bool MatrixMul(float* M1,int m1,int n1, float* M2, int m2, int n2, float* R);

float GetDistance(float* P1,float* P2,int len);

double GetDistanceD(double* P1,double* P2,int len);

#endif
