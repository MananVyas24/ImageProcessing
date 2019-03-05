/**
    Name:	Manan Vyas
	USC ID:	7483-8632-00
	email:	mvyas@usc.edu
*/

#include "IO.h"

#include <stdio.h>
#include <string.h>

void ReadImageFromFile(const char* filename,unsigned char *ImageData,int SizeI,int SizeJ,int BytesPerPixel)
{
    FILE *file;
    if (!(file=fopen(filename,"rb")))
    {
        printf("Cannot open file: %s\n",filename);
        return;
    }
    fread(ImageData, sizeof(unsigned char), SizeI*SizeJ*BytesPerPixel, file);
    fclose(file);
}

void WriteImageToFile(const char* filename,const unsigned char *ImageData,int SizeI,int SizeJ,int BytesPerPixel)
{
    FILE *file;
    if (!(file=fopen(filename,"wb")))
    {
        printf("Cannot open file: %s\n",filename);
        return;
    }
    fwrite(ImageData, sizeof(unsigned char), SizeI*SizeJ*BytesPerPixel, file);
    fclose(file);
}

void GetHistogram(const unsigned char *ImageData,int SizeI,int SizeJ, int *histogram)
{
    memset(histogram,0,256*sizeof(int)); // Initialize it into 0
    for(int i=0; i<SizeI; i++)
    {
        for(int j=0; j<SizeJ; j++)
        {
            histogram[ImageData[i*SizeJ+j]] ++;
        }
    }
}

void WriteHistogramToFile(const char* filename, const int *histogram)
{
    FILE *file;
    if (!(file=fopen(filename,"wb")))
    {
        printf("Cannot open file: %s\n",filename);
        return;
    }
    for(int i=0; i<256; i++)
    {
        fprintf(file,"%d ",histogram[i]);
    }
    fclose(file);
}

void WriteTransferToFile(const char* filename, const int *Transfer)
{
    FILE *file;
    if (!(file=fopen(filename,"wb")))
    {
        printf("Cannot open file: %s\n",filename);
        return;
    }
    for(int i=0; i<256; i++)
    {
        fprintf(file,"%d ",Transfer[i]);
    }
    fclose(file);
}

void WriteDataToFile(const char* filename, const int *data, int length)
{
    FILE *file;
    if (!(file=fopen(filename,"wb")))
    {
        printf("Cannot open file: %s\n",filename);
        return;
    }
    fprintf(file,"%d ",length);
    for(int i=0; i<length; i++)
    {
        fprintf(file,"%d ",data[i]);
    }
    fclose(file);
}
